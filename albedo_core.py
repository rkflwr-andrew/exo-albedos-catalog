"""
Core calculations for geometric albedo from secondary-eclipse depth.

Equations
---------
  d_refl = Ag (Rp/a)^2                        ->  Ag = d (a/Rp)^2
  Ag = d (a/Rp)^2 - [B(Td)/B(Ts)] (a/Rs)^2    (thermal-corrected; B integrated over the bandpass)
  Td = Ts sqrt(Rs/a) [f (1 - A_B)]^(1/4)      f = 2/3 instant re-radiation, 1/4 homogeneous

Errors are propagated deterministically (no Monte Carlo): each uncertain input
is shifted by +/- 1 sigma, the result is recomputed, and the shifts are added in
quadrature separately for the upper and lower side.  This is a few dozen function
evaluations per measurement, so hundreds of rows take seconds.
"""
import numpy as np

H = 6.62607015e-34
C = 2.99792458e8
KB = 1.380649e-23
trapz = getattr(np, "trapezoid", None) or np.trapz
WL_SCALE = {"m": 1.0, "um": 1e-6, "nm": 1e-9, "angstrom": 1e-10}


def load_filter(path, wl_unit="m"):
    d = np.loadtxt(path, comments=("#", "%"))
    d = d[np.argsort(d[:, 0])]
    return d[:, 0] * WL_SCALE[wl_unit], d[:, 1]


def planck(wl, T):
    """Planck B_lambda. wl: (n,1) metres, T: (1,m) K -> (n,m)."""
    T = np.maximum(T, 1.0)
    x = H * C / (wl * KB * T)
    with np.errstate(over="ignore"):
        return 2 * H * C**2 / wl**5 / np.expm1(x)


def band_ratio(Td, Ts, wl, tr, photon=False):
    """Bandpass-integrated B(Td)/B(Ts), vectorised over equal-length Td, Ts arrays."""
    w = (wl * tr if photon else tr)[:, None]
    wl2 = wl[:, None]
    Td = np.atleast_1d(np.asarray(Td, float))[None, :]
    Ts = np.atleast_1d(np.asarray(Ts, float))[None, :]
    num = trapz(planck(wl2, Td) * w, wl, axis=0)
    den = trapz(planck(wl2, Ts) * w, wl, axis=0)
    return num / den


def Td_irr(Ts, a_rs, f, A_B):
    return Ts * a_rs**-0.5 * (f * np.clip(1 - A_B, 0, None)) ** 0.25


def zero_bond(X, Ts, a_rs, f, wl, tr, photon):
    """Ag with Bond albedo = 0. Negative values allowed, no upper limit."""
    Td = Td_irr(Ts, a_rs, f, 0.0)
    return X - band_ratio(Td, Ts, wl, tr, photon) * a_rs**2, Td


def solve_self_consistent(X, Ts, a_rs, f, wl, tr, photon, n_iter=40):
    """
    Self-consistent Ag with A_B = 1.5 Ag. For Ag <= 0, A_B is set to 0.
    Otherwise bisection on [0, 2/3] (A_B <= 1).
    Returns (Ag, Td, capped); capped marks points with no solution below Ag = 2/3.
    """
    A0, Td0 = zero_bond(X, Ts, a_rs, f, wl, tr, photon)
    neg = A0 <= 0
    lo = np.zeros_like(X)
    hi = np.full_like(X, 2.0 / 3.0)

    def g(A):
        Td = Td_irr(Ts, a_rs, f, 1.5 * A)
        return A - X + band_ratio(Td, Ts, wl, tr, photon) * a_rs**2

    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        pos = g(mid) > 0
        hi = np.where(pos, mid, hi)
        lo = np.where(pos, lo, mid)
    A = 0.5 * (lo + hi)
    capped = (X > 2.0 / 3.0) & ~neg
    A = np.where(neg, A0, A)
    Td = np.where(neg, Td0, Td_irr(Ts, a_rs, f, 1.5 * A))
    return A, Td, capped


def propagate(F, vals, sigs):
    """
    Deterministic error propagation. F(*arrays) -> array, evaluated on a batch
    holding the central point and each input shifted by +/- sigma.
    Returns (central, upper_err, lower_err), quadrature sums of the one-sided shifts.
    """
    k = len(vals)
    rows = [list(vals)]
    for i in range(k):
        for s in (+1, -1):
            r = list(vals)
            r[i] += s * sigs[i]
            rows.append(r)
    P = np.array(rows, float).T
    out = np.asarray(F(*P), float)
    f0 = out[0]
    up2 = dn2 = 0.0
    for i in range(k):
        sp, sm = out[1 + 2 * i] - f0, out[2 + 2 * i] - f0
        up2 += max(sp, sm, 0.0) ** 2
        dn2 += min(sp, sm, 0.0) ** 2
    return float(f0), float(np.sqrt(up2)), float(np.sqrt(dn2))


def compute_all(depth_ppm, err_ppm, rp_rs, a_rs, ts, ts_err, wl, tr,
                photon=False, td=None, td_err=0.0):
    """
    Returns a dict of case -> (value, +err, -err), plus 'capped_*' flags.
    Cases: no_thermal, measured_td (if td given), self_inst, self_hom,
           zero_inst, zero_hom.
    """
    a_rp = a_rs / rp_rs
    res = {}
    Xf = lambda d: d * 1e-6 * a_rp**2

    res["no_thermal"] = propagate(lambda d: Xf(d), [depth_ppm], [err_ppm])

    if td is not None:
        res["measured_td"] = propagate(
            lambda d, T, Tdd: Xf(d) - band_ratio(Tdd, T, wl, tr, photon) * a_rs**2,
            [depth_ppm, ts, td], [err_ppm, ts_err, td_err])

    for name, f in (("inst", 2 / 3), ("hom", 0.25)):
        res["zero_" + name] = propagate(
            lambda d, T, f=f: zero_bond(Xf(d), T, a_rs, f, wl, tr, photon)[0],
            [depth_ppm, ts], [err_ppm, ts_err])

        flag = {}

        def F(d, T, f=f, flag=flag):
            A, _, capped = solve_self_consistent(Xf(d), T, a_rs, f, wl, tr, photon)
            flag["c"] = bool(capped[0])
            return A

        res["self_" + name] = propagate(F, [depth_ppm, ts], [err_ppm, ts_err])
        res["capped_" + name] = flag.get("c", False)
    return res
