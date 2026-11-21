#!/usr/bin/env python3
"""Command-line geometric albedo calculator (see albedo_core.py for the physics).

Example:
  python albedo.py --depth 60 --err 15 --rp-rs 0.1 --a-rs 4.5 --ts 6000 --ts-err 100 \
      --td 2500 --td-err 150 --filter filters/example_filter.txt --wl-unit m --photon
"""
import argparse
from albedo_core import compute_all, load_filter


def show(label, v):
    print(f"  {label:<44s} {v[0]:8.4f}  +{v[1]:.4f} / -{v[2]:.4f}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--depth", type=float, required=True, help="eclipse depth (ppm)")
    p.add_argument("--err", type=float, required=True, help="1-sigma depth error (ppm)")
    p.add_argument("--rp-rs", type=float, required=True)
    p.add_argument("--a-rs", type=float, required=True)
    p.add_argument("--ts", type=float, required=True, help="stellar Teff (K)")
    p.add_argument("--ts-err", type=float, default=0.0)
    p.add_argument("--td", type=float, default=None, help="measured dayside brightness temperature (K)")
    p.add_argument("--td-err", type=float, default=0.0)
    p.add_argument("--filter", required=True, help="text file: wavelength, transmission")
    p.add_argument("--wl-unit", choices=["m", "um", "nm", "angstrom"], default="m")
    p.add_argument("--photon", action="store_true", help="filter is a photon-counting response")
    a = p.parse_args()

    wl, tr = load_filter(a.filter, a.wl_unit)
    r = compute_all(a.depth, a.err, a.rp_rs, a.a_rs, a.ts, a.ts_err, wl, tr,
                    a.photon, a.td, a.td_err)
    print("1) No thermal correction (upper limit):")
    show("Ag", r["no_thermal"])
    if "measured_td" in r:
        print("2) Thermal-corrected with measured Td:")
        show("Ag", r["measured_td"])
    print("3) Self-consistent (A_B = 1.5 Ag; A_B = 0 if Ag <= 0):")
    show("Ag, instant re-radiation (f=2/3)", r["self_inst"])
    show("Ag, homogeneous re-radiation (f=1/4)", r["self_hom"])
    if r["capped_inst"] or r["capped_hom"]:
        print("  !! no solution with A_B <= 1 (Ag > 2/3): use calculation 4")
    print("4) Bond albedo = 0:")
    show("Ag, instant re-radiation (f=2/3)", r["zero_inst"])
    show("Ag, homogeneous re-radiation (f=1/4)", r["zero_hom"])


if __name__ == "__main__":
    main()
