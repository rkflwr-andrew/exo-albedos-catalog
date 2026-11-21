#!/usr/bin/env python3
"""
Reads data/planets.csv (one row per planet) and data/measurements.csv (one row per
eclipse measurement), runs the albedo calculations, and writes site/index.html and
site/results.csv. Measurements are grouped by planet; the planet cell is merged.
"""
import csv, html, os, sys
from albedo_core import compute_all, load_filter

_filters = {}


def read(path):
    with open(path, newline="", encoding="utf-8-sig") as f:   # utf-8-sig: tolerate Excel BOM
        return [{k.strip(): (v or "").strip() for k, v in r.items() if k} for r in csv.DictReader(f)]


def num(row, key, default=None):
    v = row.get(key, "").strip()
    return float(v) if v else default


def link(text, url):
    text = html.escape(text or "")
    url = (url or "").strip()
    if not url:
        return text
    if url.lower().startswith("10."):          # bare DOI
        url = "https://doi.org/" + url
    if not url.lower().startswith(("http://", "https://")):
        return text
    return f'<a href="{html.escape(url, quote=True)}" target="_blank" rel="noopener">{text or "link"}</a>'


def pm(v):
    return f'{v[0]:.3f}<sup>+{v[1]:.3f}</sup><sub>&minus;{v[2]:.3f}</sub>'


def sort_key(m):
    try:
        return float(m.get("mean_wl_nm", ""))
    except ValueError:
        return float("inf")


OFF = ("no", "n", "false", "0", "off")


def main():
    planets = read("data/planets.csv")
    meas = read("data/measurements.csv")
    names = [p.get("planet", "") for p in planets]
    if "" in names:
        sys.exit("ERROR: a row in data/planets.csv has no planet name")
    if len(set(names)) != len(names):
        sys.exit("ERROR: duplicate planet names in data/planets.csv")
    by_planet = {n: [] for n in names}
    for i, m in enumerate(meas, start=2):
        if m.get("planet") not in by_planet:
            sys.exit(f"ERROR in data/measurements.csv line {i}: planet '{m.get('planet')}' "
                     f"is not in data/planets.csv (names must match exactly)")
        if any(m.get(k) for k in ("td_k", "td_err_k", "td_ref", "td_ref_url")):
            sys.exit(f"ERROR in data/measurements.csv line {i}: temperatures are set once per planet; "
                     f"move td_k, td_err_k, td_ref, td_ref_url to data/planets.csv")
        m["_line"] = i
        by_planet[m["planet"]].append(m)

    groups, out, notes = [], [], []
    for p in planets:
        ms = sorted(by_planet[p["planet"]], key=sort_key)       # stable: ties keep CSV order
        if not ms:
            continue
        rp_rs, a_rs, ts = num(p, "rp_rs"), num(p, "a_rs"), num(p, "ts")
        ts_err = num(p, "ts_err", 0.0)
        td, td_err = num(p, "td_k"), num(p, "td_err_k", 0.0)     # one Td per planet
        missing_pp = [k for k, v in (("rp_rs", rp_rs), ("a_rs", a_rs), ("ts", ts)) if v is None]
        rows = []
        for k, m in enumerate(ms):
            # ---- decide whether this row gets a calculation, and from which input ----
            used = next((kind for kind in ("depth", "amp")
                         if m.get(f"{kind}_ppm") and m.get(f"{kind}_err_ppm")), None)  # depth has priority
            reason, deliberate = None, False
            if m.get("calculate", "").lower() in OFF:
                reason, deliberate = "calculation turned off", True
            elif used is None:
                reason = ("no error bar on the depth/amplitude" if (m.get("depth_ppm") or m.get("amp_ppm"))
                          else "no eclipse depth or amplitude")
            elif not m.get("filter"):
                reason = "no filter bandpass"
            elif missing_pp:
                reason = "missing planet parameter(s): " + ", ".join(missing_pp)
            calc = reason is None
            if not calc and not deliberate:
                notes.append(f"line {m['_line']} ({p['planet']}, {m.get('bandpass') or '?'}): not calculated - {reason}")

            res = {}
            if calc:
                try:
                    fk = (m["filter"], m.get("wl_unit") or "m")
                    if fk not in _filters:
                        _filters[fk] = load_filter(*fk)
                    wl, tr = _filters[fk]
                    photon = m.get("photon", "").lower() in ("1", "true", "yes")
                    res = compute_all(num(m, f"{used}_ppm"), num(m, f"{used}_err_ppm"), rp_rs, a_rs,
                                      ts, ts_err, wl, tr, photon, td, td_err)
                except Exception as e:
                    sys.exit(f"ERROR in data/measurements.csv line {m['_line']} ({p['planet']}): {e}")

            def cell(v):
                return f'<td>{pm(v)}</td>'

            def pair(a, b, flag):
                dag = "&dagger;" if flag else ""
                return f'<td>inst: {pm(a)}<br>hom: {pm(b)}{dag}</td>'

            def inp(kind):
                v = m.get(f"{kind}_ppm")
                if not v:
                    return '<td>&mdash;</td>'
                e = m.get(f"{kind}_err_ppm")
                mark = '<br><small>(used for albedo)</small>' if (calc and used == kind) else ''
                return f'<td>{html.escape(v)}' + (f' &plusmn; {html.escape(e)}' if e else '') + f'{mark}</td>'

            if calc:
                calc_cells = (cell(res["no_thermal"])
                              + (cell(res["measured_td"]) if td is not None else '<td>&mdash;</td>')
                              + pair(res["self_inst"], res["self_hom"], res["capped_inst"] or res["capped_hom"])
                              + pair(res["zero_inst"], res["zero_hom"], False))
            else:
                calc_cells = f'<td colspan="4" class="na">Not calculated<br><small>{html.escape(reason)}</small></td>'

            planet_cell = ""
            if k == 0:
                parts = []
                if rp_rs is not None:
                    parts.append(f'R<sub>p</sub>/R<sub>s</sub>={rp_rs:g}')
                if a_rs is not None:
                    parts.append(f'a/R<sub>s</sub>={a_rs:g}')
                if ts is not None:
                    parts.append(f'T<sub>s</sub>={ts:g}' + (f'&plusmn;{ts_err:g}' if ts_err else '') + ' K')
                pref = link(p.get("ref"), p.get("ref_url"))
                small = ", ".join(parts) + (("<br>" if parts else "") + pref if pref else "")
                planet_cell = (f'<td rowspan="{len(ms)}" class="planet"><b>{html.escape(p["planet"])}</b>'
                               + (f'<br><small>{small}</small>' if small else '') + '</td>')
                if td is not None:
                    planet_cell += (f'<td rowspan="{len(ms)}" class="planet">{td:g}'
                                    + (f' &plusmn; {td_err:g}' if td_err else '') + ' K'
                                    + f'<br><small>{link(p.get("td_ref"), p.get("td_ref_url"))}</small></td>')
                else:
                    planet_cell += f'<td rowspan="{len(ms)}" class="planet">&mdash;</td>'
            lit = html.escape(m.get("lit_albedo", "")) or "&mdash;"
            rows.append(
                "<tr>" + planet_cell
                + f'<td>{html.escape(m.get("bandpass", "")) or "&mdash;"}</td>'
                + f'<td>{html.escape(m.get("mean_wl_nm", "")) or "&mdash;"}</td>'
                + inp("depth") + inp("amp")
                + f'<td>{lit}<br><small>{link(m.get("ref"), m.get("ref_url"))}</small></td>'
                + calc_cells + "</tr>")
            flat = {"planet": p["planet"], "bandpass": m.get("bandpass", ""), "mean_wl_nm": m.get("mean_wl_nm", ""),
                    "lit_albedo": m.get("lit_albedo", ""), "ref": m.get("ref", ""),
                    "input_used": used if calc else "",
                    "td_k": td if td is not None else "", "td_err_k": td_err if td is not None else ""}
            for key, v in res.items():
                if isinstance(v, tuple):
                    flat[key], flat[key + "_up"], flat[key + "_dn"] = v
            out.append(flat)
        groups.append("<tbody>\n" + "\n".join(rows) + "\n</tbody>")

    head = ["Planet", "Literature T<sub>d</sub> for thermal correction / reference",
            "Bandpass", "Mean &lambda; (nm)", "Eclipse depth (ppm)", "Phase curve amplitude (ppm)",
            "Literature A<sub>g</sub> / reference",
            "A<sub>g</sub>: no thermal correction", "A<sub>g</sub>: thermal-corrected (literature T<sub>d</sub>)",
            "A<sub>g</sub>: self-consistent (A<sub>B</sub>=1.5A<sub>g</sub>)", "A<sub>g</sub>: A<sub>B</sub>=0"]
    os.makedirs("site", exist_ok=True)
    with open("site/index.html", "w", encoding="utf-8") as f:
        f.write(PAGE.replace("@@HEAD@@", "".join(f"<th>{h}</th>" for h in head))
                    .replace("@@BODY@@", "\n".join(groups)))
    keys = []
    for d in out:
        keys += [k for k in d if k not in keys]
    with open("site/results.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(out)
    ncalc = sum(1 for d in out if d["input_used"])
    print(f"Built {len(out)} measurements for {len(groups)} planets ({ncalc} with calculations) -> site/index.html")
    for n in notes:
        print("NOTE:", n)


PAGE = """<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Exoplanet geometric albedos from secondary eclipses</title>
<style>
:root{--bg:#fff;--fg:#1a1a1a;--line:#d0d0d0;--head:#eee;--link:#0b5cad}
@media(prefers-color-scheme:dark){:root{--bg:#161616;--fg:#e8e8e8;--line:#444;--head:#262626;--link:#7ab7ff}}
body{font-family:system-ui,sans-serif;margin:1.5rem;background:var(--bg);color:var(--fg)}
a{color:var(--link)} .wrap{overflow-x:auto}
table{border-collapse:collapse;font-size:.9rem} td,th{border:1px solid var(--line);padding:.4rem .6rem;text-align:left;vertical-align:top}
th{background:var(--head);position:sticky;top:0} sup,sub{font-size:.7em;line-height:0}
tbody{border-top:3px solid var(--line)} td.na{opacity:.7;font-style:italic} td.planet{background:var(--head)}
input{padding:.4rem;margin-bottom:.8rem;width:18rem} small{opacity:.85}
</style>
<h1>Eclipse depths and geometric albedos</h1>
<input id="q" placeholder="Filter by planet, bandpass, reference..." aria-label="filter table">
<div class="wrap"><table id="t"><thead><tr>@@HEAD@@</tr></thead>
@@BODY@@
</table></div>
<p><small>Values are median with +/&minus; 1&sigma; (one-sided shifts of each input added in quadrature).
inst = instant re-radiation (f = 2/3), hom = homogeneous re-radiation (f = 1/4).
&dagger; = no self-consistent solution with A<sub>B</sub> &le; 1 (A<sub>g</sub> &gt; 2/3); use the A<sub>B</sub>=0 column.
Rows without a depth/amplitude, filter, or planet parameters (or with calculation turned off) are shown without a calculation. Full numbers: <a href="results.csv">results.csv</a>.</small></p>
<script>
document.getElementById('q').addEventListener('input',e=>{const s=e.target.value.toLowerCase();
 for(const b of document.querySelectorAll('#t tbody'))b.style.display=b.textContent.toLowerCase().includes(s)?'':'none';});
</script>
</html>
"""

if __name__ == "__main__":
    main()
