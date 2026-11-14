"""Extract the 13 x 3 truncation-scan summary (w0, wa, sigma) used in the paper
figure and tables, from the hardcoded arrays in plot_final_fig.py, and write
results/cpl_scan_summary.csv."""
import csv
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
PLOT = os.path.join(HERE, "plot_final_fig.py")
OUT = os.path.join(HERE, "..", "results", "cpl_scan_summary.csv")

src = open(PLOT, encoding="utf8").read()
pat = re.compile(r"^(\w+)\s*=\s*np\.array\(\[([^\]]*)\]\)", re.M)

arrays = {}
for name, body in pat.findall(src):
    arrays[name] = [float(x) for x in body.split(",") if x.strip()]

zcuts = arrays["zcuts"]
rows = []
for i, zc in enumerate(zcuts):
    for ds, w0, w0e, wa, wae, sig in [
        ("Pantheon+", "pp_w0", "pp_w0e", "pp_wa", "pp_wae", "pp_sig"),
        ("DESY5", "dy5_w0", "dy5_w0e", "dy5_wa", "dy5_wae", "dy5_sig"),
        ("DES-Dovekie", "dov_w0", "dov_w0e", "dov_wa", "dov_wae", "dov_sig"),
    ]:
        rows.append([f"{zc:.2f}", ds, arrays[w0][i], arrays[w0e][i],
                     arrays[wa][i], arrays[wae][i], arrays[sig][i]])

with open(OUT, "w", newline="", encoding="utf8") as f:
    w = csv.writer(f)
    w.writerow(["z_cut", "dataset", "w0", "w0_err", "wa", "wa_err", "sigma"])
    w.writerows(rows)
print(f"wrote {len(rows)} rows -> {os.path.normpath(OUT)}")
