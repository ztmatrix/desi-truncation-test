"""Final paper figure: Three-way CPL truncation scan with all 39 full-Planck chains.

Values from analyze_scan.py on pp_full_planck/, dy5_full_planck/, dovekie_full_planck/
(server chains) with chains/fp_dy5_z0p00.* (local 172k chain) for DY5 z0p00.
"""
import glob
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

zcuts = np.array([0.00, 0.02, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.14, 0.15, 0.20])

# ── Pantheon+ CPL (full Planck) ──
pp_w0  = np.array([-0.8613, -0.8703, -0.8404, -0.8321, -0.8143, -0.7973, -0.7888, -0.7938, -0.8067, -0.7959, -0.7764, -0.7841, -0.6555])
pp_w0e = np.array([ 0.0557,  0.0560,  0.0733,  0.0778,  0.0835,  0.0911,  0.0929,  0.0953,  0.0955,  0.0983,  0.1059,  0.1095,  0.1391])
pp_wa  = np.array([-0.4483, -0.4205, -0.4954, -0.5191, -0.5770, -0.6197, -0.6389, -0.6327, -0.5922, -0.6236, -0.6755, -0.6533, -1.0298])
pp_wae = np.array([ 0.2165,  0.2161,  0.2553,  0.2683,  0.2834,  0.3044,  0.3048,  0.3153,  0.3119,  0.3199,  0.3435,  0.3483,  0.4470])
pp_sig = np.array([2.51, 2.32, 2.18, 2.16, 2.22, 2.23, 2.27, 2.16, 2.02, 2.08, 2.11, 1.97, 2.49])

# ── DESY5 CPL (full Planck; z0p00 from local 172k chain) ──
dy5_w0  = np.array([-0.7766, -0.7747, -0.8136, -0.8427, -0.8267, -0.8055, -0.8253, -0.8330, -0.8505, -0.7861, -0.7722, -0.8029, -0.7749])
dy5_w0e = np.array([ 0.0582,  0.0596,  0.0736,  0.0792,  0.0879,  0.1016,  0.1086,  0.1042,  0.1105,  0.1101,  0.1168,  0.1097,  0.1201])
dy5_wa  = np.array([-0.6949, -0.7032, -0.6018, -0.5096, -0.5588, -0.6050, -0.5611, -0.5365, -0.4825, -0.6668, -0.7044, -0.6138, -0.6928])
dy5_wae = np.array([ 0.2365,  0.2385,  0.2681,  0.2802,  0.3017,  0.3310,  0.3514,  0.3337,  0.3549,  0.3554,  0.3814,  0.3501,  0.3804])
dy5_sig = np.array([4.03, 3.96, 2.54, 1.99, 1.97, 1.91, 1.62, 1.63, 1.37, 1.94, 1.95, 1.80, 1.88])

# ── DES-Dovekie CPL (full Planck) ──
dov_w0  = np.array([-0.8344, -0.8377, -0.8401, -0.8654, -0.8700, -0.8701, -0.8796, -0.8862, -0.9003, -0.8708, -0.8586, -0.8877, -0.8274])
dov_w0e = np.array([ 0.0571,  0.0570,  0.0688,  0.0755,  0.0847,  0.0923,  0.0934,  0.0956,  0.0961,  0.0997,  0.0997,  0.1018,  0.1093])
dov_wa  = np.array([-0.5609, -0.5543, -0.5299, -0.4547, -0.4453, -0.4434, -0.4197, -0.3976, -0.3578, -0.4381, -0.4716, -0.3969, -0.5502])
dov_wae = np.array([ 0.2318,  0.2330,  0.2612,  0.2700,  0.2993,  0.3084,  0.3151,  0.3148,  0.3170,  0.3279,  0.3328,  0.3333,  0.3508])
dov_sig = np.array([2.94, 2.87, 2.34, 1.78, 1.54, 1.45, 1.34, 1.26, 1.13, 1.34, 1.44, 1.19, 1.59])

# ═══════════════════════════════════════════════════════════
# Figure 1: Four-panel comparison
# ═══════════════════════════════════════════════════════════
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SN_DATA = {
    "Pantheon+": os.path.join(BASE_DIR, "..", "packages", "data", "sn_data", "PantheonPlus", "Pantheon+SH0ES.dat"),
    "DESY5": os.path.join(BASE_DIR, "..", "packages", "data", "sn_data", "DESY5", "DES-SN5YR_HD.csv"),
    "DES-Dovekie": os.path.join(BASE_DIR, "..", "packages", "data", "sn_data", "DES-Dovekie", "DES-Dovekie_HD.csv"),
}
fig, axes = plt.subplots(2, 2, figsize=(16, 13))

colors = {"Pantheon+": "#2166ac", "DESY5": "#b2182b", "DES-Dovekie": "#4daf4a"}
markers = {"Pantheon+": "s", "DESY5": "o", "DES-Dovekie": "D"}

# ── Panel (a): w0 vs z_cut ──
ax = axes[0, 0]
for label, w0, w0e, c, m in [
    ("Pantheon+", pp_w0, pp_w0e, colors["Pantheon+"], markers["Pantheon+"]),
    ("DESY5", dy5_w0, dy5_w0e, colors["DESY5"], markers["DESY5"]),
    ("DES-Dovekie", dov_w0, dov_w0e, colors["DES-Dovekie"], markers["DES-Dovekie"]),
]:
    ax.errorbar(zcuts, w0, yerr=w0e, color=c, marker=m, linewidth=2.2, markersize=9, capsize=4, capthick=1.5, label=label)
ax.axhline(y=-1.0, color="gray", linestyle="--", linewidth=1.5, alpha=0.6)
# Shade Dovekie transition zone
ax.axvspan(0.05, 0.08, alpha=0.06, color=colors["DES-Dovekie"])
ax.annotate(r"$\leftarrow$ Dovekie $\rightarrow \Lambda$CDM", xy=(0.065, -0.99), fontsize=10,
            color=colors["DES-Dovekie"], fontstyle="italic", fontweight="bold")
ax.set_xlabel(r"Redshift cut $z_{\rm cut}$", fontsize=13)
ax.set_ylabel(r"$w_0$", fontsize=15)
ax.set_title(r"(a)  $w_0$ stability under truncation", fontsize=14, fontweight="bold", loc="left")
ax.legend(fontsize=11, loc="lower left", framealpha=0.9)
ax.set_xlim(-0.005, 0.205)
ax.set_ylim(-1.02, -0.50)
ax.grid(True, alpha=0.2)

# ── Panel (b): wa vs z_cut ──
ax = axes[0, 1]
for label, wa, wae, c, m in [
    ("Pantheon+", pp_wa, pp_wae, colors["Pantheon+"], markers["Pantheon+"]),
    ("DESY5", dy5_wa, dy5_wae, colors["DESY5"], markers["DESY5"]),
    ("DES-Dovekie", dov_wa, dov_wae, colors["DES-Dovekie"], markers["DES-Dovekie"]),
]:
    ax.errorbar(zcuts, wa, yerr=wae, color=c, marker=m, linewidth=2.2, markersize=9, capsize=4, capthick=1.5, label=label)
ax.axhline(y=0.0, color="gray", linestyle="--", linewidth=1.5, alpha=0.6)
ax.set_xlabel(r"Redshift cut $z_{\rm cut}$", fontsize=13)
ax.set_ylabel(r"$w_a$", fontsize=15)
ax.set_title(r"(b)  $w_a$ stability under truncation", fontsize=14, fontweight="bold", loc="left")
ax.legend(fontsize=11, loc="lower left", framealpha=0.9)
ax.set_xlim(-0.005, 0.205)
ax.set_ylim(-1.6, 0.25)
ax.grid(True, alpha=0.2)

# ── Panel (c): 2D significance ──
ax = axes[1, 0]
ax.plot(zcuts, pp_sig, color=colors["Pantheon+"], marker=markers["Pantheon+"], linewidth=2.2, markersize=9, label="Pantheon+")
ax.plot(zcuts, dy5_sig, color=colors["DESY5"], marker=markers["DESY5"], linewidth=2.2, markersize=9, label="DESY5")
ax.plot(zcuts, dov_sig, color=colors["DES-Dovekie"], marker=markers["DES-Dovekie"], linewidth=2.5, markersize=10, label="DES-Dovekie")
ax.axhline(y=1.0, color="gray", linestyle=":", alpha=0.4)
ax.axhline(y=2.0, color="gray", linestyle="--", alpha=0.4)
ax.axhline(y=3.0, color="gray", linestyle="-.", alpha=0.4)
ax.text(0.205, 1.05, r"$1\sigma$", color="gray", fontsize=9, va="bottom", ha="right")
ax.text(0.205, 2.05, r"$2\sigma$", color="gray", fontsize=9, va="bottom", ha="right")
ax.text(0.205, 3.05, r"$3\sigma$", color="gray", fontsize=9, va="bottom", ha="right")
ax.set_xlabel(r"Redshift cut $z_{\rm cut}$", fontsize=13)
ax.set_ylabel(r"Deviation from $\Lambda$CDM [$\sigma$]", fontsize=12)
ax.set_title(r"(c)  2D $(w_0, w_a)$ significance against $\Lambda$CDM", fontsize=14, fontweight="bold", loc="left")
ax.legend(fontsize=11, loc="upper right", framealpha=0.9)
ax.set_xlim(-0.005, 0.205)
ax.set_ylim(bottom=-0.1)
ax.grid(True, alpha=0.2)

# ── Panel (d): low-redshift distributions of the three samples ──
ax = axes[1, 1]
ax.set_title(r"(d)  Low-redshift distributions of the three samples",
             fontsize=14, fontweight="bold", loc="left")

def load_z(label):
    f = SN_DATA[label]
    if label == "Pantheon+":
        return np.loadtxt(f, skiprows=1, usecols=(4,))          # zCMB
    return np.loadtxt(f, delimiter=",", skiprows=1, usecols=(2,))  # zCMB (DY5) / zHD (Dovekie)

zbin = np.arange(0.0, 0.31, 0.01)
for label in SN_DATA:
    z = load_z(label)
    ax.hist(z, bins=zbin, histtype="step", linewidth=2.2, color=colors[label],
            label=f"{label}  (N={len(z)})")
for zc in [0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.20]:
    ax.axvline(zc, color="gray", linestyle=":", linewidth=1.0, alpha=0.45)
ax.axvspan(0.0, 0.12, color="gray", alpha=0.06)
ax.text(0.012, 0.94, "z ≤ 0.12 excised", transform=ax.transAxes,
        fontsize=10, color="0.35", va="top")
ax.set_xlabel("Redshift $z$", fontsize=13)
ax.set_ylabel("Number of SNe", fontsize=13)
ax.set_xlim(0.0, 0.30)
ax.set_yscale("log")
ax.legend(fontsize=11, loc="upper right", framealpha=0.9)
ax.grid(True, alpha=0.2)

fig.suptitle("CPL Dark Energy Constraints Under Redshift Truncation", fontsize=17, y=1.01, fontweight="bold")
fig.tight_layout()

import os
os.makedirs("plots", exist_ok=True)
fig.savefig("plots/paper_final_fig1.pdf", dpi=200, bbox_inches="tight")
fig.savefig("plots/paper_final_fig1.png", dpi=200, bbox_inches="tight")
print("Saved: plots/paper_final_fig1.pdf")
print("Saved: plots/paper_final_fig1.png")

# ── LaTeX table ──
print()
print(r"\begin{table*}[t]")
print(r"\caption{CPL truncation scan results. All use SN+DESI DR2 BAO+full Planck 2018 TTTEEE\_lite. $\sigma$ is the 2D $(w_0,w_a)$ deviation from $\Lambda$CDM.}")
print(r"\label{tab:cpl}")
print(r"\begin{tabular}{c|ccc|ccc|ccc}")
print(r"\hline")
print(r"$z_{\rm cut}$ & PP $w_0$ & PP $w_a$ & PP $\sigma$ & DY5 $w_0$ & DY5 $w_a$ & DY5 $\sigma$ & Dov $w_0$ & Dov $w_a$ & Dov $\sigma$ \\")
print(r"\hline")
for i, zc in enumerate(zcuts):
    def fmt(v, e, prec=3):
        ep = int(round(e * 10**prec))
        return f"${v:+.{prec}f}({ep})$" if ep < 100 else f"${v:+.{prec}f}$"
    print(f"{zc:.2f} & {pp_w0[i]:+.3f}({int(round(pp_w0e[i]*1000))}) & {pp_wa[i]:+.3f}({int(round(pp_wae[i]*1000))}) & ${pp_sig[i]:.1f}\\sigma$ & {dy5_w0[i]:+.3f}({int(round(dy5_w0e[i]*1000))}) & {dy5_wa[i]:+.3f}({int(round(dy5_wae[i]*1000))}) & ${dy5_sig[i]:.1f}\\sigma$ & {dov_w0[i]:+.3f}({int(round(dov_w0e[i]*1000))}) & {dov_wa[i]:+.3f}({int(round(dov_wae[i]*1000))}) & ${dov_sig[i]:.1f}\\sigma$ \\\\")
print(r"\hline")
print(r"\end{tabular}")
print(r"\end{table*}")

# Key numbers
print()
print("=== KEY ABSTRACT NUMBERS ===")
print(f"PP w0 range (z=0 to z=0.15): {pp_w0[:12].min():+.3f} to {pp_w0[:12].max():+.3f}")
print(f"DY5 sig drop: {dy5_sig[0]:.1f}σ (z=0) → {dy5_sig[5]:.1f}σ (z=0.08)")
print(f"Dovekie sig drop: {dov_sig[0]:.1f}σ (z=0) → {dov_sig[9]:.1f}σ (z=0.10)")
print(f"PP stable sig at z>=0.05: {np.mean(pp_sig[2:12]):.1f}+-{np.std(pp_sig[2:12]):.1f}σ")
print(f"Decomposition: {1.1}+{1.8}+{1.9} = {1.1+1.8+1.9:.1f}σ")
