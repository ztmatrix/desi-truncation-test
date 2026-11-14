"""
Analyze MCMC chains from the redshift truncation scan.

Produces:
  1. Delta_chi2 vs z_cut plot (the key figure)
  2. w0-wa contour comparison across z_cuts
  3. Summary table of significance vs z_cut

Usage:
    python analyze_results.py
"""

import os
import sys
import warnings
from pathlib import Path
import numpy as np
from getdist import MCSamples, loadMCSamples
from getdist import plots as gdplt

BASE_DIR = Path(__file__).parent
CHAIN_DIR = BASE_DIR / "chains"
PLOT_DIR = BASE_DIR / "plots"

# z_cut values used in the scan
Z_CUTS = [0.00, 0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.20]

# Dataset labels
DATASETS = {
    "pp": "Pantheon+",
    "dy5": "DESY5",
}

# LCDM reference: w0=-1, wa=0
LCDM_W0 = -1.0
LCDM_WA = 0.0


def compute_delta_chi2(samples):
    """Compute sqrt(Delta chi2) between best-fit and LCDM.

    For CPL model, this measures how far the best-fit (w0, wa) is
    from (-1, 0) in units of sigma, accounting for correlation.
    """
    w0 = samples.getParams().w
    wa = samples.getParams().wa

    # Get best-fit (maximum likelihood point)
    best_w0 = w0[np.argmax(samples.loglikes)]
    best_wa = wa[np.argmax(samples.loglikes)]

    # Simple Delta chi2 = chi2_LCDM - chi2_best
    # We approximate via the chi2 at the best-fit minus chi2 at LCDM
    # Use the samples' chi2 values
    like_best = np.max(samples.loglikes)

    # Find the sample closest to LCDM
    dist_to_lcdm = (w0 - LCDM_W0)**2 + (wa - LCDM_WA)**2
    idx_lcdm = np.argmin(dist_to_lcdm)
    like_lcdm = samples.loglikes[idx_lcdm]

    delta_chi2 = -2 * (like_lcdm - like_best)
    significance = np.sqrt(max(0, delta_chi2))

    return significance, best_w0, best_wa


def collect_results():
    """Scan chain directories and collect results."""
    results = []

    for chain_dir in sorted(CHAIN_DIR.iterdir()):
        if not chain_dir.is_dir():
            continue

        # Parse directory name: pp_z0p00, dy5_z0p10, etc.
        name = chain_dir.name
        try:
            parts = name.split("_")
            dataset_label = parts[0]
            zcut_str = parts[1]  # e.g., "z0p00"
            zcut = float(zcut_str[1:].replace("p", "."))
        except (IndexError, ValueError):
            print(f"  Skipping unrecognized directory: {name}")
            continue

        # Check for .txt chain files
        chain_files = list(chain_dir.glob("*.txt"))
        if not chain_files:
            print(f"  No chain files found in {name}")
            continue

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                samples = loadMCSamples(
                    str(chain_dir),
                    settings={"ignore_rows": 0.3}
                )
        except Exception as e:
            print(f"  Failed to load {name}: {e}")
            continue

        sig, best_w0, best_wa = compute_delta_chi2(samples)

        results.append({
            "dataset": DATASETS.get(dataset_label, dataset_label),
            "dataset_label": dataset_label,
            "zcut": zcut,
            "significance": sig,
            "best_w0": best_w0,
            "best_wa": best_wa,
            "n_sn": samples.sampler.get("nsn", "?"),
        })

        print(f"  {name}: sig={sig:.2f}σ, w0={best_w0:.3f}, wa={best_wa:.3f}")

    return results


def plot_significance_vs_zcut(results):
    """Create the key figure: significance vs z_cut for each dataset."""
    os.makedirs(PLOT_DIR, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5))

    colors = {"Pantheon+": "#2166ac", "DESY5": "#b2182b"}
    markers = {"Pantheon+": "s", "DESY5": "o"}

    for dataset_name in ["Pantheon+", "DESY5"]:
        pts = [r for r in results if r["dataset"] == dataset_name]
        if not pts:
            continue

        pts.sort(key=lambda r: r["zcut"])
        zcuts = [p["zcut"] for p in pts]
        sigs = [p["significance"] for p in pts]

        ax.plot(
            zcuts, sigs,
            color=colors[dataset_name],
            marker=markers[dataset_name],
            linewidth=2,
            markersize=8,
            label=f"{dataset_name}",
        )

    # Reference lines
    ax.axhline(y=1.0, color="gray", linestyle=":", alpha=0.6)
    ax.axhline(y=2.0, color="gray", linestyle="--", alpha=0.6)
    ax.axhline(y=3.0, color="gray", linestyle="-.", alpha=0.6)

    ax.text(0.21, 1.05, "1σ", color="gray", fontsize=9, va="bottom")
    ax.text(0.21, 2.05, "2σ", color="gray", fontsize=9, va="bottom")
    ax.text(0.21, 3.05, "3σ", color="gray", fontsize=9, va="bottom")

    ax.set_xlabel("Redshift cut $z_{\\rm cut}$", fontsize=12)
    ax.set_ylabel(
        "Deviation from $\\Lambda$CDM  [$\\sigma$]",
        fontsize=12,
    )
    ax.set_title(
        "Stability of $w \\neq -1$ signal under redshift truncation",
        fontsize=13,
    )
    ax.legend(fontsize=11, loc="upper right")
    ax.set_ylim(bottom=-0.2)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(PLOT_DIR / "significance_vs_zcut.pdf", dpi=150)
    fig.savefig(PLOT_DIR / "significance_vs_zcut.png", dpi=150)
    print(f"\nSaved: {PLOT_DIR / 'significance_vs_zcut.pdf'}")
    plt.close(fig)


def plot_w0wa_contours(results):
    """Plot w0-wa contours for selected z_cuts to show the drift."""
    os.makedirs(PLOT_DIR, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    selected_zcuts = [0.00, 0.10, 0.15, 0.20]
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(selected_zcuts)))

    for dataset_label, dataset_name in DATASETS.items():
        fig, ax = plt.subplots(figsize=(7, 6))
        ax.axvline(x=-1.0, color="gray", linestyle="--", alpha=0.5)
        ax.axhline(y=0.0, color="gray", linestyle="--", alpha=0.5)
        ax.plot(-1.0, 0.0, "*", color="black", markersize=12, label="$\\Lambda$CDM")

        for i, zcut in enumerate(selected_zcuts):
            chain_dir = CHAIN_DIR / f"{dataset_label}_z{zcut:.2f}".replace(".", "p")
            if not chain_dir.exists():
                continue

            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    samples = loadMCSamples(
                        str(chain_dir),
                        settings={"ignore_rows": 0.3}
                    )
            except Exception:
                continue

            g = samples.getLikeStats()
            ax.errorbar(
                g.means["w0"], g.means["wa"],
                xerr=g.std_devs["w0"], yerr=g.std_devs["wa"],
                color=colors[i], marker="o", markersize=6,
                label=f"$z_{{\\rm cut}}$ = {zcut:.2f}",
                capsize=3,
            )

        ax.set_xlabel("$w_0$", fontsize=12)
        ax.set_ylabel("$w_a$", fontsize=12)
        ax.set_title(f"{dataset_name}: $w_0$–$w_a$ drift under $z_{{\\rm cut}}$", fontsize=12)
        ax.legend(fontsize=9, loc="upper left")
        ax.grid(True, alpha=0.3)

        fig.tight_layout()
        fig.savefig(PLOT_DIR / f"w0wa_drift_{dataset_label}.pdf", dpi=150)
        fig.savefig(PLOT_DIR / f"w0wa_drift_{dataset_label}.png", dpi=150)
        print(f"Saved: {PLOT_DIR / f'w0wa_drift_{dataset_label}.pdf'}")
        plt.close(fig)


def print_summary_table(results):
    """Print a formatted summary table."""
    print(f"\n{'='*70}")
    print("Summary: Deviation from LCDM [sigma]")
    print(f"{'='*70}")
    print(f"{'z_cut':<10}", end="")
    for ds_name in ["Pantheon+", "DESY5"]:
        print(f"{ds_name:<20}", end="")
    print()
    print("-" * 50)

    for zcut in Z_CUTS:
        print(f"{zcut:<10.2f}", end="")
        for ds_name in ["Pantheon+", "DESY5"]:
            matches = [
                r for r in results
                if r["dataset"] == ds_name and abs(r["zcut"] - zcut) < 0.001
            ]
            if matches:
                r = matches[0]
                bar = "|" * max(1, int(r["significance"] * 5))
                print(f"{r['significance']:<5.2f}σ {bar:<14}", end="")
            else:
                print(f"{'--':<20}", end="")
        print()
    print(f"{'='*70}")


def main():
    print("Analyzing MCMC chains...")
    print(f"Chain directory: {CHAIN_DIR}")
    print()

    results = collect_results()

    if not results:
        print("\nNo results found. Run the MCMC chains first with:")
        print("  python run_scan.py --run")
        return

    print_summary_table(results)
    plot_significance_vs_zcut(results)
    plot_w0wa_contours(results)

    print("\nAnalysis complete.")


if __name__ == "__main__":
    main()
