"""
Injection-recovery simulation for the redshift truncation method.

Injects a known magnitude offset into low-z Pantheon+ SNe and verifies
that the truncation scan can recover the true cosmology.

Usage:
    python sim_inject.py --setup     # Generate injected data files and configs
    python sim_inject.py --run       # Run the MCMC chains
    python sim_inject.py --plot      # Plot the recovery results
"""

import os
import sys
import shutil
import subprocess
import numpy as np
from pathlib import Path

BASE_DIR = Path(__file__).parent
PYTHON = r"C:\Users\ztmat\AppData\Local\Python\bin\python3.exe"
Z_CUTS = [0.00, 0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.20]

# Injection parameters
INJECT_Z_MAX = 0.10       # inject offset in SNe with z <= this
INJECT_DMAG = 0.04        # magnitude offset to inject (matching DESY5 ~0.043 mag)
INJECT_SEED = 42

# Data paths
SN_DATA_DIR = BASE_DIR.parent / "packages" / "data" / "sn_data"
PP_DATA_FILE = SN_DATA_DIR / "PantheonPlus" / "Pantheon+SH0ES.dat"
PP_COV_FILE = SN_DATA_DIR / "PantheonPlus" / "Pantheon+SH0ES_STAT+SYS.cov"
INJECTED_DIR = BASE_DIR / "injected_data"


def setup_injected_data():
    """Create injected copies of the Pantheon+ data."""
    os.makedirs(INJECTED_DIR, exist_ok=True)

    # Read original data line by line (handles mixed string/numeric columns)
    print(f"Reading original Pantheon+ data...")
    with open(PP_DATA_FILE) as f:
        lines = f.readlines()
    header = lines[0].split()
    col_idx = {name.lower(): i for i, name in enumerate(header)}
    mag_col = col_idx["m_b_corr"]
    zhd_col = col_idx["zhd"]

    # Parse data rows, tracking z and mag
    z_cmb_list = []
    data_rows = []
    for line in lines[1:]:
        parts = line.split()
        z_cmb_list.append(float(parts[zhd_col]))
        data_rows.append(parts)

    z_cmb = np.array(z_cmb_list)
    n_low = np.sum(z_cmb <= INJECT_Z_MAX)
    n_total = len(z_cmb)
    print(f"Total SNe: {n_total}")
    print(f"Low-z (z <= {INJECT_Z_MAX}): {n_low}")

    # Inject offset into m_b_corr for low-z SNe
    rng = np.random.RandomState(INJECT_SEED)
    for i in range(n_total):
        if z_cmb[i] <= INJECT_Z_MAX:
            orig_mag = float(data_rows[i][mag_col])
            noise = rng.normal(0, 0.002)
            new_mag = orig_mag + INJECT_DMAG + noise
            data_rows[i][mag_col] = f"{new_mag:.6f}"

    # Write injected data file
    injected_file = INJECTED_DIR / "Pantheon+SH0ES_injected.dat"
    with open(injected_file, "w") as f:
        f.write(" ".join(header) + "\n")
        for row in data_rows:
            f.write(" ".join(row) + "\n")
    print(f"Written: {injected_file}")

    # Copy covariance matrix (unchanged)
    cov_dest = INJECTED_DIR / "Pantheon+SH0ES_STAT+SYS.cov"
    shutil.copy(PP_COV_FILE, cov_dest)
    print(f"Copied covariance: {cov_dest}")

    # Create config.dataset for Cobaya
    config_file = INJECTED_DIR / "config.dataset"
    with open(config_file, "w") as f:
        f.write("name = PANTHEONPLUS_INJECTED\n")
        f.write("data_file = injected_data/Pantheon+SH0ES_injected.dat\n")
        f.write("mag_covmat_file = injected_data/Pantheon+SH0ES_STAT+SYS.cov\n")
    print(f"Written: {config_file}")

    return injected_file, cov_dest


def generate_injected_configs():
    """Generate Cobaya YAML configs using the injected data."""
    config_dir = BASE_DIR / "configs_inject"
    os.makedirs(config_dir, exist_ok=True)

    configs = []
    for zc in Z_CUTS:
        zs = f"{zc:.2f}".replace(".", "p")
        config_path = config_dir / f"inject_z{zs}_bao.yaml"
        output_dir = f"chains/inject_z{zs}_bao"

        content = f"""# Injection-recovery: offset={INJECT_DMAG} mag, z<={INJECT_Z_MAX}
theory:
  camb:
    extra_args:
      lens_potential_accuracy: 1
      dark_energy_model: ppf

likelihood:
  truncation_sn.TruncatedInjectedPP:
    zcut: {zc}
    data_path: {INJECTED_DIR.as_posix()}
  bao.desi_dr2:

params:
  ombh2:
    prior: {{min: 0.005, max: 0.04}}
    ref: 0.0224
    proposal: 0.0001
  omch2:
    prior: {{min: 0.01, max: 0.99}}
    ref: 0.12
    proposal: 0.001
  H0:
    prior: {{min: 20, max: 100}}
    ref: 67.4
    proposal: 0.5
  tau:
    prior:
      dist: norm
      loc: 0.054
      scale: 0.007
    ref: 0.054
    proposal: 0.003
  As:
    prior: {{min: 2, max: 4}}
    ref: 3.05
    proposal: 0.01
  ns:
    prior: {{min: 0.8, max: 1.2}}
    ref: 0.965
    proposal: 0.003
  w:
    prior: {{min: -3, max: 1}}
    ref: -1.0
    proposal: 0.05
    latex: w_0

sampler:
  mcmc:
    drag: true
    max_tries: 1000
    Rminus1_stop: 0.01
    Rminus1_cl_stop: 0.03
    max_samples: 10000

output: {output_dir}
"""
        with open(config_path, "w") as f:
            f.write(content)
        configs.append((config_path, output_dir))

    print(f"Generated {len(configs)} injected configs in {config_dir}/")
    return configs


def run_chains(configs):
    """Run all injected MCMC chains."""
    env = os.environ.copy()
    env["PYTHONPATH"] = str(BASE_DIR) + (
        os.pathsep + env.get("PYTHONPATH", "")
    )

    for config_path, output_dir in configs:
        print(f"\nRunning: {config_path.name}")
        subprocess.run(
            [PYTHON, "-m", "cobaya.run", str(config_path)],
            cwd=str(BASE_DIR), env=env,
        )


def plot_results():
    """Plot injection-recovery results."""
    from getdist import loadMCSamples
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    zcuts = Z_CUTS
    inject_w = []
    inject_err = []

    for zc in zcuts:
        zs = f"{zc:.2f}".replace(".", "p")
        prefix = f"chains/inject_z{zs}_bao"
        try:
            samples = loadMCSamples(prefix, settings={"ignore_rows": 0.3})
            w = samples.getParams().w
            inject_w.append(w.mean())
            inject_err.append(w.std())
        except Exception as e:
            print(f"  Failed: {prefix}: {e}")
            inject_w.append(np.nan)
            inject_err.append(np.nan)

    inject_w = np.array(inject_w)
    inject_err = np.array(inject_err)

    # Compare with clean PP results
    pp_bao_w = np.array([-0.9146, -0.9182, -0.9191, -0.9092, -0.9135,
                         -0.9139, -0.9161, -0.8961])
    pp_bao_err = np.array([0.0393, 0.0409, 0.0488, 0.0555, 0.0540,
                          0.0550, 0.0563, 0.0603])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    # Left: w vs z_cut
    ax1.errorbar(zcuts, pp_bao_w, yerr=pp_bao_err, color="#2166ac",
                 marker="s", linewidth=2, markersize=8, capsize=4,
                 label="Pantheon+ (clean)")
    ax1.errorbar(zcuts, inject_w, yerr=inject_err, color="#d6604d",
                 marker="^", linewidth=2, markersize=8, capsize=4,
                 label=f"Pantheon+ injected ({INJECT_DMAG} mag, z<{INJECT_Z_MAX})")
    ax1.axhline(y=-1.0, color="gray", linestyle="--", linewidth=1.5, alpha=0.7)
    ax1.axvline(x=INJECT_Z_MAX, color="red", linestyle=":", alpha=0.5, linewidth=1)
    ax1.set_xlabel(r"Redshift cut $z_{\rm cut}$", fontsize=12)
    ax1.set_ylabel(r"$w$", fontsize=12)
    ax1.set_title("Injection-Recovery: $w$ vs $z_{\\rm cut}$", fontsize=13)
    ax1.legend(fontsize=10, loc="lower left")
    ax1.set_xlim(-0.01, 0.21)
    ax1.grid(True, alpha=0.3)

    # Right: significance
    inject_sig = np.abs(inject_w + 1.0) / inject_err
    pp_sig = np.abs(pp_bao_w + 1.0) / pp_bao_err
    ax2.plot(zcuts, pp_sig, color="#2166ac", marker="s", linewidth=2,
             markersize=8, label="Pantheon+ (clean)")
    ax2.plot(zcuts, inject_sig, color="#d6604d", marker="^", linewidth=2,
             markersize=8, label=f"Injected ({INJECT_DMAG} mag)")
    ax2.axvline(x=INJECT_Z_MAX, color="red", linestyle=":", alpha=0.5, linewidth=1)
    ax2.axhline(y=1.0, color="gray", linestyle=":", alpha=0.6)
    ax2.axhline(y=2.0, color="gray", linestyle="--", alpha=0.6)
    ax2.set_xlabel(r"Redshift cut $z_{\rm cut}$", fontsize=12)
    ax2.set_ylabel(r"Deviation from $\Lambda$CDM [$\sigma$]", fontsize=12)
    ax2.set_title("Significance: Injection-Recovery", fontsize=13)
    ax2.legend(fontsize=10, loc="upper right")
    ax2.set_xlim(-0.01, 0.21)
    ax2.set_ylim(bottom=-0.1)
    ax2.grid(True, alpha=0.3)

    fig.suptitle(
        f"Method Validation: Inject {INJECT_DMAG} mag offset, recover via $z_{{\\rm cut}}$",
        fontsize=14, y=1.01)
    fig.tight_layout()
    fig.savefig("plots/sim_inject.pdf", dpi=150, bbox_inches="tight")
    fig.savefig("plots/sim_inject.png", dpi=150, bbox_inches="tight")
    print("Saved: plots/sim_inject.pdf")

    # Print table
    print(f"\n{'z_cut':<8} {'Clean w':<22} {'Injected w':<22} {'Injected sig':<12}")
    print("-" * 65)
    for i in range(len(zcuts)):
        print(f"{zcuts[i]:<8.2f} {pp_bao_w[i]:+.4f} +/- {pp_bao_err[i]:.4f}   "
              f"{inject_w[i]:+.4f} +/- {inject_err[i]:.4f}   {inject_sig[i]:.2f}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--setup", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()

    if args.setup:
        setup_injected_data()
        configs = generate_injected_configs()
        print(f"\nReady. Run with --run to execute {len(configs)} chains.")
    elif args.run:
        configs = generate_injected_configs()
        run_chains(configs)
    elif args.plot:
        plot_results()
    else:
        print("Usage: python sim_inject.py --setup | --run | --plot")
