"""Profile-likelihood cross-check of the 2D (w0, wa) significances.

For a given chain: fix (w0, wa) on a grid, optimize the remaining 7
cosmological params, build the profile Delta-chi2 surface, and compare
sigma_profile = sqrt(min Delta-chi2 at LambdaCDM - global min) with the
MCMC-covariance sigma used in the paper.

Usage:
  python profile_scan.py test <prefix>            # 3 points, timing check
  python profile_scan.py run  <prefix> [--grid N] # full grid + LambdaCDM
  python profile_scan.py all                      # the 3 representative chains
"""
import json
import glob
import re
import sys
import time
import os

import numpy as np
import scipy.optimize as opt

FREE = ["ombh2", "omch2", "H0", "tau", "As", "ns", "A_planck"]
FREE_BOUNDS = {
    "ombh2": (0.005, 0.04), "omch2": (0.05, 0.4), "H0": (50.0, 90.0),
    "tau": (0.001, 0.5), "As": (5e-10, 5e-9), "ns": (0.8, 1.2),
    "A_planck": (0.8, 1.3),
}


def load_chain(prefix):
    files = sorted(f for f in glob.glob(prefix + ".*.txt")
                   if re.search(r"\.\d+\.txt$", f))
    blocks = []
    for f in files:
        with open(f) as fh:
            cols = fh.readline().split()[1:]
        blocks.append(np.loadtxt(f, comments="#"))
    d = np.concatenate(blocks)
    return d, cols


def chain_stats(prefix):
    d, cols = load_chain(prefix)
    wi, wai = cols.index("w"), cols.index("wa")
    w, wa = d[:, wi], d[:, wai]
    raw = {c: d[:, i].mean() for i, c in enumerate(cols)}
    raw_std = {c: d[:, i].std() for i, c in enumerate(cols)}
    min_mlp = d[:, cols.index("minuslogpost")].min()
    map_row = d[np.argmin(d[:, cols.index("minuslogpost")])]
    map_vals = {c: map_row[i] for i, c in enumerate(cols)}
    means = {k: raw[k] for k in FREE}
    sigmas = {k: raw_std[k] for k in FREE}
    cov_ww = np.cov(np.vstack([w, wa]))
    return (means, sigmas, w.mean(), w.std(), wa.mean(), wa.std(), cov_ww,
            min_mlp, map_vals)


def build_model(updated_yaml, w0_val, wa_val):
    import yaml
    from cobaya.model import get_model
    cfg = yaml.safe_load(open(updated_yaml, encoding="utf-8"))
    for k in ["sampler", "output", "resume", "force", "version",
              "checkpoint", "covmat", "post", "debug", "stop_at_error"]:
        cfg.pop(k, None)
    params = cfg.setdefault("params", {})
    # cobaya 3.6.2 fixed-param ("value") transfer to CAMB is buggy; emulate
    # fixing with an ultra-narrow box prior (constant logprior offset across
    # all grid points, so relative Delta-chi2 is unaffected)
    if w0_val is not None:
        params["w"] = {"prior": {"min": w0_val - 5e-5, "max": w0_val + 5e-5},
                       "latex": "w_0"}
    if wa_val is not None:
        params["wa"] = {"prior": {"min": wa_val - 5e-5, "max": wa_val + 5e-5},
                        "latex": "w_a"}
    return get_model(cfg)


def logpost(model, pdict):
    """Cobaya 3.6 model.logpost returns the log-posterior (larger = better)."""
    res = model.logpost(pdict)
    return res[0] if isinstance(res, (tuple, list)) else res


LOG_PARAMS = {"As", "A_planck"}


def log_bounds(k):
    lo, hi = FREE_BOUNDS[k]
    if k in LOG_PARAMS:
        return (np.log10(lo), np.log10(hi))
    return (lo, hi)


def unscale(s, lo, hi):
    return lo + s * (hi - lo)


def _norm_bounds(x0, means, sigmas, n_sig=3):
    """Narrow bounds: chain mean +/- n_sig sigma (clipped to physical bounds)."""
    lo = np.array([max(FREE_BOUNDS[k][0], means[k] - n_sig * sigmas[k])
                   for k in FREE])
    hi = np.array([min(FREE_BOUNDS[k][1], means[k] + n_sig * sigmas[k])
                   for k in FREE])
    xl = np.array([np.log10(v) if FREE[i] in LOG_PARAMS else v
                   for i, v in enumerate(x0)])
    for i, k in enumerate(FREE):
        if k in LOG_PARAMS:
            lo[i], hi[i] = np.log10(lo[i]), np.log10(hi[i])
            xl[i] = np.log10(x0[i])
    return lo, hi, xl


def profile_one(model, x0, means, sigmas, min_mlp, w0v=None, wav=None,
                maxiter=400):
    """Nelder-Mead in normalized [0,1]^7 space maximizing logpost;
    penalty outside bounds and on non-physical solutions (CAMB failures
    return absurdly negative logpost)."""
    lo, hi, xl = _norm_bounds(x0, means, sigmas)
    s0 = (xl - lo) / (hi - lo)
    # chain minimum minuslogpost -> maximum logpost floor
    lo_threshold = -(min_mlp + 100.0)

    def f(s):
        # project into the box: the 3-sigma bounds are the feasible region
        sc = np.clip(s, 0.0, 1.0)
        x = lo + sc * (hi - lo)
        vals = {FREE[i]: (10 ** x[i] if FREE[i] in LOG_PARAMS else x[i])
                for i in range(len(FREE))}
        if w0v is not None:
            vals["w"] = w0v
        if wav is not None:
            vals["wa"] = wav
        r = logpost(model, vals)
        if not np.isfinite(r) or r < lo_threshold:
            return 1e6
        return -r

    t0 = time.time()
    res = opt.minimize(f, s0, method="Nelder-Mead",
                       options={"maxiter": maxiter, "xatol": 1e-3, "fatol": 1e-2})
    return res.fun, res.success, time.time() - t0, res.nit


REPRESENTATIVE = {
    "dy5_z0p00":   r"D:\research\cosmology\truncation_sn\chains\fp_dy5_z0p00",
    "pp_z0p00":    r"D:\research\cosmology\truncation_sn\chains\fp_pp_z0p00",
    "dovekie_z0p08": r"D:\research\cosmology\truncation_sn\dovekie_full_planck\fp_dovekie_z0p08",
}


def yaml_for(prefix):
    return prefix + ".updated.yaml"


def _point_worker(args):
    prefix, w0v, wav, means, sigmas, min_mlp = args
    # start from the chain sample nearest to the grid point: chain samples are
    # guaranteed feasible (CAMB evaluated fine), avoiding MAP-params traps
    d, cols = load_chain(prefix)
    wi, wai = cols.index("w"), cols.index("wa")
    dist = (d[:, wi] - w0v) ** 2 + (d[:, wai] - wav) ** 2
    row = d[np.argmin(dist)]
    x0 = np.array([row[cols.index(k)] for k in FREE])
    model = build_model(yaml_for(prefix), w0v, wav)
    mlp, ok, dt, nit = profile_one(model, x0, means, sigmas, min_mlp, w0v, wav)
    return w0v, wav, mlp, ok, dt, nit


def run_chain(prefix, grid_n=5, test=False, nproc=6):
    name = os.path.basename(prefix)
    means, sigmas, w0m, w0s, wam, was, cov_ww, min_mlp, map_vals = chain_stats(prefix)
    print(f"[{name}] N={sum(1 for _ in glob.glob(prefix + '.*.txt'))} files "
          f"w0={w0m:.4f}+-{w0s:.4f} wa={wam:.4f}+-{was:.4f}")

    if test:
        # center, center+1sigma(w0), LambdaCDM
        pts = [(w0m, wam), (w0m + w0s, wam), (-1.0, 0.0)]
    else:
        w0g = np.linspace(w0m - 1.75 * w0s, w0m + 1.75 * w0s, grid_n)
        wag = np.linspace(wam - 1.75 * was, wam + 1.75 * was, grid_n)
        pts = [(w0, wa) for w0 in w0g for wa in wag] + [(-1.0, 0.0)]

    jobs = [(prefix, w0v, wav, means, sigmas, min_mlp) for w0v, wav in pts]
    results = {}
    if test:
        for j in jobs:
            w0v, wav, mlp, ok, dt, nit = _point_worker(j)
            key = f"{w0v:.4f}_{wav:.4f}"
            results[key] = {"w0": w0v, "wa": wav, "minlogpost": mlp,
                            "success": ok, "sec": dt, "nit": nit}
            print(f"  ({w0v:+.3f},{wav:+.3f}) mlp={mlp:.3f} ok={ok} "
                  f"{dt:.0f}s/{nit}it", flush=True)
    else:
        import multiprocessing as mp
        with mp.Pool(nproc) as pool:
            for w0v, wav, mlp, ok, dt, nit in pool.imap_unordered(_point_worker, jobs):
                key = f"{w0v:.4f}_{wav:.4f}"
                results[key] = {"w0": w0v, "wa": wav, "minlogpost": mlp,
                                "success": ok, "sec": dt, "nit": nit}
                print(f"  ({w0v:+.3f},{wav:+.3f}) mlp={mlp:.3f} ok={ok} "
                      f"{dt:.0f}s/{nit}it", flush=True)
        pool.close()

    return results, name


def analyze(name, prefix, results, w0m, wam):
    is_lcdm = lambda k: abs(results[k]["w0"] + 1) < 1e-6 and abs(results[k]["wa"]) < 1e-6
    best = min((k for k in results if not is_lcdm(k)),
               key=lambda k: results[k]["minlogpost"])
    chi2 = {k: 2 * (results[k]["minlogpost"] - results[best]["minlogpost"])
            for k in results}
    lcdm_key = None
    for k in results:
        if abs(results[k]["w0"] + 1) < 1e-6 and abs(results[k]["wa"]) < 1e-6:
            lcdm_key = k
    if lcdm_key is None:
        return None
    sig_profile = np.sqrt(max(chi2[lcdm_key], 0.0))
    d, cols = load_chain(prefix)
    wi, wai = cols.index("w"), cols.index("wa")
    w, wa = d[:, wi], d[:, wai]
    cov = np.cov(np.vstack([w, wa]))
    delta = np.array([w.mean() + 1, wa.mean()])
    sig_mcmc = np.sqrt(delta @ np.linalg.inv(cov) @ delta)
    print(f"\n[{name}] MCMC sigma = {sig_mcmc:.3f} | "
          f"PROFILE sigma = {sig_profile:.3f} | delta = {sig_mcmc - sig_profile:+.3f}")
    return {"name": name, "w0_mean": w0m, "wa_mean": wam,
            "sig_mcmc": sig_mcmc, "sig_profile": sig_profile,
            "delta_chi2_lcdm": chi2[lcdm_key],
            "grid": {k: {"w0": results[k]["w0"], "wa": results[k]["wa"],
                         "delta_chi2": chi2[k]} for k in results}}


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    out = []
    if mode == "all":
        jobs = list(REPRESENTATIVE.items())
    elif mode == "test":
        jobs = [(sys.argv[2], sys.argv[2])] if len(sys.argv) > 2 else [("dy5_z0p00", REPRESENTATIVE["dy5_z0p00"])]
    else:
        prefix = sys.argv[1]
        jobs = [(os.path.basename(prefix), prefix)]

    for short, prefix in jobs:
        grid_n = 5
        results, name = run_chain(prefix, grid_n=grid_n, test=(mode == "test"))
        d, cols = load_chain(prefix)
        wi, wai = cols.index("w"), cols.index("wa")
        entry = analyze(name, prefix, results, d[:, wi].mean(), d[:, wai].mean())
        if entry:
            out.append(entry)

    with open(r"D:\research\cosmology\truncation_sn\profile_results.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\nSaved: profile_results.json")
