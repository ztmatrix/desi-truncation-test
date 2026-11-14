"""
Redshift-truncated SN likelihoods for Cobaya.

Usage in YAML config:
    likelihood:
      truncation_sn.TruncatedPantheonPlus:
        zcut: 0.1
      truncation_sn.TruncatedDESy5:
        zcut: 0.15
"""

import os
import numpy as np
from cobaya.likelihood import Likelihood
from cobaya.likelihoods.sn.pantheonplus import PantheonPlus
from cobaya.likelihoods.sn.desy5 import DESy5


class TruncatedPantheonPlus(PantheonPlus):
    """PantheonPlus with configurable redshift truncation.

    SNe with zcmb <= zcut are excluded before computing the likelihood.
    Set zcut=0.0 (default) to include all SNe (same as original).
    """

    zcut: float = 0.0

    def configure(self):
        zmask = self.zcmb > max(self.zcut, 0.01)
        self._apply_mask(zmask=zmask)
        self.pre_vars = 0.0


class TruncatedDESy5(DESy5):
    """DESY5 with configurable redshift truncation.

    SNe with zcmb <= zcut are excluded before computing the likelihood.
    Set zcut=0.0 (default) to include all SNe (same as original DESY5,
    which applies zcmb > 0.01 internally).
    """

    zcut: float = 0.0

    def configure(self):
        zmask = self.zcmb > max(self.zcut, 0.01)
        self._apply_mask(zmask=zmask)
        self.pre_vars = self.mag_err**2


class TruncatedInjectedPP(PantheonPlus):
    """PantheonPlus reading from a custom injected data directory.

    Used for injection-recovery simulations. Set data_path in YAML
    to point to the directory containing the injected data files.

    YAML usage:
        likelihood:
          truncation_sn.TruncatedInjectedPP:
            zcut: 0.10
            data_path: /path/to/injected_data
    """

    zcut: float = 0.0
    data_path: str = ""

    def init_params(self, ini):
        self.twoscriptmfit = False
        # Use custom data path
        custom_path = self.data_path or self.path
        data_file = custom_path + "/Pantheon+SH0ES_injected.dat"
        cov_file = custom_path + "/Pantheon+SH0ES_STAT+SYS.cov"
        self._read_data_file(data_file)
        self.covs = {}
        self.log.debug("Reading covmat for: mag")
        self.covs["mag"] = self._read_covmat(cov_file)
        self.alphabeta_covmat = False
        self.configure()
        self.inverse_covariance_matrix()
        if not self.use_abs_mag:
            self._marginalize_abs_mag()
        self.marginalize = False

    def configure(self):
        zmask = self.zcmb > max(self.zcut, 0.01)
        self._apply_mask(zmask=zmask)
        self.pre_vars = 0.0


class TruncatedDESDovekie(Likelihood):
    """DES-Dovekie Hubble diagram with full inverse covariance matrix.

    Uses the packed upper-triangle inverse covariance (covtot_inv_000.npz)
    with analytic M_B marginalization — same statistical treatment as DESY5.

    YAML usage:
        likelihood:
          truncation_sn.TruncatedDESDovekie:
            zcut: 0.10
            data_file: /path/to/DES-Dovekie_HD.csv
            cov_file: /path/to/covtot_inv_000.npz
    """

    zcut: float = 0.0
    data_file: str = ""
    cov_file: str = ""
    use_abs_mag: bool = False

    def initialize(self):
        import pandas as pd
        # Resolve relative paths against CWD first, then the packages data
        # layout (works both on servers: <cwd>/data/DES-Dovekie, and locally:
        # <module>/../packages/data/sn_data/DES-Dovekie).
        data_dir = os.path.join(
            os.path.dirname(__file__), "..", "packages", "data", "sn_data")

        def _resolve(path, fallback):
            if not path:
                return fallback
            if os.path.isabs(path) or os.path.exists(path):
                return path
            cand = os.path.join(data_dir, *path.split("/")[-2:])
            return cand if os.path.exists(cand) else path

        self.data_file = _resolve(
            self.data_file,
            os.path.join(data_dir, "DES-Dovekie", "DES-Dovekie_HD.csv"))
        self.cov_file = _resolve(
            self.cov_file,
            os.path.join(data_dir, "DES-Dovekie", "covtot_inv_000.npz"))

        df = pd.read_csv(self.data_file)
        self.z_full = df.zHD.values.astype(np.float64)
        self.mu_full = df.MU.values.astype(np.float64)

        # Filter bad measurements (keep if MUERR is reasonable)
        good = df.MUERR.values < 1.0
        self.good_mask = good

        # Build inverse covariance from packed upper triangle
        cov_data = np.load(self.cov_file)
        nsn = int(cov_data["nsn"][0])
        packed = cov_data["cov"]
        full_invcov = np.zeros((nsn, nsn))
        idx = 0
        for i in range(nsn):
            for j in range(i, nsn):
                full_invcov[i, j] = packed[idx]
                full_invcov[j, i] = packed[idx]
                idx += 1

        # Sub-select to good SNe only
        self.invcov_full = full_invcov[np.ix_(good, good)]
        self.nsn = np.sum(good)
        self.log.debug(f"DES-Dovekie: {self.nsn} SNe (from {nsn}) after quality cuts")

    def get_requirements(self):
        reqs = {"angular_diameter_distance": {"z": self.z_full[self.good_mask]}}
        if self.use_abs_mag:
            reqs["Mb"] = None
        return reqs

    def logp(self, **params_values):
        zmask = self.z_full[self.good_mask] > max(self.zcut, 0.01)
        z_use = self.z_full[self.good_mask][zmask]
        mu_obs = self.mu_full[self.good_mask][zmask]

        if np.sum(zmask) == 0:
            return -np.inf

        D_A = self.provider.get_angular_diameter_distance(
            self.z_full[self.good_mask])
        d_L = D_A[zmask] * (1 + z_use)**2
        mu_th = 5.0 * np.log10(d_L) + 25.0

        # Sub-select inverse covariance for truncated SNe
        invcov = self.invcov_full[np.ix_(zmask, zmask)]

        # Analytic M_B marginalization using full inverse covariance
        residuals = mu_th - mu_obs
        A = residuals @ invcov @ residuals
        B = np.sum(invcov @ residuals)
        C = np.sum(invcov)

        chi2 = A - B**2 / C
        return -0.5 * chi2
