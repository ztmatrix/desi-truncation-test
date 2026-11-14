"""
Download SN data files required by Cobaya's built-in likelihoods.

The data is downloaded from the CobayaSampler/sn_data GitHub repository.
If you have a VPN or proxy, set the HTTP_PROXY / HTTPS_PROXY env vars
before running this script.

After downloading, the data will be at:
    <repository>/packages/data/sn_data/

Set COBAYA_PACKAGES_PATH to <repository>/packages before running Cobaya,
so that its built-in SN likelihoods find the data.

Then you can use Cobaya's built-in sn.pantheonplus and sn.desy5
likelihoods (or the Truncated* wrappers in truncation_sn.py).
"""

import os
import sys
import urllib.request
import tarfile
import shutil
import socket

# ── Config ──────────────────────────────────────────────────────────
DATA_URL = (
    "https://github.com/CobayaSampler/sn_data/archive/refs/tags/v1.8.tar.gz"
)
PACKAGES_DIR = os.path.join(os.path.dirname(__file__), "..", "packages")
DATA_DIR = os.path.join(PACKAGES_DIR, "data", "sn_data")
# ─────────────────────────────────────────────────────────────────────

os.makedirs(DATA_DIR, exist_ok=True)

# Check if already installed
required_files = ["pantheonplus.dat", "pantheonplus.mag_covmat",
                  "desy5.dat", "desy5.mag_covmat"]
existing = [f for f in required_files
            if os.path.exists(os.path.join(DATA_DIR, f))]
if len(existing) == len(required_files):
    print("All data files already present. Nothing to download.")
    print(f"Data directory: {DATA_DIR}")
    for f in required_files:
        full = os.path.join(DATA_DIR, f)
        size_mb = os.path.getsize(full) / (1024*1024)
        print(f"  {f}  ({size_mb:.1f} MB)")
    sys.exit(0)
elif existing:
    print(f"Partial data found ({len(existing)}/{len(required_files)} files).")
    print("Re-downloading all data...")
    shutil.rmtree(DATA_DIR)
    os.makedirs(DATA_DIR, exist_ok=True)

print(f"Downloading SN data from GitHub...")
print(f"URL: {DATA_URL}")
print(f"Target: {DATA_DIR}")
print()

socket.setdefaulttimeout(120)

tarball_path = os.path.join(PACKAGES_DIR, "data", "sn_data.tar.gz")

try:
    print("Downloading (this may take several minutes in China)...")
    urllib.request.urlretrieve(DATA_URL, tarball_path)
    file_size_mb = os.path.getsize(tarball_path) / (1024*1024)
    print(f"Downloaded {file_size_mb:.1f} MB")
except Exception as e:
    print(f"\nAutomatic download failed: {e}")
    print("\n── Manual Download Instructions ──")
    print(f"1. Open this URL in your browser (with VPN if needed):")
    print(f"   {DATA_URL}")
    print(f"2. Save the file to:")
    print(f"   {tarball_path}")
    print(f"3. Re-run this script to extract the data.")
    sys.exit(1)

# Extract
print("Extracting...")
with tarfile.open(tarball_path, "r:gz") as tar:
    # The tarball contains a top-level directory like "sn_data-1.8/"
    members = tar.getmembers()
    top_dir = members[0].name.split("/")[0]

    # Extract to temp location
    temp_dir = os.path.join(PACKAGES_DIR, "data", "_temp_extract")
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    tar.extractall(path=temp_dir)

    # Move files from the extracted top-level dir to DATA_DIR
    extracted_content = os.path.join(temp_dir, top_dir)
    for item in os.listdir(extracted_content):
        src = os.path.join(extracted_content, item)
        dst = os.path.join(DATA_DIR, item)
        if os.path.isfile(src):
            shutil.move(src, dst)

    shutil.rmtree(temp_dir)

# Clean up tarball
os.remove(tarball_path)

# Verify
print("\nInstalled data files:")
for f in required_files:
    full = os.path.join(DATA_DIR, f)
    if os.path.exists(full):
        size_mb = os.path.getsize(full) / (1024*1024)
        print(f"  OK  {f}  ({size_mb:.1f} MB)")
    else:
        print(f"  MISSING  {f}")

print("\nSetup complete. You can now run Cobaya with the SN likelihoods.")
