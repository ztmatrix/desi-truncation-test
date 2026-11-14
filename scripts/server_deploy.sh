#!/bin/bash
# ============================================================
# Deployment script for truncation_sn on Ubuntu 22.04 server
# Run AFTER extracting the tarball: bash server_deploy.sh
# ============================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== Step 1: System packages ==="
sudo apt-get update -y
sudo apt-get install -y build-essential gfortran libopenmpi-dev openmpi-bin \
    python3 python3-pip python3-venv wget tmux

echo ""
echo "=== Step 2: Python environment ==="
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install numpy scipy matplotlib pyyaml camb cobaya getdist mpi4py pandas

echo ""
echo "=== Step 3: Install clipy ==="
CLIPY_PATH="packages/code/planck_2018/clipy"
if [ -d "$CLIPY_PATH" ] && { [ -f "$CLIPY_PATH/setup.py" ] || [ -f "$CLIPY_PATH/pyproject.toml" ]; }; then
    pip install -e "$CLIPY_PATH"
    python -c "import clipy; print(f'clipy OK')"
else
    echo "ERROR: clipy not found at $CLIPY_PATH"
    exit 1
fi

echo ""
echo "=== Step 4: Check Planck data ==="
PLANCK_CLIK="packages/data/planck_2018/baseline/plc_3.0/hi_l/plik_lite/plik_lite_v22_TTTEEE.clik"
if [ -d "$PLANCK_CLIK" ]; then
    echo "Planck data OK: $PLANCK_CLIK"
else
    echo "ERROR: Planck clik file not found at $PLANCK_CLIK"
    exit 1
fi

echo ""
echo "=== Step 5: Check Dovekie data ==="
if [ -f "data/DES-Dovekie/DES-Dovekie_HD.csv" ] && [ -f "data/DES-Dovekie/covtot_inv_000.npz" ]; then
    echo "Dovekie data OK"
else
    echo "WARNING: Dovekie data missing"
fi

echo ""
echo "=== Step 6: Verify installation ==="
python -c "import camb; print(f'CAMB {camb.__version__} OK')"
python -c "import cobaya; print(f'Cobaya OK')"

echo ""
echo "=== Done! ==="
echo ""
echo "To run all chains:"
echo "  tmux new -s chains"
echo "  source venv/bin/activate"
echo "  bash server_run.sh"
echo ""
echo "To test a single chain first:"
echo "  source venv/bin/activate"
echo "  python -m cobaya.run configs/fp_pp_z0p00.yaml --test"
