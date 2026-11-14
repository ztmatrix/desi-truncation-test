#!/bin/bash
# ============================================================
# Server 3 (16-core): DY5 tail + DES-Dovekie (16 total)
# Usage: source venv/bin/activate && bash server_run3.sh
# ============================================================
set -e

PROJECT_DIR="$HOME/truncation_sn"
cd "$PROJECT_DIR"
export PYTHONPATH="$PROJECT_DIR"

MPI_BLOCKS=4          # MPI processes per chain
PARALLEL_CHAINS=4     # chains running simultaneously (16-core: 16/4=4)
SHUTDOWN_WHEN_DONE=true
MAX_RETRIES=3

export OMP_NUM_THREADS=1
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

CHAINS=(
    # DES-Dovekie (reverse: high z_cut first)
    "fp_dovekie_z0p20" "fp_dovekie_z0p15" "fp_dovekie_z0p14" "fp_dovekie_z0p12" "fp_dovekie_z0p11"
    "fp_dovekie_z0p10" "fp_dovekie_z0p09" "fp_dovekie_z0p08" "fp_dovekie_z0p07" "fp_dovekie_z0p06"
    "fp_dovekie_z0p05" "fp_dovekie_z0p02" "fp_dovekie_z0p00"
    # DESY5 tail (reverse)
    "fp_dy5_z0p20" "fp_dy5_z0p15" "fp_dy5_z0p14"
)

run_chain() {
    local cfg=$1
    echo "[$(date '+%m-%d %H:%M')] Starting $cfg ..."
    mpirun -np $MPI_BLOCKS --allow-run-as-root --oversubscribe python -m cobaya.run \
        "configs/${cfg}.yaml" --resume > "chains/${cfg}.log" 2>&1
    local rc=$?
    if [ $rc -eq 0 ]; then
        echo "[$(date '+%m-%d %H:%M')] $cfg CONVERGED."
    else
        echo "[$(date '+%m-%d %H:%M')] $cfg stopped (rc=$rc)."
        echo "  Failed chain: $cfg (rc=$rc)" >> chains/_failed.txt
    fi
}

echo "===================================="
echo "Server 3 (16-core): DY5 tail + Dovekie"
echo "Started: $(date)"
echo "Cores: $(nproc), MPI blocks/chain: $MPI_BLOCKS, Parallel: $PARALLEL_CHAINS"
echo "Total chains: ${#CHAINS[@]}"
echo "===================================="
echo ""

mkdir -p chains

QUEUE=(); RUNNING=0; TOTAL=${#CHAINS[@]}; DONE=0

for cfg in "${CHAINS[@]}"; do
    while [ $RUNNING -ge $PARALLEL_CHAINS ]; do
        for i in "${!QUEUE[@]}"; do
            pid=${QUEUE[$i]}
            if ! kill -0 $pid 2>/dev/null; then
                wait $pid 2>/dev/null || true
                unset "QUEUE[$i]"
                RUNNING=$((RUNNING - 1)); DONE=$((DONE + 1))
            fi
        done
        QUEUE=("${QUEUE[@]}")
        sleep 5
    done

    if [ -f "configs/${cfg}.yaml" ]; then
        run_chain "$cfg" &
        QUEUE+=($!); RUNNING=$((RUNNING + 1))
    else
        echo "[$(date '+%H:%M')] SKIP $cfg (no config)"
        DONE=$((DONE + 1))
    fi
done

for pid in "${QUEUE[@]}"; do
    wait $pid 2>/dev/null || true; DONE=$((DONE + 1))
done

# Retry failed chains
RETRY=0
while [ -f chains/_failed.txt ] && [ $RETRY -lt $MAX_RETRIES ]; do
    FAILED=$(sort -u chains/_failed.txt); rm -f chains/_failed.txt
    RETRY=$((RETRY + 1))
    echo ""; echo "=== Retry batch $RETRY/$MAX_RETRIES ==="
    QUEUE=(); RUNNING=0
    while IFS= read -r line; do
        [ -z "$line" ] && continue
        cfg=$(echo "$line" | awk -F':' '{print $1}' | xargs)
        while [ $RUNNING -ge $PARALLEL_CHAINS ]; do
            for i in "${!QUEUE[@]}"; do
                pid=${QUEUE[$i]}
                if ! kill -0 $pid 2>/dev/null; then
                    wait $pid 2>/dev/null || true
                    unset "QUEUE[$i]"; RUNNING=$((RUNNING - 1))
                fi
            done
            QUEUE=("${QUEUE[@]}"); sleep 5
        done
        [ -f "configs/${cfg}.yaml" ] && { run_chain "$cfg" & QUEUE+=($!); RUNNING=$((RUNNING + 1)); }
    done <<< "$FAILED"
    for pid in "${QUEUE[@]}"; do wait $pid 2>/dev/null || true; done
done

echo ""
echo "===================================="
echo "ALL CHAINS COMPLETE: $(date)"
echo "===================================="

if [ "$SHUTDOWN_WHEN_DONE" = true ]; then
    echo "Shutting down in 60s..."; sleep 60; sudo shutdown -h now
fi
