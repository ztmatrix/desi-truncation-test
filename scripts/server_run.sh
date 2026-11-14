#!/bin/bash
# ============================================================
# Run all full-Planck truncation chains with auto-shutdown
# Usage: source venv/bin/activate && bash server_run.sh
# ============================================================
set -e  # stop on critical errors (not background jobs)

PROJECT_DIR="$HOME/truncation_sn"
cd "$PROJECT_DIR"
export PYTHONPATH="$PROJECT_DIR"

# ---- Settings (adjust for your instance) ----
MPI_BLOCKS=4          # MPI processes per chain
PARALLEL_CHAINS=4     # chains running simultaneously (16-core: 16/4=4)
SHUTDOWN_WHEN_DONE=true  # auto power-off when all done
MAX_RETRIES=3         # retry failed chains up to N times

export OMP_NUM_THREADS=1
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# ---- Chain list (39 total, 13 zcut × 3 datasets) ----
CHAINS=(
    # Pantheon+
    "fp_pp_z0p00" "fp_pp_z0p02" "fp_pp_z0p05" "fp_pp_z0p06" "fp_pp_z0p07"
    "fp_pp_z0p08" "fp_pp_z0p09" "fp_pp_z0p10" "fp_pp_z0p11" "fp_pp_z0p12"
    "fp_pp_z0p14" "fp_pp_z0p15" "fp_pp_z0p20"
    # DESY5
    "fp_dy5_z0p00" "fp_dy5_z0p02" "fp_dy5_z0p05" "fp_dy5_z0p06" "fp_dy5_z0p07"
    "fp_dy5_z0p08" "fp_dy5_z0p09" "fp_dy5_z0p10" "fp_dy5_z0p11" "fp_dy5_z0p12"
    "fp_dy5_z0p14" "fp_dy5_z0p15" "fp_dy5_z0p20"
    # DES-Dovekie
    "fp_dovekie_z0p00" "fp_dovekie_z0p02" "fp_dovekie_z0p05" "fp_dovekie_z0p06" "fp_dovekie_z0p07"
    "fp_dovekie_z0p08" "fp_dovekie_z0p09" "fp_dovekie_z0p10" "fp_dovekie_z0p11" "fp_dovekie_z0p12"
    "fp_dovekie_z0p14" "fp_dovekie_z0p15" "fp_dovekie_z0p20"
)

# ---- Helper: run a single chain (with resume support) ----
run_chain() {
    local cfg=$1
    echo "[$(date '+%m-%d %H:%M')] Starting $cfg (PID=$$) ..."
    mpirun -np $MPI_BLOCKS --allow-run-as-root --oversubscribe python -m cobaya.run \
        "configs/${cfg}.yaml" > "chains/${cfg}.log" 2>&1
    local rc=$?
    if [ $rc -eq 0 ]; then
        echo "[$(date '+%m-%d %H:%M')] $cfg CONVERGED."
    else
        echo "[$(date '+%m-%d %H:%M')] $cfg stopped (rc=$rc). Checkpoint saved for resume."
        echo "  Failed chain: $cfg (rc=$rc)" >> chains/_failed.txt
    fi
}

# ---- Parallel dispatcher ----
echo "===================================="
echo "Server run started: $(date)"
echo "Cores: $(nproc), MPI blocks/chain: $MPI_BLOCKS, Parallel: $PARALLEL_CHAINS"
echo "Total chains to run: ${#CHAINS[@]}"
echo "===================================="
echo ""

mkdir -p chains

QUEUE=()
RUNNING=0
TOTAL=${#CHAINS[@]}
DONE=0

for cfg in "${CHAINS[@]}"; do
    # Wait if full
    while [ $RUNNING -ge $PARALLEL_CHAINS ]; do
        # Check for finished jobs
        for i in "${!QUEUE[@]}"; do
            pid=${QUEUE[$i]}
            if ! kill -0 $pid 2>/dev/null; then
                wait $pid 2>/dev/null || true
                unset "QUEUE[$i]"
                RUNNING=$((RUNNING - 1))
                DONE=$((DONE + 1))
            fi
        done
        QUEUE=("${QUEUE[@]}")  # compact
        sleep 5
    done

    if [ -f "configs/${cfg}.yaml" ]; then
        run_chain "$cfg" &
        QUEUE+=($!)
        RUNNING=$((RUNNING + 1))
    else
        echo "[$(date '+%H:%M')] SKIP $cfg (no config)"
        DONE=$((DONE + 1))
    fi
done

# Wait for remaining
for pid in "${QUEUE[@]}"; do
    wait $pid 2>/dev/null || true
    DONE=$((DONE + 1))
done

# ---- Retry failed chains ----
RETRY=0
while [ -f chains/_failed.txt ] && [ $RETRY -lt $MAX_RETRIES ]; do
    FAILED=$(cat chains/_failed.txt | sort -u)
    rm -f chains/_failed.txt
    RETRY=$((RETRY + 1))
    FAILED_COUNT=$(echo "$FAILED" | grep -c ':' || true)
    echo ""
    echo "=== Retry batch $RETRY/$MAX_RETRIES ($FAILED_COUNT failed chains) ==="

    QUEUE=(); RUNNING=0
    while IFS= read -r line; do
        [ -z "$line" ] && continue
        cfg=$(echo "$line" | awk -F':' '{print $1}' | xargs)
        [ -z "$cfg" ] && continue

        while [ $RUNNING -ge $PARALLEL_CHAINS ]; do
            for i in "${!QUEUE[@]}"; do
                pid=${QUEUE[$i]}
                if ! kill -0 $pid 2>/dev/null; then
                    wait $pid 2>/dev/null || true
                    unset "QUEUE[$i]"
                    RUNNING=$((RUNNING - 1))
                fi
            done
            QUEUE=("${QUEUE[@]}")
            sleep 5
        done

        if [ -f "configs/${cfg}.yaml" ]; then
            echo "  Retrying $cfg ..."
            run_chain "$cfg" &
            QUEUE+=($!)
            RUNNING=$((RUNNING + 1))
        fi
    done <<< "$FAILED"

    for pid in "${QUEUE[@]}"; do
        wait $pid 2>/dev/null || true
    done
done

if [ -f chains/_failed.txt ]; then
    echo "WARNING: Some chains still failed after $MAX_RETRIES retries:"
    cat chains/_failed.txt
fi

echo ""
echo "===================================="
echo "ALL ${DONE}/${TOTAL} CHAINS COMPLETE"
echo "Finished: $(date)"
echo "Results in: $PROJECT_DIR/chains/"
echo "===================================="

# ---- Auto shutdown ----
if [ "$SHUTDOWN_WHEN_DONE" = true ]; then
    echo ""
    echo "Shutting down server in 60 seconds..."
    echo "Ctrl+C NOW to cancel!"
    sleep 60
    sudo shutdown -h now
fi
