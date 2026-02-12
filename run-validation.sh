#!/bin/bash
#
# Wrapper script for partis-zenodo validation pipeline
# Usage: ./run-validation.sh --actions <actions> [options]
#

set -e

# Defaults
BASEDIR="${BASEDIR:-/home/drich/matsen-lab/bcr-larch/partis-zenodo-1}"
ACTIONS=""
N_SUB_PROCS=15
N_MAX_PROCS=5
N_REPLICATES=3
OBS_TIMES="15:20:30:40:50"
N_SIM_EVENTS=70

usage() {
    cat <<EOF
Usage: $(basename "$0") --actions <actions> [options]

Required:
  --actions ACTIONS       Colon-separated actions (e.g., gctree:igphyml or tree-perf)

Options:
  --base-outdir DIR       Output directory (default: \$BASEDIR or $BASEDIR)
  --n-sub-procs N         Number of sub-processes (default: $N_SUB_PROCS)
  --n-max-procs N         Max parallel processes (default: $N_MAX_PROCS)
  --n-replicates N        Number of replicates/seeds (default: $N_REPLICATES)
  --obs-times LIST        Colon-separated obs-times (default: $OBS_TIMES)
  --n-sim-events N        Number of simulated events (default: $N_SIM_EVENTS)
  --dry-run               Print command without executing
  -h, --help              Show this help

Actions:
  Phase 1: simu
  Phase 2: cache-parameters, partition, write-fake-paired-annotations,
           iqtree, raxml, gctree, gctree-mut-mult, gctree-no-dag, igphyml, bcrlarch-pars
  Phase 3: tree-perf (runs all *-tree-perf actions)
           Or use individual: iqtree-tree-perf, raxml-tree-perf, etc.
  Phase 4: plot, combine-plots

Action Shortcuts:
  --tree-perf-basic       Run tree-perf for: iqtree, raxml, gctree, gctree-mut-mult, bcrlarch-pars
                          (skips gctree-no-dag and igphyml)

Examples:
  $(basename "$0") --actions gctree:gctree-mut-mult:igphyml
  $(basename "$0") --actions tree-perf
  $(basename "$0") --tree-perf-basic                         # Skip gctree-no-dag and igphyml
  $(basename "$0") --actions plot:combine-plots
  $(basename "$0") --actions gctree --base-outdir /path/to/output --dry-run
EOF
    exit 0
}

# Parse arguments
DRY_RUN=false
while [[ $# -gt 0 ]]; do
    case $1 in
        --actions)
            ACTIONS="$2"
            shift 2
            ;;
        --tree-perf-basic)
            ACTIONS="iqtree-tree-perf:raxml-tree-perf:gctree-tree-perf:gctree-mut-mult-tree-perf:bcrlarch-pars-tree-perf"
            shift
            ;;
        --base-outdir)
            BASEDIR="$2"
            shift 2
            ;;
        --n-sub-procs)
            N_SUB_PROCS="$2"
            shift 2
            ;;
        --n-max-procs)
            N_MAX_PROCS="$2"
            shift 2
            ;;
        --n-replicates)
            N_REPLICATES="$2"
            shift 2
            ;;
        --obs-times)
            OBS_TIMES="$2"
            shift 2
            ;;
        --n-sim-events)
            N_SIM_EVENTS="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

if [[ -z "$ACTIONS" ]]; then
    echo "Error: --actions is required"
    usage
fi

# Build the command
CMD="./test/cf-paired-loci.py \\
    --label gct-valid \\
    --version v6 \\
    --n-replicates $N_REPLICATES \\
    --obs-times-list $OBS_TIMES \\
    --n-sim-events-list $N_SIM_EVENTS \\
    --carry-cap-list 1000 \\
    --simu-type bcr-phylo \\
    --perf-metrics coar:rf:mrca \\
    --calc-antns \\
    --inference-extra-args=\"--no-indels --simultaneous-true-clonal-seqs\" \\
    --plot-metrics tree-perf \\
    --final-plot-xvar obs-times \\
    --simu-extra-args=\"--target-distance 10 --context-depend 1 --tdist-weights random-uniform --min-target-distance 2 --n-sim-seqs-per-generation 89 --parameter-variances n-sim-seqs-per-generation,23 --aa-paratope-positions N=60 --aa-struct-positions N=100 --leaf-sampling-scheme high-affinity --n-naive-seq-copies 100\" \\
    --n-sub-procs $N_SUB_PROCS \\
    --n-max-procs $N_MAX_PROCS \\
    --single-light-locus igk \\
    --base-outdir $BASEDIR \\
    --actions $ACTIONS"

echo "Running: $CMD"
echo ""

if [[ "$DRY_RUN" == "true" ]]; then
    echo "(dry-run mode - command not executed)"
    exit 0
fi

# Execute
eval "$CMD"
