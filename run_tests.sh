#!/usr/bin/env bash
#
# run_tests.sh - build foap4 with chosen options and run its tests.
#
# The script has two phases:
#   1) build   - invoke `make` with the requested compiler/build options
#   2) run     - execute the selected test binaries (via mpirun) and report
#                a pass/fail summary.
#
# Run `./run_tests.sh --help` for the full option list.

set -uo pipefail

# ---------------------------------------------------------------------------
# Locate the repository root (directory holding this script / the Makefile).
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
REPO_ROOT="$SCRIPT_DIR"
cd "$REPO_ROOT" || { echo "error: cannot enter $REPO_ROOT" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
F90C="${F90C:-mpif90}"          # Fortran compiler (Makefile default: mpif90)
OFFLOAD="acc"                   # acc | omp | ompcpu | none
FLOAT_BITS="64"                 # 32 | 64
DEBUG=0                         # DEBUG=1
SAFE=0                          # SAFE=1
BUILDDIR=""                     # empty -> let the Makefile pick build-<tag>
JOBS="$( (nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4) )"

DIMS="both"                     # 2 | 3 | both
NP=1                            # MPI ranks per test
OMP_THREADS=""                  # exported as OMP_NUM_THREADS when set
TIMEOUT=600                     # per-test wall-clock limit in seconds (0 = none)
TEST_OVERRIDE=""                # explicit space separated test list
MPI_RUNNER="mpirun"             # set to "" / "serial" to run without mpirun

NO_BUILD=0
DO_CLEAN=0
FAIL_FAST=0
LIST_ONLY=0
DRY_RUN=0

MAKE_ARGS=()                    # extra NAME=VALUE pairs for make
MPI_ARGS=()                     # extra arguments for mpirun
TEST_ARGS=()                    # extra arguments for every test executable

# Curated test sets.
MAIN_2D="test_refinement_2d test_advection_2d test_xdmf_writer_2d test_euler_2d test_shallow_water_2d"
MAIN_3D="test_refinement_3d test_advection_3d test_xdmf_writer_3d test_euler_3d"

# ---------------------------------------------------------------------------
# Help
# ---------------------------------------------------------------------------
usage() {
cat <<'EOF'
Usage: ./run_tests.sh [OPTIONS] [-- TEST_ARGS...]

Build foap4 with the given options and run the tests.

Build options:
  -c, --compiler CMD      Fortran compiler (F90C)        [default: mpif90]
  -o, --offload MODEL     acc | omp | ompcpu | none      [default: acc]
  -f, --float-bits N      Floating point precision: 32|64 [default: 64]
  -d, --debug             Enable debug flags (DEBUG=1)
  -s, --safe              Use -O2 instead of -Ofast (SAFE=1)
  -B, --builddir DIR      Explicit build directory (BUILDDIR)
  -j, --jobs N            Parallel build jobs            [default: nproc]
      --make-arg ARG      Extra argument passed to make (repeatable),
                          e.g. --make-arg FFLAGS_USER=-mavx2
      --clean             Remove the build directory before building
      --no-build          Skip compilation; run the existing binaries
      --dry-run           Only print the build/run commands

Test selection / execution:
  -D, --dims DIM          2 | 3 | both                   [default: both]
  -n, --np N              MPI ranks per test             [default: 1]
      --threads N         Export OMP_NUM_THREADS=N for the tests
      --tests "T1 T2 ..." Explicit list of test names to run
      --timeout SEC       Per-test timeout, 0 to disable  [default: 600]
      --serial            Run the executables directly (no mpirun)
      --mpirun CMD        MPI launcher command            [default: mpirun]
      --mpi-arg ARG       Extra argument for mpirun (repeatable)
      --fail-fast         Stop at the first failing test
      --list              List the tests that would run and exit
  -h, --help              Show this help

Anything after `--` is forwarded to every test executable, e.g.
    ./run_tests.sh -D 2 -- -viewer=paraview
EOF
}

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -c|--compiler)   F90C="$2"; shift 2 ;;
    -o|--offload)    OFFLOAD="$2"; shift 2 ;;
    -f|--float-bits) FLOAT_BITS="$2"; shift 2 ;;
    -d|--debug)      DEBUG=1; shift ;;
    -s|--safe)       SAFE=1; shift ;;
    -B|--builddir)   BUILDDIR="$2"; shift 2 ;;
    -j|--jobs)       JOBS="$2"; shift 2 ;;
    --make-arg)      MAKE_ARGS+=("$2"); shift 2 ;;
    --clean)         DO_CLEAN=1; shift ;;
    --no-build)      NO_BUILD=1; shift ;;
    --dry-run)       DRY_RUN=1; shift ;;
    -D|--dims)       DIMS="$2"; shift 2 ;;
    -n|--np)         NP="$2"; shift 2 ;;
    --threads)       OMP_THREADS="$2"; shift 2 ;;
    --tests)         TEST_OVERRIDE="$2"; shift 2 ;;
    --timeout)       TIMEOUT="$2"; shift 2 ;;
    --serial)        MPI_RUNNER=""; shift ;;
    --mpirun)        MPI_RUNNER="$2"; shift 2 ;;
    --mpi-arg)       MPI_ARGS+=("$2"); shift 2 ;;
    --fail-fast)     FAIL_FAST=1; shift ;;
    --list)          LIST_ONLY=1; shift ;;
    -h|--help)       usage; exit 0 ;;
    --)              shift; TEST_ARGS=("$@"); break ;;
    -*)              echo "error: unknown option '$1'" >&2; usage >&2; exit 2 ;;
    *)               echo "error: unexpected argument '$1'" >&2; usage >&2; exit 2 ;;
  esac
done

# ---------------------------------------------------------------------------
# Validate options
# ---------------------------------------------------------------------------
case "$OFFLOAD" in
  acc|omp|ompcpu|none) ;;
  *) echo "error: --offload must be one of acc, omp, ompcpu, none" >&2; exit 2 ;;
esac
case "$FLOAT_BITS" in
  32|64) ;;
  *) echo "error: --float-bits must be 32 or 64" >&2; exit 2 ;;
esac
case "$DIMS" in
  2|3|both) ;;
  *) echo "error: --dims must be 2, 3 or both" >&2; exit 2 ;;
esac
if ! [[ "$NP" =~ ^[0-9]+$ ]] || [[ "$NP" -lt 1 ]]; then
  echo "error: --np must be a positive integer" >&2; exit 2
fi
if ! [[ "$TIMEOUT" =~ ^[0-9]+$ ]]; then
  echo "error: --timeout must be a non-negative integer" >&2; exit 2
fi

# ---------------------------------------------------------------------------
# Prerequisite check
# ---------------------------------------------------------------------------
missing_tools=()
need() { command -v "$1" >/dev/null 2>&1 || missing_tools+=("$1"); }

if [[ "$LIST_ONLY" -eq 0 ]]; then
  need make
  need "$F90C"
fi
if [[ "$NO_BUILD" -eq 0 && "$DRY_RUN" -eq 0 ]]; then
  need mpicc
  need fypp
fi
if [[ "$LIST_ONLY" -eq 0 && "$DRY_RUN" -eq 0 && -n "$MPI_RUNNER" ]]; then
  need "$MPI_RUNNER"
fi

if [[ "${#missing_tools[@]}" -gt 0 ]]; then
  echo "error: required tool(s) not found: ${missing_tools[*]}" >&2
  for m in "${missing_tools[@]}"; do
    case "$m" in
      fypp)    echo "  - fypp: install via 'pip install --user fypp' (or conda)" >&2 ;;
      mpicc|mpif90) echo "  - $m: part of an MPI toolchain (OpenMPI/MPICH)" >&2 ;;
      mpirun)  echo "  - mpirun: required to run tests; pass --serial to run without it" >&2 ;;
    esac
  done
  exit 2
fi

if [[ "$NO_BUILD" -eq 0 && "$DRY_RUN" -eq 0 && ! -f p4est/build/local/lib/libp4est.a ]]; then
  echo "warning: p4est library not found under p4est/build/local/lib;" >&2
  echo "         build it first with 'git submodule update --init --recursive && bash build_p4est.sh'" >&2
fi

# ---------------------------------------------------------------------------
# Assemble the make variables and target
# ---------------------------------------------------------------------------
MAKE_VARS=()
[[ -n "$BUILDDIR" ]] && MAKE_VARS+=("BUILDDIR=$BUILDDIR")
MAKE_VARS+=("F90C=$F90C")
MAKE_VARS+=("OFFLOAD=$OFFLOAD")
MAKE_VARS+=("FLOAT_BITS=$FLOAT_BITS")
[[ "$DEBUG" -eq 1 ]] && MAKE_VARS+=("DEBUG=1")
[[ "$SAFE"  -eq 1 ]] && MAKE_VARS+=("SAFE=1")
if [[ "${#MAKE_ARGS[@]}" -gt 0 ]]; then
  MAKE_VARS+=("${MAKE_ARGS[@]}")
fi

MAKE_TARGET="all"
case "$DIMS" in
  2)    MAKE_TARGET="2d" ;;
  3)    MAKE_TARGET="3d" ;;
  both) MAKE_TARGET="all" ;;
esac

# ---------------------------------------------------------------------------
# Determine the build directory (ask the Makefile so logic stays single-source)
# ---------------------------------------------------------------------------
resolve_builddir() {
  if [[ -n "$BUILDDIR" ]]; then
    printf '%s\n' "$BUILDDIR"
    return
  fi
  local out
  out="$(make --no-print-directory "${MAKE_VARS[@]}" build-summary 2>/dev/null \
         | sed -n 's/.*Build directory *: *//p')"
  if [[ -z "$out" ]]; then
    echo "error: could not determine BUILDDIR from the Makefile" >&2
    return 1
  fi
  printf '%s\n' "$out"
}

# ---------------------------------------------------------------------------
# Compile the test list
# ---------------------------------------------------------------------------
build_test_list() {
  local -a tests=() sel=()
  if [[ -n "$TEST_OVERRIDE" ]]; then
    read -r -a tests <<< "$TEST_OVERRIDE"
  else
    if [[ "$DIMS" == "2" || "$DIMS" == "both" ]]; then
      read -r -a sel <<< "$MAIN_2D"; tests+=("${sel[@]}")
    fi
    if [[ "$DIMS" == "3" || "$DIMS" == "both" ]]; then
      read -r -a sel <<< "$MAIN_3D"; tests+=("${sel[@]}")
    fi
  fi
  printf '%s\n' "${tests[@]}"
}

TEST_LIST=()
while IFS= read -r line; do
  [[ -n "$line" ]] && TEST_LIST+=("$line")
done < <(build_test_list)

# ---------------------------------------------------------------------------
# --list
# ---------------------------------------------------------------------------
if [[ "$LIST_ONLY" -eq 1 ]]; then
  echo "Build options: F90C=$F90C OFFLOAD=$OFFLOAD FLOAT_BITS=$FLOAT_BITS DEBUG=$DEBUG SAFE=$SAFE"
  echo "Tests (DIMS=$DIMS, nps=$NP):"
  for t in "${TEST_LIST[@]}"; do echo "  $t"; done
  exit 0
fi

# ---------------------------------------------------------------------------
# Phase 1: build
# ---------------------------------------------------------------------------
BD=""
if [[ "$NO_BUILD" -eq 0 ]]; then
  echo "================================================================"
  echo " Building: compiler=$F90C  offload=$OFFLOAD  float=$FLOAT_BITS" \
       "debug=$DEBUG safe=$SAFE"
  echo " make target: $MAKE_TARGET   vars: ${MAKE_VARS[*]}"
  echo "================================================================"

  if [[ "$DO_CLEAN" -eq 1 ]]; then
    echo ">> cleaning build directory"
    make --no-print-directory clean "${MAKE_VARS[@]}" || true
  fi

  MAKE_CMD=(make --no-print-directory -j"$JOBS")
  [[ "$DRY_RUN" -eq 1 ]] && MAKE_CMD+=(-n)
  MAKE_CMD+=("$MAKE_TARGET" "${MAKE_VARS[@]}")

  if [[ "$DRY_RUN" -eq 1 ]]; then
    printf '+'; printf ' %q' "${MAKE_CMD[@]}"; printf '\n'
  fi

  if ! "${MAKE_CMD[@]}"; then
    echo "error: build failed" >&2
    exit 1
  fi
fi

# ---------------------------------------------------------------------------
# Locate binaries
# ---------------------------------------------------------------------------
if ! BD="$(resolve_builddir)"; then
  exit 1
fi
BINDIR="$BD/bin"
echo "Build directory : $BD"
echo "Binaries in     : $BINDIR"

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo ">> dry run: tests will not be executed"
  exit 0
fi

mkdir -p "$REPO_ROOT/output"
LOGDIR="$BD/test-logs"
mkdir -p "$LOGDIR"

if [[ -n "$OMP_THREADS" ]]; then
  export OMP_NUM_THREADS="$OMP_THREADS"
fi

# ---------------------------------------------------------------------------
# Phase 2: run tests
# ---------------------------------------------------------------------------
echo
echo "================================================================"
echo " Running tests (np=$NP, timeout=${TIMEOUT}s, logs in $LOGDIR)"
echo "================================================================"

declare -a PASSED=() FAILED=() MISSING=()
rc_overall=0

for t in "${TEST_LIST[@]}"; do
  bin="$BINDIR/$t"
  log="$LOGDIR/$t.log"

  if [[ ! -x "$bin" ]]; then
    printf 'SKIP  %-32s (binary not found: %s)\n' "$t" "$bin"
    MISSING+=("$t")
    rc_overall=1
    continue
  fi

  # The XDMF writer tests assert a single MPI rank.
  local_np="$NP"
  if [[ "$t" == test_xdmf_writer_* && "$NP" -ne 1 ]]; then
    echo "note: $t requires a single rank; using np=1" >&2
    local_np=1
  fi

  cmd=()
  [[ -n "$MPI_RUNNER" ]] && cmd+=("$MPI_RUNNER")
  [[ ${#MPI_ARGS[@]} -gt 0 ]] && cmd+=("${MPI_ARGS[@]}")
  [[ -n "$MPI_RUNNER" ]] && cmd+=(-np "$local_np")
  # OpenMP runs generally want the threads free to migrate across cores.
  if [[ -n "$MPI_RUNNER" && "$OFFLOAD" == "ompcpu" ]]; then
    cmd+=(--bind-to none)
  fi
  cmd+=("$bin")
  [[ ${#TEST_ARGS[@]} -gt 0 ]] && cmd+=("${TEST_ARGS[@]}")

  start="$(date +%s.%N)"
  if [[ "$TIMEOUT" -gt 0 ]]; then
    timeout "$TIMEOUT" "${cmd[@]}" >"$log" 2>&1
  else
    "${cmd[@]}" >"$log" 2>&1
  fi
  rc=$?
  end="$(date +%s.%N)"
  dur="$(awk -v a="$start" -v b="$end" 'BEGIN{printf "%.1f", b-a}')"

  if [[ "$rc" -eq 0 ]]; then
    printf 'PASS  %-32s %6ss\n' "$t" "$dur"
    PASSED+=("$t")
  else
    if [[ "$rc" -eq 124 ]]; then
      printf 'FAIL  %-32s %6ss  (timeout)\n' "$t" "$dur"
    else
      printf 'FAIL  %-32s %6ss  (exit %d)\n' "$t" "$dur" "$rc"
    fi
    FAILED+=("$t")
    rc_overall=1
    echo "      see $log"
    if [[ "$FAIL_FAST" -eq 1 ]]; then
      echo "      --fail-fast: stopping"
      break
    fi
  fi
done

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo
echo "================================================================"
echo " Summary"
echo "================================================================"
printf ' passed : %d\n' "${#PASSED[@]}"
printf ' failed : %d\n' "${#FAILED[@]}"
printf ' skipped: %d\n' "${#MISSING[@]}"
if [[ "${#FAILED[@]}" -gt 0 ]]; then
  echo " failed tests:"
  for t in "${FAILED[@]}"; do echo "   - $t"; done
fi
if [[ "${#MISSING[@]}" -gt 0 ]]; then
  echo " missing binaries:"
  for t in "${MISSING[@]}"; do echo "   - $t"; done
fi

exit "$rc_overall"
