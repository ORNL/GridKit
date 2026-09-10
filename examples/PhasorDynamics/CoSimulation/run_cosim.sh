#!/bin/sh

set -eu

example_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
dynamic=DynamicSimulation
analysis=
case_name=three-bus
interval=
endpoint=
client_study=
server_study=

usage() {
    cat <<EOF
Usage: sh run_cosim.sh [options]
  --case two-bus|three-bus|ieee39|ieee39-contingencies
  --dynamic PATH             DynamicSimulation executable
  --analysis PATH            ContingencyAnalysis executable (if specified, run
                             contingencies instead)
  --dt SECONDS               Override the exchange interval
  --endpoint ADDRESS         Override the temporary IPC endpoint
  --client-study PATH        Override the client study
  --server-study PATH        Override the server study
EOF
}

while [ "$#" -gt 0 ]; do
    case "$1" in
        --help|-h) usage; exit 0 ;;
        --case) case_name=${2-} ;;
        --dynamic) dynamic=${2-} ;;
        --analysis) analysis=${2-} ;;
        --dt) interval=${2-} ;;
        --endpoint) endpoint=${2-} ;;
        --client-study) client_study=${2-} ;;
        --server-study) server_study=${2-} ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
    if [ "$#" -lt 2 ]; then
        echo "Missing value for $1" >&2
        exit 1
    fi
    shift 2
done

case "$case_name" in
    two-bus|three-bus)
        case "$case_name" in
          two-bus) prefix=TwoBus ;;
          three-bus) prefix=ThreeBus ;;
        esac

        toy_dir=$example_dir/Toy
        client_study=${client_study:-$toy_dir/${prefix}CoSimClient.solver.json}
        server_study=${server_study:-$toy_dir/${prefix}CoSimServer.solver.json}
        interval=${interval:-0.0010416666666666667}
        client_vr=1; client_vi=2; client_ir=3; client_ii=4
        server_vr=3; server_vi=4; server_ir=1; server_ii=2 ;;
    ieee39|ieee39-contingencies)
        suffix=solver.json
        if [ "$case_name" = ieee39-contingencies ]; then
            suffix=contingency.solver.json
            analysis=${analysis:-ContingencyAnalysis}
        fi

        ieee39_dir=$example_dir/IEEE39
        client_study=${client_study:-$ieee39_dir/IEEE39Grid.$suffix}
        server_study=${server_study:-$ieee39_dir/IEEE39Generator.$suffix}
        interval=${interval:-0.0005208333333333333}
        client_vr=30; client_vi=31; client_ir=32; client_ii=33
        server_vr=30; server_vi=31; server_ir=32; server_ii=33 ;;
    *) echo "Unknown case: $case_name" >&2; exit 1 ;;
esac

work_dir=$(mktemp -d "${TMPDIR:-/tmp}/gridkit-cosim.XXXXXXXX")
server_pid=
cleanup() {
    if [ -n "$server_pid" ]; then
        kill "$server_pid" 2>/dev/null || :
        wait "$server_pid" 2>/dev/null || :
    fi
    rm -rf -- "$work_dir"
}
trap cleanup 0
trap 'exit 1' HUP INT TERM
endpoint=${endpoint:-ipc://$work_dir/peer}

set -- "$server_study" --cosim-role server \
    --cosim-endpoint "$endpoint" --cosim-dt "$interval" \
    --cosim-input "vr=$server_vr" --cosim-input "vi=$server_vi" \
    --cosim-output "ir=$server_ir" --cosim-output "ii=$server_ii"
if [ -n "$analysis" ]; then
    set -- "$@" --cosim-serve
fi
"$dynamic" "$@" &
server_pid=$!

"${analysis:-$dynamic}" "$client_study" --cosim-role client \
    --cosim-endpoint "$endpoint" --cosim-dt "$interval" \
    --cosim-input "ir=$client_ir" --cosim-input "ii=$client_ii" \
    --cosim-output "vr=$client_vr" --cosim-output "vi=$client_vi"

wait "$server_pid"
server_pid=
