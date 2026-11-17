#!/usr/bin/env sh

set -eu
readonly _cosim_dir=$(dirname "$(realpath "$0")")
cd "${_cosim_dir}"

case "${1:-three-bus}" in
    three-bus)
        _server=ThreeBusCoSimServer.case.json
        _client=ThreeBusCoSimClient.case.json
        _server_bus=1 _client_bus=2
        _vr=1 _vi=2 _ir=3 _ii=4
        _srv_vr=3 _srv_vi=4 _srv_ir=1 _srv_ii=2
        _dt=0.0010416666666666667
        ;;
    ieee39)
        _server=IEEE39Generator.case.json
        _client=IEEE39Grid.case.json
        _server_bus=39 _client_bus=39
        _vr=30 _vi=31 _ir=32 _ii=33
        _srv_vr=30 _srv_vi=31 _srv_ir=32 _srv_ii=33
        # determined through experimentation as a timestep with small amounts
        # of error but there are more options. time steps smaller than about
        # 7.93ms appear to fail uniformly. those above intermittently fail
        # but are more likely to pass. time steps near 7.93ms fail with
        # large amounts of error.
        _dt=0.125
        ;;
    *) echo "Usage: $0 [three-bus|ieee39] [bus|signals]" >&2; exit 2 ;;
esac

case "${2:-bus}" in
    bus) set -- --output-bus "$_server_bus" ;;
    signals) set -- --output-signals "$_srv_ir" "$_srv_ii" --input-signals "$_srv_vr" "$_srv_vi" ;;
    *) echo "Output selection must be bus or signals" >&2; exit 2 ;;
esac
_output_mode=$1
./CoSimServer "$_server" "$@" --input-tags "${_client_bus}_vr" "${_client_bus}_vi" --output-tags "${_server_bus}_ir" "${_server_bus}_ii" --time-step "$_dt" &
_srv_pid=$!
trap 'kill "$_srv_pid" 2>/dev/null || true' EXIT HUP INT TERM

if [ "$_output_mode" = --output-bus ]; then
    set -- --output-bus "$_client_bus"
else
    set -- --output-signals "$_vr" "$_vi" --input-signals "$_ir" "$_ii"
fi
./CoSimClient "$_client" "$@" --input-tags "${_server_bus}_ir" "${_server_bus}_ii" --output-tags "${_client_bus}_vr" "${_client_bus}_vi" --time-step "$_dt"
wait "$_srv_pid"
trap - EXIT HUP INT TERM
