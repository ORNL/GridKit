#!/usr/bin/env julia
"""
pm_solve.jl
-----------
PowerModels.jl-based AC power flow driver for MATPOWER .m files.

Mirrors the CLI contract of the GridKit `solve_pf` binary (see
uq-usecase/pf-solver/solve_pf.cpp), so the same Python parsing utilities
(pf_utils.py: `_parse_pf_stdout`, `diff_vs_base`) can consume its stdout,
and the same downstream tooling (m_to_case.py) can consume its solved .m
output.

Usage:
    julia112 --project=pm-solver pm_solve.jl <input.m> [--output-m <solved.m>] [--solver ipopt] [--pf-type ac|dc]

Behavior:
    - Parses <input.m> via PowerModels.parse_file
    - Solves power flow via PowerModels.solve_ac_pf (default) or solve_dc_pf
      (--pf-type dc). AC: fixed dispatch at PV buses, slack absorbs mismatch,
      reactive dispatch solved at all online gens. DC: linearized real-power-only
      approximation; voltage magnitude is not modeled (fixed at 1.0 p.u.), only
      angle is solved.
    - Prints one line per bus to stdout, in the same format as solve_pf:
          bus <i>  V=<pu>  theta_deg=<deg>  type=<1|2|3|4>
      (sorted by bus number)
    - Prints diagnostics to stderr (parse counts, solver, termination status,
      solve time)
    - Exit code 0 if converged (termination_status is LOCALLY_SOLVED,
      OPTIMAL, or ALMOST_LOCALLY_SOLVED), 1 otherwise
    - If --output-m is given and the solve converges, writes a solved .m file:
      mpc.bus columns 8-9 (Vm, Va) and mpc.gen columns 2-3 (Pg, Qg, online
      gens only) are replaced with the PF solution; everything else (branch,
      gencost, genfuel, bus_name, comments, offline gen rows) is copied
      through unchanged. This is a superset of solve_pf's --output-m, which
      only patches mpc.bus (see pm_helper.md for why: m_to_case.py needs the
      solved reactive dispatch Qg to patch Genrou q0 correctly).

Solver selection:
    --solver ipopt (default). The flag exists so other NLP/PF backends can be
    added later without changing the CLI contract; only "ipopt" is
    implemented today.
"""

using PowerModels
using Ipopt
using JuMP
using JSON

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

function parse_args(argv)
    if isempty(argv)
        println(stderr, "usage: pm_solve.jl <input.m> [--output-m <solved.m>] [--solver ipopt] [--pf-type ac|dc|ac-opf|dc-opf] [--tol <float>] [--max-iter <int>]")
        exit(2)
    end
    input_m = argv[1]
    output_m = nothing
    solver = "ipopt"
    pf_type = "ac"
    tol = 1e-8
    max_iter = 300
    flat_start = false
    i = 2
    while i <= length(argv)
        if argv[i] == "--output-m" && i < length(argv)
            output_m = argv[i+1]
            i += 2
        elseif argv[i] == "--solver" && i < length(argv)
            solver = argv[i+1]
            i += 2
        elseif argv[i] == "--pf-type" && i < length(argv)
            pf_type = argv[i+1]
            if !(pf_type in ("ac", "dc", "ac-opf", "dc-opf"))
                println(stderr, "--pf-type must be 'ac', 'dc', 'ac-opf' or 'dc-opf'")
                exit(2)
            end
            i += 2
        elseif argv[i] == "--tol" && i < length(argv)
            tol = parse(Float64, argv[i+1])
            i += 2
        elseif argv[i] == "--max-iter" && i < length(argv)
            max_iter = parse(Int, argv[i+1])
            i += 2
        elseif argv[i] == "--flat-start"
            flat_start = true
            i += 1
        else
            println(stderr, "unrecognized argument: $(argv[i])")
            exit(2)
        end
    end
    return input_m, output_m, solver, pf_type, tol, max_iter, flat_start
end

function parse_batch_args(argv)
    manifest = nothing
    results = nothing
    solver = "ipopt"
    tol = 1e-8
    max_iter = 300
    flat_start = false
    i = 1
    while i <= length(argv)
        a = argv[i]
        if a == "--batch" && i < length(argv)
            manifest = argv[i+1]; i += 2
        elseif a == "--results" && i < length(argv)
            results = argv[i+1]; i += 2
        elseif a == "--solver" && i < length(argv)
            solver = argv[i+1]; i += 2
        elseif a == "--tol" && i < length(argv)
            tol = parse(Float64, argv[i+1]); i += 2
        elseif a == "--max-iter" && i < length(argv)
            max_iter = parse(Int, argv[i+1]); i += 2
        elseif a == "--flat-start"
            flat_start = true; i += 1
        else
            println(stderr, "unrecognized batch argument: $(a)")
            exit(2)
        end
    end
    if manifest === nothing
        println(stderr, "usage: pm_solve.jl --batch <manifest.json> --results <results.json>")
        exit(2)
    end
    if results === nothing
        results = manifest * ".results.json"
    end
    return manifest, results, solver, tol, max_iter, flat_start
end

function make_optimizer(solver::String, tol::Float64, max_iter::Int)
    if solver == "ipopt"
        # print_level=3 gives a concise convergence summary (iteration count,
        # final NLP error) without the per-iteration table. stdout is redirected
        # to stderr during the solve so this output doesn't pollute the bus lines
        # that Python reads from stdout. sb=yes suppresses the Ipopt banner.
        return optimizer_with_attributes(
            Ipopt.Optimizer,
            "print_level" => 3,
            "sb" => "yes",
            "tol" => tol,
            "max_iter" => max_iter,
        )
    else
        error("pm_solve.jl: unsupported --solver '$(solver)' (only 'ipopt' is implemented)")
    end
end

# ---------------------------------------------------------------------------
# Solved .m writer (minimal diff: patch mpc.bus Vm/Va + mpc.gen Pg/Qg only)
# ---------------------------------------------------------------------------

"""
    write_solved_m(src_path, dst_path, bus_solution, gen_solution)

`bus_solution`: Dict{Int,Tuple{Float64,Float64}}  bus_i => (Vm_pu, Va_deg)
`gen_solution`: Dict{Int,Tuple{Float64,Float64}}  mpc.gen row (1-based) => (Pg_MW, Qg_MVAr)

Line-based patch: detects `mpc.bus = [ ... ];` and `mpc.gen = [ ... ];` blocks
and rewrites only the relevant numeric columns of each row, preserving the
original whitespace/formatting of every other line (comments, other blocks,
non-matched rows) byte-for-byte.
"""
function write_solved_m(
    src_path::String,
    dst_path::String,
    bus_solution::Dict{Int,Tuple{Float64,Float64}},
    gen_solution::Dict{Int,Tuple{Float64,Float64}},
)
    lines = readlines(src_path; keep=true)
    out = IOBuffer()

    in_bus = false
    in_gen = false
    gen_row = 0

    for line in lines
        stripped = strip(split(line, '%')[1])

        if !in_bus && !in_gen
            if occursin(r"^\s*mpc\.bus\s*=\s*\[", line)
                in_bus = true
                write(out, line)
                continue
            elseif occursin(r"^\s*mpc\.gen\s*=\s*\[", line)
                in_gen = true
                gen_row = 0
                write(out, line)
                continue
            end
            write(out, line)
            continue
        end

        if in_bus
            if occursin("];", stripped) || stripped == "]"
                in_bus = false
                write(out, line)
                continue
            end
            if isempty(stripped)
                write(out, line)
                continue
            end
            tokens = split(stripped, r"\s+")
            bus_i = tryparse(Int, tokens[1])
            if bus_i === nothing || !haskey(bus_solution, bus_i)
                write(out, line)
                continue
            end
            vm, va_deg = bus_solution[bus_i]
            tokens[8] = string(round(vm; digits=7))
            tokens[9] = string(round(va_deg; digits=6))
            write(out, "\t" * join(tokens, "\t") * "\n")
            continue
        end

        if in_gen
            if occursin("];", stripped) || stripped == "]"
                in_gen = false
                write(out, line)
                continue
            end
            if isempty(stripped)
                write(out, line)
                continue
            end
            gen_row += 1
            tokens = split(stripped, r"\s+")
            if haskey(gen_solution, gen_row)
                pg, qg = gen_solution[gen_row]
                tokens[2] = string(round(pg; digits=6))
                tokens[3] = string(round(qg; digits=6))
            end
            write(out, "\t" * join(tokens, "\t") * ";\n")
            continue
        end
    end

    open(dst_path, "w") do io
        write(io, String(take!(out)))
    end
end

# ---------------------------------------------------------------------------
# Solve one case
# ---------------------------------------------------------------------------

"""
    solve_one(input_m, output_m, solver, pf_type, tol, max_iter, flat_start; emit_bus_lines)

Solve a single case. Returns a Dict describing the outcome; it does NOT exit,
so a batch can continue after a failure.

`emit_bus_lines` prints the one-line-per-bus block on stdout that the single
case CLI contract promises. Batch mode turns it off: with many cases the lines
are ambiguous, and the caller reads the solved `.m` instead.
"""
function solve_one(input_m, output_m, solver, pf_type, tol, max_iter, flat_start;
                   emit_bus_lines::Bool=true)
    data = PowerModels.parse_file(input_m)

    n_bus = length(data["bus"])
    n_gen = length(data["gen"])
    n_off = count(g -> g["gen_status"] == 0, values(data["gen"]))
    baseMVA = data["baseMVA"]

    if flat_start
        for (_, b) in data["bus"]
            b["vm"] = 1.0
            b["va"] = 0.0
        end
    end

    println(stderr, "[pm_solve] parsed $(input_m): $(n_bus) buses, $(n_gen) gens ($(n_off) offline), baseMVA=$(baseMVA)")
    println(stderr, "[pm_solve] solver=$(solver)  pf_type=$(pf_type)  tol=$(tol)  max_iter=$(max_iter)  flat_start=$(flat_start)")

    optimizer = make_optimizer(solver, tol, max_iter)

    # Redirect stdout -> stderr during solve so Ipopt's print_level=3 summary
    # (iteration count, final NLP error) goes to stderr, not stdout.
    local result
    let prev_stdout = stdout
        redirect_stdout(stderr)
        try
            result = if pf_type == "dc"
                solve_dc_pf(data, optimizer)
            elseif pf_type == "ac-opf"
                solve_ac_opf(data, optimizer)
            elseif pf_type == "dc-opf"
                solve_dc_opf(data, optimizer)
            else
                solve_ac_pf(data, optimizer)
            end
        finally
            redirect_stdout(prev_stdout)
        end
    end

    status = result["termination_status"]
    converged = status in (
        JuMP.LOCALLY_SOLVED,
        JuMP.OPTIMAL,
        JuMP.ALMOST_LOCALLY_SOLVED,
        JuMP.ALMOST_OPTIMAL,
    )

    println(stderr, "[pm_solve] termination_status=$(status)")
    println(stderr, "[pm_solve] solve_time=$(get(result, "solve_time", NaN))s")
    println(stderr, "[pm_solve] converged=$(converged)")
    # OPF minimises generation cost against mpc.gencost; PF has no objective.
    if endswith(pf_type, "opf")
        println(stderr, "[pm_solve] objective=$(get(result, "objective", NaN))")
    end

    out = Dict{String,Any}(
        "input" => input_m,
        "output" => output_m === nothing ? "" : output_m,
        "pf_type" => pf_type,
        "ok" => converged,
        "status" => string(status),
        "solve_time" => get(result, "solve_time", NaN),
        "objective" => endswith(pf_type, "opf") ? get(result, "objective", NaN) : NaN,
        "wrote_output" => false,
    )
    if !converged
        return out
    end

    # --- bus solution: bus_i => (Vm, Va_deg) ---
    # Note: for pf_type="dc", vm is fixed at 1.0 p.u. by the DC approximation
    # (voltage magnitude is not modeled); only va is meaningful.
    bus_solution = Dict{Int,Tuple{Float64,Float64}}()
    for (k, sol) in result["solution"]["bus"]
        bus_i = data["bus"][k]["bus_i"]
        vm = sol["vm"]
        va_deg = sol["va"] * 180.0 / pi
        bus_solution[bus_i] = (vm, va_deg)
    end

    # --- gen solution: mpc.gen row (1-based) => (Pg_MW, Qg_MVAr) ---
    # DC has no reactive dispatch, so qg is absent from its solution. AC PF
    # solves qg at fixed pg; AC OPF optimises both.
    gen_solution = Dict{Int,Tuple{Float64,Float64}}()
    if pf_type in ("ac", "ac-opf", "dc-opf")
        for (k, sol) in result["solution"]["gen"]
            row = data["gen"][k]["source_id"][2]
            pg_mw = get(sol, "pg", NaN) * baseMVA
            qg_mvar = get(sol, "qg", 0.0) * baseMVA
            gen_solution[row] = (pg_mw, qg_mvar)
        end
    end

    if emit_bus_lines
        for bus_i in sort(collect(keys(bus_solution)))
            vm, va_deg = bus_solution[bus_i]
            btype = 0
            for (k, b) in data["bus"]
                if b["bus_i"] == bus_i
                    btype = b["bus_type"]
                    break
                end
            end
            println("bus $(bus_i)  V=$(vm)  theta_deg=$(va_deg)  type=$(btype)")
        end
    end

    if output_m !== nothing && output_m != ""
        if pf_type == "dc"
            println(stderr, "[pm_solve] WARNING: --output-m with pf_type=dc will write vm=1.0 for all buses (DC approximation); use pf_type=ac to produce a physically meaningful solved .m")
        end
        write_solved_m(input_m, output_m, bus_solution, gen_solution)
        println(stderr, "[pm_solve] wrote solved .m -> $(output_m)")
        out["wrote_output"] = true
    end

    return out
end

# ---------------------------------------------------------------------------
# Batch
# ---------------------------------------------------------------------------

"""
    run_batch(manifest_path, results_path, solver, tol, max_iter, flat_start)

Solve many cases in one Julia process, which is the point: interpreter startup
dominates a single-case run (~27 s of a 30 s condition), while a warm solve is
0.03-3 s.

**A failure in one case never aborts the batch and never stops at that case.**
Each entry gets its own `ok` flag in the results file so the caller can skip
the dynamic simulation for exactly the cases that failed. The process exits 0
whenever the manifest itself was processed; a nonzero exit means the batch
could not be run at all, not that some case failed to converge.
"""
function run_batch(manifest_path, results_path, solver, tol, max_iter, flat_start)
    manifest = JSON.parsefile(manifest_path)
    results = Vector{Dict{String,Any}}()
    n_ok = 0
    for (i, item) in enumerate(manifest)
        inp = item["input"]
        outp = get(item, "output", "")
        pft = get(item, "pf_type", "ac")
        println(stderr, "[pm_solve] === batch $(i)/$(length(manifest)): $(inp) ($(pft)) ===")
        rec = try
            if !(pft in ("ac", "dc", "ac-opf", "dc-opf"))
                error("invalid pf_type '$(pft)'")
            end
            solve_one(inp, outp, solver, pft, tol, max_iter, flat_start;
                      emit_bus_lines=false)
        catch e
            # Parse errors, missing files, solver blowups: record and continue.
            msg = sprint(showerror, e)
            println(stderr, "[pm_solve] case FAILED with exception: $(msg)")
            Dict{String,Any}(
                "input" => inp, "output" => outp, "pf_type" => pft,
                "ok" => false, "status" => "EXCEPTION", "error" => msg,
                "solve_time" => NaN, "objective" => NaN, "wrote_output" => false,
            )
        end
        rec["index"] = i - 1
        n_ok += rec["ok"] ? 1 : 0
        push!(results, rec)
    end
    # JSON has no NaN literal; JSON.jl refuses to write one. Non-finite numbers
    # here only ever mean "not applicable" (no objective for a PF, no solve_time
    # for a case that threw), so null is the honest encoding.
    for rec in results, (k, v) in rec
        if v isa AbstractFloat && !isfinite(v)
            rec[k] = nothing
        end
    end
    open(results_path, "w") do io
        JSON.print(io, results)
    end
    println(stderr, "[pm_solve] batch done: $(n_ok)/$(length(manifest)) converged -> $(results_path)")
    return 0
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    if "--batch" in ARGS
        manifest, results, solver, tol, max_iter, flat_start = parse_batch_args(ARGS)
        PowerModels.silence()
        exit(run_batch(manifest, results, solver, tol, max_iter, flat_start))
    end

    input_m, output_m, solver, pf_type, tol, max_iter, flat_start = parse_args(ARGS)
    PowerModels.silence()
    out = solve_one(input_m, output_m, solver, pf_type, tol, max_iter, flat_start)
    exit(out["ok"] ? 0 : 1)
end

main()
