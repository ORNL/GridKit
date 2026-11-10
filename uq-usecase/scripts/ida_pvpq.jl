#!/usr/bin/env julia
# ida_pvpq.jl: AC power flow with PV->PQ generator Q-limit switching.
#
# Separate from pm_solve.jl (the production driver, not modified). Usage:
#
#     julia --project=pm-solver ida_pvpq.jl manifest.json results.json
#
# manifest.json: [{"tag", "input", "flat", "start"}]. results.json: one record per case.
#
# Per case: solve_ac_pf. If a non-slack generator's Qg is outside [qmin, qmax],
# switch it to PQ with Qg fixed at the violated limit, then re-solve from the
# previous solution. Repeat until no violation or MAXIT rounds.
#
# PowerModels' PF formulation requires bus_type 2 on any bus with an active
# generator, so a switched generator is deactivated (gen_status 0). Its fixed
# injection (Pg, Qlim) is added as a load of (-Pg, -Qlim) at the same bus, and the
# bus becomes bus_type 1 only when no active generator remains on it.

#
# Switching is one-way (no PQ->PV reversal). The slack bus is never switched.
using PowerModels
using Ipopt
using JuMP
using JSON

const OK_STATUS = ("LOCALLY_SOLVED", "OPTIMAL", "ALMOST_LOCALLY_SOLVED")
const MAXIT = 30
const QTOL = 1e-6

function solve_pf(data)
    opt = optimizer_with_attributes(Ipopt.Optimizer,
        "print_level" => 0, "tol" => 1e-8, "max_iter" => 300)
    return solve_ac_pf(data, opt)
end

# Non-slack generators whose Qg is outside their limits: (gen key, bus key, violated limit).
function q_violations(data, res)
    out = Tuple{String,String,Float64}[]
    for (gk, g) in data["gen"]
        g["gen_status"] == 1 || continue
        bk = string(g["gen_bus"])
        data["bus"][bk]["bus_type"] == 3 && continue
        qg = res["solution"]["gen"][gk]["qg"]
        if qg > g["qmax"] + QTOL
            push!(out, (gk, bk, Float64(g["qmax"])))
        elseif qg < g["qmin"] - QTOL
            push!(out, (gk, bk, Float64(g["qmin"])))
        end
    end
    return out
end

function active_gens_at(data, bk)
    return [gk for (gk, g) in data["gen"] if g["gen_status"] == 1 && string(g["gen_bus"]) == bk]
end

function switch!(data, viols, switched)
    next_id = maximum(parse.(Int, collect(keys(data["load"])))) + 1
    touched = Set{String}()
    for (gk, bk, lim) in viols
        g = data["gen"][gk]
        g["gen_status"] = 0
        lid = next_id
        next_id += 1
        data["load"][string(lid)] = Dict{String,Any}(
            "index" => lid, "source_id" => Any["load", lid],
            "load_bus" => g["gen_bus"], "pd" => -Float64(g["pg"]), "qd" => -lim,
            "status" => 1)
        push!(switched, (gk, bk))
        push!(touched, bk)
    end
    for bk in touched
        isempty(active_gens_at(data, bk)) && (data["bus"][bk]["bus_type"] = 1)
    end
end

function copy_solution!(data, res)
    for (bk, b) in res["solution"]["bus"]
        data["bus"][bk]["vm"] = b["vm"]
        data["bus"][bk]["va"] = b["va"]
    end
end

function run_case(input, flat)
    data = PowerModels.parse_file(input)
    if flat
        for (_, b) in data["bus"]
            b["vm"] = 1.0
            b["va"] = 0.0
        end
    end
    switched = Tuple{String,String}[]
    plain_ok = false
    plain_qviol = -1
    res = nothing
    status = ""
    final_viols = Tuple{String,String,Float64}[]
    it_done = 0
    for it in 0:MAXIT
        res = solve_pf(data)
        status = string(res["termination_status"])
        if !(status in OK_STATUS)
            return Dict{String,Any}("plain_converged" => it == 0 ? false : plain_ok,
                "plain_qviol" => plain_qviol, "converged" => false, "status" => status,
                "n_switched" => length(switched), "iterations" => it, "final_qviol" => -1,
                "vmin" => nothing, "vmax" => nothing, "switched" => switched)
        end
        copy_solution!(data, res)
        viols = q_violations(data, res)
        if it == 0
            plain_ok = true
            plain_qviol = length(viols)
        end
        it_done = it
        final_viols = viols
        isempty(viols) && break
        it == MAXIT && break
        switch!(data, viols, switched)
    end
    vms = [b["vm"] for b in values(res["solution"]["bus"])]
    return Dict{String,Any}("plain_converged" => plain_ok, "plain_qviol" => plain_qviol,
        "converged" => true, "status" => status, "n_switched" => length(switched),
        "iterations" => it_done, "final_qviol" => length(final_viols),
        "vmin" => minimum(vms), "vmax" => maximum(vms), "switched" => switched)
end

function main(args)
    manifest = JSON.parsefile(args[1])
    out = Any[]
    for item in manifest
        rec = run_case(item["input"], item["flat"])
        rec["tag"] = item["tag"]
        rec["start"] = item["start"]
        push!(out, rec)
        println(stderr, "[ida_pvpq] $(item["tag"]): converged=$(rec["converged"]) switched=$(rec["n_switched"]) final_qviol=$(rec["final_qviol"])")
    end
    open(args[2], "w") do io
        JSON.print(io, out, 2)
    end
end

main(ARGS)
