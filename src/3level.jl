using Dates
using Printf
using TOML
using JSON
using PowerModels
using LinearAlgebra
using DelimitedFiles

using JuMP
using Gurobi
const GRB_ENV = Gurobi.Env()

const battery_multiplier = Dict(
    "ieee_14" => 0.5,
    "ieee_73" => 1.0,
    "pegase_89" => 1.0,
    "ieee_118" => 2.0,
    "ieee_162" => 2.0,
    "ieee_300" => 6,
    "pegase_1354" => 20.0,
    "rte_1888" => 20.0,
)
const time_limit = 3600 * 6

include("network.jl")

# user input
network_name   = ARGS[1]
b              = parse(Int, ARGS[2])
k              = parse(Int, ARGS[3])
simulation     = parse(Int, ARGS[4])
method         = ARGS[5]

demand_path     = "./data/demand"
network_path    = "./data/network"
output_path     = "./"
if !isdir(output_path)
    mkdir(output_path)
end

# network data
network = generate_network_instance(joinpath(network_path, "$(network_name).m"))
T = 24  # number of timestamps
N = network.N
E = network.E
L = network.L
Pmax = network.Pmax
Pmin = network.Pmin
Fmax = network.Fmax
B    = network.B
δin  = network.δin
δout = network.δout

# demand scenarios
if b == 2 && k == 3
    demand = readdlm(joinpath(demand_path, "$(network_name)_demand_$(simulation).txt"))
elseif b == 2 && k == 5
    demand = readdlm(joinpath(demand_path, "$(network_name)_demand_$(simulation+10).txt"))
elseif b == 3 && k == 5
    demand = readdlm(joinpath(demand_path, "$(network_name)_demand_$(simulation+20).txt"))
elseif b == 5 && k == 10
    demand = readdlm(joinpath(demand_path, "$(network_name)_demand_$(simulation+30).txt"))
end

# battery parameters
battery_params = TOML.parsefile("./data/battery_params.toml")
multiplier = battery_multiplier[network_name]
PCmax = battery_params["PCmax"] * multiplier
PCmin = battery_params["PCmin"] * multiplier
PDmax = battery_params["PDmax"] * multiplier
PDmin = battery_params["PDmin"] * multiplier
Emax = battery_params["Emax"] * multiplier
Emin = battery_params["Emin"] * multiplier
E0 = battery_params["E0"] * multiplier
e = battery_params["e"]

function solve_upper_bound_vch(x;
        penalty::Real = 0.0,
        ϵ_gap::Float64=5e-3,
        verbose = true,
    )
    model = JuMP.direct_model(MOI.instantiate(() -> Gurobi.Optimizer(GRB_ENV)))
    @variable(model, y[1:E], Bin)
    @variable(model, zp[1:T, 1:E] ≥ 0)
    @variable(model, zn[1:T, 1:E] ≥ 0)
    @variable(model, α[1:T, 1:N])
    @variable(model, βp[1:T, 1:E] ≥ 0)
    @variable(model, βn[1:T, 1:E] ≥ 0)
    @variable(model, γp[1:T, 1:N] ≥ 0)
    @variable(model, γn[1:T, 1:N] ≥ 0)
    @variable(model, τ[1:T, 1:N])
    @variable(model, τ0[1:N])
    @variable(model, νp[1:T, 1:N] ≥ 0)
    @variable(model, νn[1:T, 1:N] ≥ 0)
    @variable(model, ωp[1:T, 1:N] ≥ 0)
    @variable(model, ωn[1:T, 1:N] ≥ 0)
    @variable(model, μp[1:T, 1:N] ≥ 0)
    @variable(model, μn[1:T, 1:N] ≥ 0)
    @variable(model, ϕ[1:T, 1:N] ≥ 0)

    # vch variables
    @variable(model, vch1[1:T, 1:N])
    @variable(model, vch2[1:T, 1:N])
    @variable(model, vch3[1:T, 1:N])
    @variable(model, vch4[1:T, 1:N])

    # dual constraints
    pex_dual  = @constraint(model, [t in 1:T, i in 1:N], α[t,i] ≤ 1)
    pls_dual  = @constraint(model, [t in 1:T, i in 1:N], - α[t,i] ≤ 1)
    f_dual    = @constraint(model, [e in 1:E, t in 1:T], α[t, L[e][1]] - α[t, L[e][2]] + βp[t,e] - βn[t,e] == 0)
    pg_dual   = @constraint(model, [i in 1:N, t in 1:T], - α[t, i] + γp[t, i] - γn[t, i] == 0)
    pc_dual   = @constraint(model, [i in 1:N, t in 1:T], α[t, i] - e * τ[t, i] + νp[t, i] - νn[t, i] + vch1[t,i] == penalty)
    pd_dual   = @constraint(model, [i in 1:N, t in 1:T], - α[t, i] + 1/e * τ[t, i] + ωp[t, i] - ωn[t, i] + vch2[t,i] == penalty)
    ps_dual_t = @constraint(model, [i in 1:N, t in 1:T-1], τ[t,i] - τ[t+1,i] + μp[t,i] - μn[t,i] + vch3[t,i] == 0)
    ps_dual_T = @constraint(model, [i in 1:N], τ[T,i] + μp[T,i] - μn[T,i] + vch3[T,i] == 0)
    ps_dual_0 = @constraint(model, [i in 1:N], τ0[i] - τ[1,i] == 0)
    u_dual    = @constraint(model, [i in 1:N, t in 1:T], -PCmin * νp[t,i] + PCmax * νn[t,i] + PDmin * ωp[t,i] - PDmax * ωn[t,i] - ϕ[t,i] ≤ 0)

    # dual of vch constraints
    hat_pc = [0, PCmax, PCmax, 0, 0, 0]
    hat_pd = [0, 0, 0, 0, PDmax, PDmax]
    hat_e  = [Emin, Emin, Emax - PCmax * e, Emax, Emax, Emin + PDmax / e]

    lambda_dual = @constraint(model, [i in 1:N, t in 1:T, k_ in 1:6], - hat_pc[k_] * vch1[t,i] - hat_pd[k_] * vch2[t,i] - hat_e[k_] * vch3[t,i] + vch4[t,i] ≤ 0)

    # N-K constraint
    @constraint(model, sum(y) ≤ k)

    # McCormick envelope constraints
    @constraint(model, [t in 1:T, e in 1:E], zp[t, e] ≥ βp[t, e] + 2 * y[e] - 2)
    @constraint(model, [t in 1:T, e in 1:E], zp[t, e] ≤ βp[t, e])
    @constraint(model, [t in 1:T, e in 1:E], zp[t, e] ≤ 2 * y[e])
    @constraint(model, [t in 1:T, e in 1:E], zn[t, e] ≥ βn[t, e] + 2 * y[e] - 2)
    @constraint(model, [t in 1:T, e in 1:E], zn[t, e] ≤ βn[t, e])
    @constraint(model, [t in 1:T, e in 1:E], zn[t, e] ≤ 2 * y[e])

    @objective(model, Max, (
        - sum(demand[t,i] * α[t,i] for t in 1:T, i in 1:N) 
        - sum(Fmax[e] * (βp[t,e] + βn[t,e]) for t in 1:T, e in 1:E)
        + sum(Fmax[e] * (zp[t,e] + zn[t,e]) for t in 1:T, e in 1:E)
        + sum(Pmin[i] * γp[t,i] for t in 1:T, i in 1:N)
        - sum(Pmax[i] * γn[t,i] for t in 1:T, i in 1:N)
        + sum(PDmin * x[i] * ωp[t,i] for t in 1:T, i in 1:N)
        - sum(PDmax * x[i] * ωn[t,i] for t in 1:T, i in 1:N)
        + sum(Emin * x[i] * μp[t,i] for t in 1:T, i in 1:N)
        - sum(Emax * x[i] * μn[t,i] for t in 1:T, i in 1:N)
        + sum(E0 * x[i] * τ0[i] for i in 1:N)
        - sum(x[i] * ϕ[t,i] for t in 1:T, i in 1:N)
        + sum(vch4[t,i] for t in 1:T, i in 1:N)
    ))
    set_optimizer_attribute(model, MOI.Silent(), !verbose)
    set_optimizer_attribute(model, "MIPGap", ϵ_gap)
    optimize!(model)

    return model
end

function solve_upper_bound(x;
        penalty::Real = 0.0,
        ϵ_gap::Float64=5e-3,
        verbose = true,
    )
    model = JuMP.direct_model(MOI.instantiate(() -> Gurobi.Optimizer(GRB_ENV)))
    @variable(model, y[1:E], Bin)
    @variable(model, zp[1:T, 1:E] ≥ 0)
    @variable(model, zn[1:T, 1:E] ≥ 0)
    @variable(model, α[1:T, 1:N])
    @variable(model, βp[1:T, 1:E] ≥ 0)
    @variable(model, βn[1:T, 1:E] ≥ 0)
    @variable(model, γp[1:T, 1:N] ≥ 0)
    @variable(model, γn[1:T, 1:N] ≥ 0)
    @variable(model, τ[1:T, 1:N])
    @variable(model, τ0[1:N])
    @variable(model, νp[1:T, 1:N] ≥ 0)
    @variable(model, νn[1:T, 1:N] ≥ 0)
    @variable(model, ωp[1:T, 1:N] ≥ 0)
    @variable(model, ωn[1:T, 1:N] ≥ 0)
    @variable(model, μp[1:T, 1:N] ≥ 0)
    @variable(model, μn[1:T, 1:N] ≥ 0)
    @variable(model, ϕ[1:T, 1:N] ≥ 0)

    # dual constraints
    pex_dual  = @constraint(model, [t in 1:T, i in 1:N], α[t,i] ≤ 1)
    pls_dual  = @constraint(model, [t in 1:T, i in 1:N], - α[t,i] ≤ 1)
    f_dual    = @constraint(model, [e in 1:E, t in 1:T], α[t, L[e][1]] - α[t, L[e][2]] + βp[t,e] - βn[t,e] == 0)
    pg_dual   = @constraint(model, [i in 1:N, t in 1:T], - α[t, i] + γp[t, i] - γn[t, i] == 0)
    pc_dual   = @constraint(model, [i in 1:N, t in 1:T], α[t, i] - e * τ[t, i] + νp[t, i] - νn[t, i] == penalty)
    pd_dual   = @constraint(model, [i in 1:N, t in 1:T], - α[t, i] + 1/e * τ[t, i] + ωp[t, i] - ωn[t, i] == penalty)
    ps_dual_t = @constraint(model, [i in 1:N, t in 1:T-1], τ[t,i] - τ[t+1,i] + μp[t,i] - μn[t,i] == 0)
    ps_dual_T = @constraint(model, [i in 1:N], τ[T,i] + μp[T,i] - μn[T,i] == 0)
    ps_dual_0 = @constraint(model, [i in 1:N], τ0[i] - τ[1,i] == 0)
    u_dual    = @constraint(model, [i in 1:N, t in 1:T], -PCmin * νp[t,i] + PCmax * νn[t,i] + PDmin * ωp[t,i] - PDmax * ωn[t,i] - ϕ[t,i] ≤ 0)

    # N-K constraint
    @constraint(model, sum(y) ≤ k)

    # McCormick envelope constraints
    @constraint(model, [t in 1:T, e in 1:E], zp[t, e] ≥ βp[t, e] + 2 * y[e] - 2)
    @constraint(model, [t in 1:T, e in 1:E], zp[t, e] ≤ βp[t, e])
    @constraint(model, [t in 1:T, e in 1:E], zp[t, e] ≤ 2 * y[e])
    @constraint(model, [t in 1:T, e in 1:E], zn[t, e] ≥ βn[t, e] + 2 * y[e] - 2)
    @constraint(model, [t in 1:T, e in 1:E], zn[t, e] ≤ βn[t, e])
    @constraint(model, [t in 1:T, e in 1:E], zn[t, e] ≤ 2 * y[e])

    @objective(model, Max, (
        - sum(demand[t,i] * α[t,i] for t in 1:T, i in 1:N) 
        - sum(Fmax[e] * (βp[t,e] + βn[t,e]) for t in 1:T, e in 1:E)
        + sum(Fmax[e] * (zp[t,e] + zn[t,e]) for t in 1:T, e in 1:E)
        + sum(Pmin[i] * γp[t,i] for t in 1:T, i in 1:N)
        - sum(Pmax[i] * γn[t,i] for t in 1:T, i in 1:N)
        + sum(PDmin * x[i] * ωp[t,i] for t in 1:T, i in 1:N)
        - sum(PDmax * x[i] * ωn[t,i] for t in 1:T, i in 1:N)
        + sum(Emin * x[i] * μp[t,i] for t in 1:T, i in 1:N)
        - sum(Emax * x[i] * μn[t,i] for t in 1:T, i in 1:N)
        + sum(E0 * x[i] * τ0[i] for i in 1:N)
        - sum(x[i] * ϕ[t,i] for t in 1:T, i in 1:N)
    ))
    set_optimizer_attribute(model, MOI.Silent(), !verbose)
    set_optimizer_attribute(model, "MIPGap", ϵ_gap)
    optimize!(model)

    return model
end

function solve_lower_bound(Y;
        verbose::Bool = false,
        ϵ_gap::Float64=5e-3,
    )
    S = length(Y)

    model = JuMP.direct_model(MOI.instantiate(() -> Gurobi.Optimizer(GRB_ENV)))
    @variable(model, x[1:N], Bin)
    @variable(model, v[1:S] .≥ 0)
    @variable(model, u)

    @objective(model, Min, u + sum(v))

    for (k, (y, zp, zn, α, βp, βn, γp, γn, τ, τ0, νp, νn, ωp, ωn, μp, μn, ϕ)) in enumerate(Y)
        @constraint(model, u + v[k] ≥ 
            - sum(demand[t,i] * α[t,i] for t in 1:T, i in 1:N) 
            - sum(Fmax[e] * (βp[t,e] + βn[t,e]) for t in 1:T, e in 1:E)
            + sum(Fmax[e] * (zp[t,e] + zn[t,e]) for t in 1:T, e in 1:E)
            + sum(Pmin[i] * γp[t,i] for t in 1:T, i in 1:N)
            - sum(Pmax[i] * γn[t,i] for t in 1:T, i in 1:N)
            + sum(PDmin * x[i] * ωp[t,i] for t in 1:T, i in 1:N)
            - sum(PDmax * x[i] * ωn[t,i] for t in 1:T, i in 1:N)
            + sum(Emin * x[i] * μp[t,i] for t in 1:T, i in 1:N)
            - sum(Emax * x[i] * μn[t,i] for t in 1:T, i in 1:N)
            + sum(E0 * x[i] * τ0[i] for i in 1:N)
            - sum(x[i] * ϕ[t,i] for t in 1:T, i in 1:N)
        )
    end
    @constraint(model, sum(x) ≤ b)

    set_optimizer_attribute(model, MOI.Silent(), !verbose)
    set_optimizer_attribute(model, "MIPGap", ϵ_gap)
    optimize!(model)
    return model
end

function solve_lower_bound_vch(Y;
        verbose::Bool = false,
        ϵ_gap::Float64=5e-3,
    )
    S = length(Y)

    model = JuMP.direct_model(MOI.instantiate(() -> Gurobi.Optimizer(GRB_ENV)))
    @variable(model, x[1:N], Bin)
    @variable(model, v[1:S] .≥ 0)
    @variable(model, u)

    @objective(model, Min, u + sum(v))

    for (k, (y, zp, zn, α, βp, βn, γp, γn, τ, τ0, νp, νn, ωp, ωn, μp, μn, ϕ, vch1, vch2, vch3, vch4)) in enumerate(Y)
        @constraint(model, u + v[k] ≥ 
            - sum(demand[t,i] * α[t,i] for t in 1:T, i in 1:N) 
            - sum(Fmax[e] * (βp[t,e] + βn[t,e]) for t in 1:T, e in 1:E)
            + sum(Fmax[e] * (zp[t,e] + zn[t,e]) for t in 1:T, e in 1:E)
            + sum(Pmin[i] * γp[t,i] for t in 1:T, i in 1:N)
            - sum(Pmax[i] * γn[t,i] for t in 1:T, i in 1:N)
            + sum(PDmin * x[i] * ωp[t,i] for t in 1:T, i in 1:N)
            - sum(PDmax * x[i] * ωn[t,i] for t in 1:T, i in 1:N)
            + sum(Emin * x[i] * μp[t,i] for t in 1:T, i in 1:N)
            - sum(Emax * x[i] * μn[t,i] for t in 1:T, i in 1:N)
            + sum(E0 * x[i] * τ0[i] for i in 1:N)
            - sum(x[i] * ϕ[t,i] for t in 1:T, i in 1:N)
            + sum(vch4[t,i] for t in 1:T, i in 1:N)
        )
    end
    @constraint(model, sum(x) ≤ b)

    set_optimizer_attribute(model, MOI.Silent(), !verbose)
    set_optimizer_attribute(model, "MIPGap", ϵ_gap)
    optimize!(model)
    return model
end

function solve_trilevel_battery(;
        penalty::Real=0.0,
        verbose::Bool=false,
        max_iter::Int=100,
        ϵ_gap::Float64=5e-3,
    )
    # initialize
    ub = 1e6
    lb = -1e6
    Zub = []
    Zlb = []
    x_hist = []
    y_hist = []
    s_hist = []

    x_ = zeros(N)       # initialize with some solution
    push!(x_hist, x_)

    @printf "%4s \t %4s \t %20s \t %40s \t %10s\n" "iter" "time" "battery location" "interdiction" "gap"
    t_start = now()
    for iter = 1:max_iter
        t_now = (now() - t_start).value / 1000
        (t_now < time_limit) || continue

        # solve upper bound
        if method == "vch"
            model_ub = solve_upper_bound_vch(x_; penalty, verbose, ϵ_gap)
        elseif method == "reg" || method == "lpr"
            model_ub = solve_upper_bound(x_; penalty, verbose, ϵ_gap)
        else
            @error "Invalid method $(method)"
        end

        primal_status(model_ub) != MOI.FEASIBLE_POINT && break
        ub = min(ub, objective_bound(model_ub))
        y_ = value.(model_ub[:y])
        if method == "vch"
            s_ = (
                value.(model_ub[:y]), value.(model_ub[:zp]), value.(model_ub[:zn]), 
                value.(model_ub[:α]), value.(model_ub[:βp]), value.(model_ub[:βn]), 
                value.(model_ub[:γp]), value.(model_ub[:γn]), value.(model_ub[:τ]), value.(model_ub[:τ0]), 
                value.(model_ub[:νp]), value.(model_ub[:νn]), value.(model_ub[:ωp]), value.(model_ub[:ωn]),
                value.(model_ub[:μp]), value.(model_ub[:μn]), value.(model_ub[:ϕ]),
                value.(model_ub[:vch1]), value.(model_ub[:vch2]), value.(model_ub[:vch3]), value.(model_ub[:vch4])
            )
        else
            s_ = (
                value.(model_ub[:y]), value.(model_ub[:zp]), value.(model_ub[:zn]), 
                value.(model_ub[:α]), value.(model_ub[:βp]), value.(model_ub[:βn]), 
                value.(model_ub[:γp]), value.(model_ub[:γn]), value.(model_ub[:τ]), value.(model_ub[:τ0]), 
                value.(model_ub[:νp]), value.(model_ub[:νn]), value.(model_ub[:ωp]), value.(model_ub[:ωn]),
                value.(model_ub[:μp]), value.(model_ub[:μn]), value.(model_ub[:ϕ])
            )
        end
        push!(y_hist, y_)
        push!(s_hist, s_)

        # check termination
        gap = round.(abs.(ub - lb) / abs.(ub) * 100, digits = 2)
        push!(Zlb, lb)
        push!(Zub, ub)

        battery_node = findall(x->x>0,x_)
        interdiction = [L[e] for e in findall(x->x>0, y_)]
        @printf "%4s \t %4.0f \t %20s \t %40s \t %10s\n" iter t_now battery_node interdiction gap

        gap < 0.5 && break

        # solve lower bound
        if method == "vch"
            model_lb = solve_lower_bound_vch(s_hist; verbose, ϵ_gap)
        elseif method == "reg" || method == "lpr"
            model_lb = solve_lower_bound(s_hist; verbose, ϵ_gap)
        else
            @error "Invalid method $(method)"
        end

        primal_status(model_lb) != MOI.FEASIBLE_POINT && break
        lb = max(lb, objective_value(model_lb))
        x_ = value.(model_lb[:x])[:,end]
        if x_ ∈ x_hist
            x_ = zeros(N)
            x_[rand(1:N,b)] .= 1
        end
        push!(x_hist, x_)
    end

    # record the end time
    t_end = (now() - t_start).value / 1000

    return x_hist, y_hist, Zlb, Zub, t_end
end

function save(path, dict)
    data = JSON.json(dict, 4)
    open(path, "w") do f
        write(f, data)
    end
    return nothing
end

# penalty
if method == "relaxation" || method == "vch"
    penalty = 0.0
elseif method == "regularized"
    penalty = round((1-e^2) / (1+e^2), digits =5)
end
x_hist, y_hist, Zlb, Zub, t_solve = solve_trilevel_battery(;penalty)

result = Dict{Any,Any}(
    "battery_location" => findall(x->x>1e-3,x_hist[end]),
    "interdicted_lines" => [L[e] for e in findall(x->x>1e-3, y_hist[end])],
    "cost" => Zub[end],
    "solve_time" => t_solve,
    "x_hist" => x_hist,
    "y_hist" => y_hist,
    "Zlb" => Zlb,
    "Zub" => Zub,
)
result["demand"] = demand
result["budget"] = b
result["contingency"] = k
result["battery_params"] = Dict(
    "PCmax" => PCmax,
    "PCmin" => PCmin,
    "PDmax" => PDmax,
    "PDmin" => PDmin,
    "Emax" => Emax,
    "Emin" => Emin,
    "E0" => E0,
    "e" => e,
)
result["method"] = method
result["network"] = network_name
save(joinpath(output_path, "$(network_name)_$(method).json"), result)
