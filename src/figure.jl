using JuMP
using Gurobi
const GRB_ENV = Gurobi.Env()

using Plots
using Measures

T = 2
N = 2
E = 1
L = [(1,2)]

Pmin = [2, 2]
Pmax = [4, 4]
Fmax = [1]
PCmin = 0
PDmin = 0
PCmax = 1
PDmax = 1
Emin = 0
Emax = 6
E0 = 0

D = [5 1; 8 4]
x = [0, 1]

range_gen = 0.001:0.001:1.0
for denom in [1.9, 2.0, 2.1]
    eta = 1/sqrt(denom)
    objval = Any[]
    for λ in range_gen
        λc = λ
        λd = λ
        model = JuMP.Model(() -> Gurobi.Optimizer(GRB_ENV))
        @variable(model, u[1:T, 1:N], Bin)
        @variable(model, f[1:T, 1:E])
        @variable(model, pg[1:T, 1:N] ≥ 0)
        @variable(model, pc[1:T, 1:N] ≥ 0)
        @variable(model, pd[1:T, 1:N] ≥ 0)
        @variable(model, pls[1:T, 1:N] ≥ 0)
        @variable(model, pex[1:T, 1:N] ≥ 0)
        @variable(model, ps[0:T, 1:N] ≥ 0)
    
        flow_constraints = @constraint(model, [t in 1:T, i in 1:N], sum(f[t,e] for (e,l) in enumerate(L) if l[1] == i) - sum(f[t,e] for (e,l) in enumerate(L) if l[2] == i) == pg[t,i] - D[t,i] - pc[t,i] + pd[t,i] + pls[t,i] - pex[t,i])
        @constraint(model, [t in 1:T, e in 1:E], -Fmax[e] ≤ f[t,e])
        @constraint(model, [t in 1:T, e in 1:E], Fmax[e] ≥ f[t,e])
        @constraint(model, [t in 1:T, i in 1:N], ps[t,i] == ps[t-1,i] + eta * pc[t,i] - 1/eta * pd[t,i])
        @constraint(model, [t in 1:T, i in 1:N], ps[t,i] ≥ Emin * x[i])
        @constraint(model, [t in 1:T, i in 1:N], ps[t,i] ≤ Emax * x[i])
        @constraint(model, [i in 1:N], ps[0,i] == E0 * x[i])
        @constraint(model, [t in 1:T, i in 1:N], pc[t,i] ≥ PCmin * u[t,i])
        @constraint(model, [t in 1:T, i in 1:N], pc[t,i] ≤ PCmax * u[t,i])
        @constraint(model, [t in 1:T, i in 1:N], pd[t,i] ≥ PDmin * (x[i]-u[t,i]))
        @constraint(model, [t in 1:T, i in 1:N], pd[t,i] ≤ PDmax * (x[i]-u[t,i]))
        @constraint(model, [t in 1:T, i in 1:N], u[t,i] ≤ x[i])
        @constraint(model, [t in 1:T, i in 1:N], pg[t,i] ≥ Pmin[i])
        @constraint(model, [t in 1:T, i in 1:N], pg[t,i] ≤ Pmax[i])
    
        @objective(model, Min, sum(pls) + sum(pex) + λc * sum(pc) + λd * sum(pd))
        set_optimizer_attribute(model, MOI.Silent(), true)
        optimize!(model)
        if termination_status(model) == MOI.OPTIMAL
            push!(objval, sum(value.(pls)) + sum(value.(pex)))
        end
    end

    t_stop = findfirst(x->x>3.6, objval)
    plot(range_gen[1:t_stop], vcat(objval[1:t_stop-1], objval[t_stop-1]),
        xlabel = "λ",
        ylabel = "c(p)",
        label = "λ ∈ [0.0, 1.0]",
        xlims = (0.0,1.0),
        ylims = (3.2,4.3),
        xtickfontsize=15,
        ytickfontsize=15,
        xguidefontsize=15,
        yguidefontsize=15,
        linewidth = 8,
        legend = :bottomright,
        legendfontsize=15,
        margin=5mm,
    )
    plot!(range_gen[t_stop:end], objval[t_stop:end],
        label = false,
        linecolor = palette(:default)[1],
        linewidth = 8,
    )
    eta = 1/sqrt(denom)
    lambda = (1-eta^2)/(1+eta^2)
    if lambda > range_gen[t_stop-1]
        plot!([lambda], [objval[end]], seriestype=:scatter,
            label = "λ = (1-e²)/(1+e²)",
            color = palette(:default)[2],
            legendfontsize=15,
            markersize=12,
        )
    else
        plot!([lambda], [objval[1]], seriestype=:scatter,
            label = "λ = (1-e²)/(1+e²)",
            color = palette(:default)[2],
            legendfontsize=15,
            markersize=12,
        )
    end
    savefig("./figure2_$(denom).png")
end

diff_est = [T * sum(x) * max(PCmax * (1-eta^2)/(1+eta^2), PDmax * (1-eta^2)/(1+eta^2)) for eta in range_gen]
diff_tru = []
for eta in range_gen
    λc = (1-eta^2)/(1+eta^2)
    λd = λc
    model = JuMP.Model(() -> Gurobi.Optimizer(GRB_ENV))
    @variable(model, u[1:T, 1:N], Bin)
    @variable(model, f[1:T, 1:E])
    @variable(model, pg[1:T, 1:N] ≥ 0)
    @variable(model, pc[1:T, 1:N] ≥ 0)
    @variable(model, pd[1:T, 1:N] ≥ 0)
    @variable(model, pls[1:T, 1:N] ≥ 0)
    @variable(model, pex[1:T, 1:N] ≥ 0)
    @variable(model, ps[0:T, 1:N] ≥ 0)

    flow_constraints = @constraint(model, [t in 1:T, i in 1:N], sum(f[t,e] for (e,l) in enumerate(L) if l[1] == i) - sum(f[t,e] for (e,l) in enumerate(L) if l[2] == i) == pg[t,i] - D[t,i] - pc[t,i] + pd[t,i] + pls[t,i] - pex[t,i])
    @constraint(model, [t in 1:T, e in 1:E], -Fmax[e] ≤ f[t,e])
    @constraint(model, [t in 1:T, e in 1:E], Fmax[e] ≥ f[t,e])
    @constraint(model, [t in 1:T, i in 1:N], ps[t,i] == ps[t-1,i] + eta * pc[t,i] - 1/eta * pd[t,i])
    @constraint(model, [t in 1:T, i in 1:N], ps[t,i] ≥ Emin * x[i])
    @constraint(model, [t in 1:T, i in 1:N], ps[t,i] ≤ Emax * x[i])
    @constraint(model, [i in 1:N], ps[0,i] == E0 * x[i])
    @constraint(model, [t in 1:T, i in 1:N], pc[t,i] ≥ PCmin * u[t,i])
    @constraint(model, [t in 1:T, i in 1:N], pc[t,i] ≤ PCmax * u[t,i])
    @constraint(model, [t in 1:T, i in 1:N], pd[t,i] ≥ PDmin * (x[i]-u[t,i]))
    @constraint(model, [t in 1:T, i in 1:N], pd[t,i] ≤ PDmax * (x[i]-u[t,i]))
    @constraint(model, [t in 1:T, i in 1:N], u[t,i] ≤ x[i])
    @constraint(model, [t in 1:T, i in 1:N], pg[t,i] ≥ Pmin[i])
    @constraint(model, [t in 1:T, i in 1:N], pg[t,i] ≤ Pmax[i])

    @objective(model, Min, sum(pls) + sum(pex) + λc * sum(pc) + λd * sum(pd))
    set_optimizer_attribute(model, MOI.Silent(), true)
    optimize!(model)
    if termination_status(model) == MOI.OPTIMAL
        obj_reg = value.(sum(pls) + sum(pex))
    end
    @objective(model, Min, sum(pls) + sum(pex))
    optimize!(model)
    if termination_status(model) == MOI.OPTIMAL
        obj_ori = objective_value(model)
    end
    push!(diff_tru, obj_reg - obj_ori)
end
plot(range_gen, diff_tru,
    xlabel = "η",
    ylabel = "c(p̂) - c(p̄)",
    label = "True difference",
    xlims = (0.0,1.0),
    ylims = (-0.05,2.05),
    xtickfontsize=15,
    ytickfontsize=15,
    xguidefontsize=15,
    yguidefontsize=15,
    linewidth = 8,
    legend = :topright,
    legendfontsize=15,
    margin=5mm,
    linecolor = palette(:default)[1],
)
plot!(range_gen, diff_est,
    label = "Worst bound",
    linecolor = palette(:default)[2],
    linewidth = 8,
    linestyle = :dash,
)
savefig("./figure2_diff.png")
