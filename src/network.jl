struct network_structure
    name::String
    N::Int
    E::Int
    L::Vector{Tuple}  # edges

    Pmax::Vector{Float64}   # generator upper limit
    Pmin::Vector{Float64}   # generator lower limit
    Fmax::Vector{Float64}   # transmission line limit
    B::Vector{Float64}      
    D::Vector{Float64}      
    δin::Dict{Int,Vector}
    δout::Dict{Int,Vector}
end

function generate_network_instance(network_file)   
    network_data = PowerModels.parse_file(network_file)

    N = length(network_data["bus"])
    E = length(network_data["branch"])
    T = 24

    V = sort!(parse.(Int, collect(keys(network_data["bus"]))))
    VtoN = Dict(V[i] => i for i in 1:N)

    L = [(VtoN[branch["f_bus"]], VtoN[branch["t_bus"]]) for (b_id, branch) in network_data["branch"]]

    Pmax = zeros(N)
    Pmin = zeros(N)
    Fmax = zeros(E)
    B = zeros(E)
    D = zeros(N)
    for (id, gdat) in network_data["gen"]
        v = gdat["gen_bus"]
        i = VtoN[v]
        Pmax[i] = Pmax[i] + gdat["pmax"]
        Pmin[i] = Pmin[i] + (gdat["pmin"] < 1e-3 ? gdat["pmax"] / 3 : gdat["pmin"])
    end
    for (id, edat) in network_data["branch"]
        e = parse(Int, id)
        Fmax[e] = edat["rate_a"]
        B[e] = - edat["br_x"] / (edat["br_x"]^2 + edat["br_r"]^2)
    end
    for (id, ldat) in network_data["load"]
        v = ldat["load_bus"]
        i = VtoN[v]
        D[i] = D[i] + max(ldat["pd"], 0)
    end
    δin  = Dict(i => [] for i in 1:N)
    δout = Dict(i => [] for i in 1:N)
    for (e, l) in enumerate(L)
        i = l[1]
        j = l[2]
        push!(δout[i], e)
        push!(δin[j], e)
    end

    return network_structure(network_name,N,E,L,Pmax,Pmin,Fmax,B,D,δin,δout)
end
