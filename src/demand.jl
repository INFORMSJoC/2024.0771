using Random
using LinearAlgebra
using DelimitedFiles
using PowerModels

const demand_benchmark = [386836, 375487, 367811, 364067, 361245, 361645, 365598, 374857, 386862, 396082, 398125, 394677, 389635, 384745, 379314, 375475, 376334, 382493, 398515, 411869, 419882, 422048, 417696, 407309]

include("./network.jl")

# user input
network_name = ARGS[1]
noise_level = parse(Float64, ARGS[2])
output_path = "./"
if !isdir(output_path)
    mkdir(output_path)
end

network_path    = "./data/network"

# read in network
network = generate_network_instance(joinpath(network_path, "$(network_name).m"))
N = network.N 
T = 24
Pmax = network.Pmax

# generate and save scenarios
demand_nominal = network.D
if sum(demand_nominal) < sum(Pmax) * 0.8
    scaler = ceil(sum(Pmax) * 0.8 / sum(demand_nominal), digits = 1)
else
    scaler = 1.0
end
demand = zeros(T, N)
for t in 1:T, i in 1:N
    demand[t,i] = demand_nominal[i] * demand_benchmark[t] / demand_benchmark[1] * scaler
end
demand = round.(demand .+ randn(T,N) .* demand .* noise_level, digits = 2)
writedlm(joinpath(output_path, "$(network_name)_demand.txt"), demand)
