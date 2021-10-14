include("confs.jl")
using Plots
using JLD
using Dates
using Statistics

"""
Deciding whether to accept a proposed change or not.

Parameter:
micstate::Microstate : The current state.
numbering::Int64 : A tuple of spin index to flip.

Return:
True of false. True means to flip false means not to.
"""
function accept(micstate::Microstate, numbering::Int64, Sp_new, T::Float64)
    ΔE = Energy_Diff(micstate, numbering, Sp_new)

    p = exp(-ΔE / T) / (1 + exp(-ΔE / T))
    r = rand()
    if r < p
        return true
    else
        return false
    end
end

"""
Container to store value of measurement for each microstate in the Markov Chain.
"""
struct Sample
    E::Array{Float64,1}
    S::Array{Array{Float64,1},1}
    function Sample()
        new(Float64[], Array{Float64,1}[])
    end
end

"""
Function to generate and store samples for a single temperature.

Parameter:
T::Float64 : Temperature
steps::Int64 : steps of the Monte Carlo Simulation
cutoff::Int64 : After the cutoff do we start to take down the samples.
σ::Float64 : Variation of the proposed new value for the Potts value each step in the simulation.

Output:
Samples::Array{Microstate,1} : An Array of samples, each entry is of type Microstate, representing one Microstate.
"""
function sampling(;
    T::Float64,
    sweeps::Int64 = 10^4,
    thermal_cutoff::Int64 = Int(floor(0.4 * sweeps)),
    σ::Float64 = 1.0,
    init_micstate::Microstate,
)
    micstate = init_micstate
    L1 = init_micstate.config.L1
    L2 = init_micstate.config.L2
    samples = Sample()
    for i = 1:sweeps
        for numbering = 1:L1*L2
            cos_θ_new = rand(Uniform(-1, 1))
            θ_new = acos(cos_θ_new)
            ϕ_new = rand(Uniform(0.0, 2 * π))
            sp_new = Spin(θ = θ_new, ϕ = ϕ_new)
            #val = mod1(conf.conf[m,n] + Int(ceil(rand(Normal(0, σ)))), conf.Q)
            if accept(micstate, numbering, sp_new, T)
                update(micstate, numbering, sp_new)
            end
        end

        i < thermal_cutoff + 1 && continue

        push!(samples.E, micstate.E)
        push!(samples.S, copy(micstate.S))
    end
    #println("length of sample is ", length(samples))
    return (samples, micstate.config)
end

@time result = sampling(
    T = 0.01,
    sweeps = 10^6,
    σ = 1.0,
    init_micstate = Microstate(J=1.0, Dz=0.1,config=Lattice(4,4)),
    )
