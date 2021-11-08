include("MC_high_low_spin.jl")

function Eavg_site(E_sample::Array; L1, L2)
    samplelength = length(E_sample)
    Eavg = 0
    for i = 1:samplelength
        Eavg += E_sample[i]
    end
    Eavg = Eavg / samplelength/(L1*L2)
    return Eavg
end

function C_site(E_sample::Array;T, L1, L2)
    samplelength = length(E_sample)
    Eavg = Eavg_site(E_sample, L1=L1, L2=L2)
    E2avg = 0
    for i = 1:samplelength
        E2avg += (E_sample[i]/(L1*L2))^2
    end
    E2avg = E2avg / samplelength
    C = 1 / T * (E2avg - Eavg^2) * (L1*L2)
    return C
end

function S2xyz(S_sample::Array; L1, L2)
    samplelength = length(S_sample)
    S2xyz_avg = [0,0,0]
    for i = 1:samplelength
        S2xyz_avg += S_sample[i].*S_sample[i]
    end
    S2xyz_avg = S2xyz_avg / samplelength / (L1*L2)^2
    return S2xyz_avg
end

function χ_xy(S_sample::Array; L1, L2, Temp)
    samplelength = length(S_sample)
    S2xyz_avg = S2xyz(S_sample; L1=L1, L2=L2)
    S_parallel = 0 #parallel part of S on lattice
    for i = 1:samplelength
        S_parallel += sqrt(S_sample[i][1]^2 + S_sample[i][2]^2)
    end
    S_parallel = S_parallel / samplelength / (L1*L2)^2
    χ_xy = 1/Temp * (S2xyz_avg[1] + S2xyz_avg[2] - S_parallel^2)
    return χ_xy
end

function S2xyz_stag(S_stag_sample::Array; L1, L2)
    samplelength = length(S_stag_sample)
    S2xyz_stag_avg = [0,0,0]
    for i = 1:samplelength
        S2xyz_stag_avg += S_stag_sample[i].*S_stag_sample[i]
    end
    S2xyz_stag_avg = S2xyz_stag_avg / samplelength / (L1*L2)^2
    return S2xyz_stag_avg
end

function B(S_sample::Array;S2xyz, L1, L2)
    samplelength = length(S_sample)
    S4_avg = 0
    for i = 1:samplelength
        S4_avg += sum(S_sample[i].*S_sample[i])^2
    end
    S4_avg = S4_avg / samplelength / (L1*L2)^4
    S2_avg = sum(S2xyz)
    B = 1 - S4_avg / ((1+2/(3-1)) * S2_avg^2)
    return S4_avg
end

"""
Run the simulation for a series of different temps.

Output:

Simulation::Array{Array,1} : It's an array that is of length 5.
    1st entry : An array of ranges of temperatures that is simulated for different J1 values. Each entry correspond to one Δ value.
    2nd entry : An array of Average E values for the simulated ranges of temperatures.
    3rd entry : An array of Average C values for the simulated ranges of temperatures.
    4th entry : An Array of Average [Sx^2,Sy^2,Sz^2] values for the simulated ranges of temperatures(of different Δ values).
    5th entry : An array of Average Binder's ratios for the simulated ranges of temperatures.
    6th entry : An array of Average staggered [Sx^2,Sy^2,Sz^2] for the simulated ranges of temperatures.
    7th entry : An array of Average norm n for the simulated ranges of temperatures.
    8th entry : An array of Average χ on the lattice plane(x,y plane) for the simulated ranges of temperatures.
"""
function driver(;
    L1 = 6,
    L2 = 6,
    Δ = [0.1],
    Dz = 0.1,
    sweeps = 10^6,
    Temp = range(0.01,0.7,length=50),
    norm_update = True
)
    simulation = []
    Temperatures = []
    #L = 16##more variables can be specified here
    Esim = []
    Csim = []
    S2xyzsim = []
    Bsim = []
    S2xyz_stag_sim = []
    n_sim = []
    χ_xy_sim = []
    for (k,Δ) in enumerate(Δ)
        Esamples = []
        Csamples = []
        Bsamples = []
        S2xyzsamples = []
        S2xyz_stag_samples = []
        n_samples = []
        χ_xy_samples = []
        init_config = Lattice(L1,L2)
        for T in Temp
            samples = sampling(
                T = T,
                sweeps = sweeps,
                σ = 1.0,
                init_micstate = Microstate(J=1.0, Δ=Δ, Dz=Dz,config=init_config),
                norm_update = norm_update
                )
            Eavg_site_T = Eavg_site(samples[1].E; L1=L1, L2=L2)
            C_site_T = C_site(samples[1].E;T=T, L1=L1, L2=L2)
            S2xyz_T = S2xyz(samples[1].S; L1=L1, L2=L2)
            S2xyz_stag_T = S2xyz_stag(samples[1].S_stag; L1=L1, L2=L2)
            B_T = B(samples[1].S; S2xyz=S2xyz_T, L1=L1, L2=L2)
            n_T = mean(samples[1].n)
            χ_xy_T = χ_xy(samples[1].S; L1=L1, L2=L2, Temp=T)

            push!(Esamples, Eavg_site_T)
            push!(Csamples, C_site_T)
            push!(S2xyzsamples, S2xyz_T)
            push!(Bsamples, B_T)
            push!(S2xyz_stag_samples, S2xyz_stag_T)
            push!(n_samples, n_T)
            push!(χ_xy_samples, χ_xy_T)
            #output the final configuration as the initial conf for the next temperature.
            init_config = samples[2]
        end
        push!(Temperatures, Temp)
        push!(Esim, Esamples)
        push!(Csim, Csamples)
        push!(S2xyzsim, S2xyzsamples)
        push!(Bsim, Bsamples)
        push!(S2xyz_stag_sim, S2xyz_stag_samples)
        push!(n_sim, n_samples)
        push!(χ_xy_sim, χ_xy_samples)

    end
    push!(simulation, Temperatures)
    push!(simulation, Esim)
    push!(simulation, Csim)
    push!(simulation, S2xyzsim)
    push!(simulation, Bsim)
    push!(simulation, S2xyz_stag_sim)
    push!(simulation, n_sim)
    push!(simulation, χ_xy_sim)
    return simulation
end
L1 = 6;L2 = 6;
Δ = [0.0, 0.5, 1.0, 1.5, 2.0]; sweeps = 10^6; Temp = [range(0.01, 1.5, length=40)...,range(1.5, 0.01, length=40)...]; norm_update = true; Dz = 0.1

content = "high_low_spin_L1=$(L1)_L2=$(L2)_Dz=$(Dz)_Delta=$(Δ)_sweeps=$(sweeps)"
@time sim = driver(L1 = L1, L2 = L2, Dz = Dz, Δ = Δ, sweeps = sweeps, Temp = Temp, norm_update = norm_update)

my_time = Dates.now()
save_path = "/nfs/home/zyt329/Research/high_low_spin/results/"

time_finished = "Date_$(Dates.format(my_time, "e_dd_u_yyyy_HH_MM_SS"))"
save_name = save_path*content*time_finished*".jld"
save(save_name, "sim", [sim,(content, L1, L2, Dz, sweeps, Δ, norm_update)])

println(content*" finished at "*time_finished)
