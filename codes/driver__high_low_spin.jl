include("MC_newM.jl")
include("Temp_range_gen.jl")
include("error_analysis.jl")

"""
Run the simulation for a series of different temps.

Output:

Simulation::Array{Array,1} : It's an array that is of length 5.
    1st entry : An array of Average E values for the simulated ranges of temperatures.
    2nd entry : An array of Average C values for the simulated ranges of temperatures.
    3rd entry : An array of Average Binder's ratios for the simulated ranges of temperatures.
    4th entry : An array of ranges of temperatures that is simulated for different J1 values. Each entry correspond to one J1 value.
    5th entry : An array of simulated J1 values.
    6th entry : An Array of Average M2 values for the simulated ranges of temperatures(of different J1 values).
    7th entry : An Array of Average Distances along synthetic dimension between the two sites whos indexes are set by "Syn_Dist_inds". Each array corresponds to a simulated range of temperature(of the corresponding J1 valueg).

"""
function driver(;
    L = 8,
    Q = 16,
    J1 = [range(0.0, 1.0, length = 5); 0.4; 0.45],
    sweeps = 10^6,
    Temp = range(0.3,0.7,length=50),
    Syn_Dist_inds=((1,1),(1,2))
)
    #Temp = [0.4, 0.6]
    #L = 16##more variables can be specified here

    Esim = []
    Csim = []
    Bsim = []
    simulation = []
    M2sim = []
    Syn_Distsim = []
    #=Temps = [
        [0.598,0.603],
        [0.573,0.588],
        [0.552,0.557],
        [0.5899,0.5904],
        [0.55,0.60],
        [0.556,0.561],
        [0.555,0.560],
    ]=#
    Temperatures = []
    for k = 1:length(J1)
        #=Temp = range(Temps[k][1], Temps[k][2], length = 5)=#
        #Temp = Temp_range_gen(J = J1[k])
        Esamples = []
        Csamples = []
        Bsamples = []
        M2samples = []
        Syn_Distsamples = []
        init_conf::Array{Int64,2} = ones(Int64, L, L)
        for T in Temp
            samples = sampling(
                T = T,
                J = [1.0, J1[k]],
                L = L,
                Q = Q,
                sweeps = sweeps,
                init_conf = init_conf,
                Syn_Dist_inds = Syn_Dist_inds
            )
            samplelength = length(samples[1].E)
            Eavg = 0
            E2avg = 0
            Syn_Distavg = 0
            for i = 1:samplelength
                Eavg += samples[1].E[i]
                E2avg += (samples[1].E[i])^2
                Syn_Distavg += samples[1].Syn_Dist[i]
            end
            Eavg = Eavg / samplelength
            E2avg = E2avg / samplelength
            Syn_Distavg = Syn_Distavg / (samplelength*Q)
            C = 1 / T^2 * (E2avg - Eavg^2) / L^2
            Eavg = Eavg / L^2
            push!(Esamples, Eavg)
            push!(Csamples, C)
            push!(Syn_Distsamples, Syn_Distavg)
            # Calculating B
            M2 = 0
            M4 = 0
            for M in samples[1].M
                Δ = cal_M2(M, Q)
                M2 += Δ
                M4 += Δ^2
            end
            M2 = M2 / samplelength
            push!(M2samples, M2)
            M4 = M4 / samplelength
            B = 1 - M4 / ((1+2/(Q-1)) * M2^2)
            push!(Bsamples, B)
            #output the final configuration as the initial conf for the next temperature.
            init_conf = samples[2]
        end
        push!(Esim, Esamples)
        push!(Csim, Csamples)
        push!(Bsim, Bsamples)
        push!(Temperatures, Temp)
        push!(M2sim, M2samples)
        push!(Syn_Distsim, Syn_Distsamples)
    end
    push!(simulation, Esim)
    push!(simulation, Csim)
    push!(simulation, Bsim)
    push!(simulation, Temperatures)
    push!(simulation, J1)
    push!(simulation, M2sim)
    push!(simulation, Syn_Distsim)
    #plot(Temp, Esamples)
    #plot!(Temp, Csamples)

    return simulation

end
L = 16; Q = 36
J1 = [0.8];sweeps = 4*10^6;Temp = range(0.01, 1, length=50)

content = "newM_nonperiodic_M2added_avgdist_J1$(J1[1])_$(J1[end])_Q$(Q)_sweep$(sweeps)"
@time sim = driver(L = L, Q = Q, J1 = J1, sweeps = sweeps, Temp = Temp)

my_time = Dates.now()
save_path = "/nfs/home/zyt329/Research/Synthetic_dim_code/paper_result/"
#"E:/UC Davis/Research/Synthetic Dimensions/Synthetic_dim_code/Bcrossing_result_newM_nonperiodic/"

time_finished = "Date_$(Dates.format(my_time, "e_dd_u_yyyy_HH_MM_SS"))"
save_name = save_path*content*"_L_$(L)__Q_$(Q)__sweeps_$(sweeps)_"*time_finished*".jld"
save(save_name, "sim", [sim,(content, L, Q, sweeps)])
#=save(
    "/nfs/home/zyt329/Research/Synthetic_dim_code/Bcrossing_result/Metropolis_Wolff_fineTemp4_L_$(L)__Q_$(Q)__sweeps_10e7_Date_$(Dates.format(my_time, "e_dd_u_yyyy_HH_MM_SS")).jld",
    "sim",
    sim,
)
println("L=$L, Q=$Q, Wolff_fineTemp4_10e7 finished at $(Dates.format(my_time, "e_dd_u_yyyy_HH_MM_SS"))")=#
#=save(
    "E:/UC Davis/Research/Synthetic Dimensions/Synthetic_dim_code/test/"*content*"_L_$(L)__Q_$(Q)__sweeps_$(sweeps)_Date_$(Dates.format(my_time, "e_dd_u_yyyy_HH_MM_SS")).jld",
    "sim",
    sim,
)=#
println("L=$L, Q=$Q,"*content*" finished at "*time_finished)
