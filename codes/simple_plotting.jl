using JLD
using Plots

save_path = "E:/UC Davis/Research/high_low_spin/results/"

files = [
    "E:/UC Davis/Research/high_low_spin/results/high_low_spin_L1=10_L2=10_Dz=0.1_Delta=[0.0, 0.5, 1.0, 1.5, 2.0]_sweeps=1000000Date_Sat_06_Nov_2021_07_46_58.jld",
]

file = files[1]
result = load(file, "sim")[1]
params = load(file, "sim")[2]

Quant_name = "E"
Dz_num = 5

#=plot(
    title = "$Quant_name : Dz = $(params[4]), Sweeps = $(params[5])", #ylims=(0,0.03)
)=#

Temp = result[1][Dz_num]
Quantity = []
sample_len = length(result[2][1])
for ind = 1:sample_len
    #push!(Quantity, (params[2]*params[3])*(result[4][Dz_num][ind][1]+result[4][Dz_num][ind][2])/Temp[ind])
    push!(Quantity, result[2][Dz_num][ind])
end
#(params[2]*params[3])*(result[4][Dz_num][ind][1]+result[4][Dz_num][ind][2])

plot!(
    Temp,
    Quantity,
    dpi = 800,
    xlabel = "T",
    ylabel = "$Quant_name",
    legend = :bottomright,
    label = "$Quant_name ($(params[2]) by $(params[3])), Delta = $(params[6][Dz_num])",
)

#savefig("E:/UC Davis/Research/high_low_spin/plots/$(Quant_name)_$(params[1])_antiferro_vary_n_Tmax=5.png")
