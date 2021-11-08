using JLD
using DelimitedFiles

include("confs.jl")

function spherical2xyz(θ,ϕ;norm=1)
    x = round(norm*sin(θ)*cos(ϕ),digits=3)
    y = round(norm*sin(θ)*sin(ϕ),digits=3)
    z = round(norm*cos(θ),digits=3)
    return [x,y,z]
end

function print_vec(sp_conf::Lattice)
    for (i,sp) in enumerate(sp_conf.spin_conf)
        print("$i : \n")
        spherical2xyz(sp.θ,sp.ϕ;norm=0.5)
        print("\n")
    end
end

print_vec(result[2])


function printing(;result = result, name = "test")

    open(
        "E:/UC Davis/Research/high_low_spin/results_txt/"*name*".txt",
        "w",
    ) do io
        for (i,sp) in enumerate(result[2].spin_conf)
            start = coordinate(i; L1=result[2].L1, L2=result[2].L2)
            vec = spherical2xyz(sp.θ,sp.ϕ;norm=sp.norm*1)
            writedlm(io, [start... 0 vec...])
        end
    end

end

printing(name="test")
