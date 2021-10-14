include("confs.jl")

function spherical2xyz(θ,ϕ;norm=1)
    x = round(norm*sin(θ)*cos(ϕ),digits=3)
    y = round(norm*sin(θ)*sin(ϕ),digits=3)
    z = round(norm*cos(θ),digits=3)
    println("[$x, $y, $z]")
end

function print_vec(sp_conf::Lattice)
    for (i,sp) in enumerate(sp_conf.spin_conf)
        print("$i : \n")
        spherical2xyz(sp.θ,sp.ϕ;norm=0.5)
        print("\n")
    end
end

print_vec(result[2])

#spherical2xyz(1.5095325598940483, 2.961253128145517)

function print_lattice(L=4)
    for i in 1:L
        for j in 1:L
            print("($i, $j, 0)")
        end
    end
end
print_lattice()
