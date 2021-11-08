using Random, Distributions

include("lattice.jl")
println()

"""
Calculate energy of a certain configuration. The convention for the cross product is that spin at odd sites cross spin at even site.
"""
function Energy(;J=1, Δ, Dz, config::Lattice)
    spin_conf = config.spin_conf
    nearest_neib_list = config.nearest_neib_list
    N = config.L1 * config.L2
    ####start calculating energy
    E = 0
    for sp in spin_conf
        E += Δ * sp.norm
    end
    for i in config.odd_sites
        sp_i = spin_conf[i]
        for j in nearest_neib_list[i]
            sp_j = spin_conf[j]
            E += J*sp_i.norm*sp_j.norm*((sp_i*sp_j) + Dz*(sp_i ⨱ sp_j))
        end
    end

    return E
end

@debug E = Energy(;J=1, Δ=1, Dz=0.1, config=Lattice(2,2,spin_conf=[Spin(ϕ=0),Spin(ϕ=π/2),Spin(ϕ=π/2),Spin(ϕ=0)]))

function cal_Sx(config::Lattice)
    Sx = 0
    for sp in config.spin_conf
        Sx += sp.norm*sin(sp.θ)*cos(sp.ϕ)
    end
    return Sx
end

function cal_Sy(config::Lattice)
    Sy = 0
    for sp in config.spin_conf
        Sy += sp.norm*sin(sp.θ)*sin(sp.ϕ)
    end
    return Sy
end

function cal_Sz(config::Lattice)
    Sz = 0
    for sp in config.spin_conf
        Sz += sp.norm*cos(sp.θ)
    end
    return Sz
end
@debug Sx=cal_Sx(Lattice(2,2))
@debug Sy=cal_Sy(Lattice(2,2))
@debug Sz=cal_Sz(Lattice(2,2))

function cal_Sx_stag(config::Lattice)
    Sx = 0
    for (ind,sp) in enumerate(config.spin_conf)
        coord = coordinate(ind; L1=config.L1, L2=config.L2)
        Sx += (-1)^(coord[1]+coord[2])*sp.norm*sin(sp.θ)*cos(sp.ϕ)
    end
    return Sx
end

function cal_Sy_stag(config::Lattice)
    Sy = 0
    for (ind,sp) in enumerate(config.spin_conf)
        coord = coordinate(ind; L1=config.L1, L2=config.L2)
        Sy += (-1)^(coord[1]+coord[2])*sp.norm*sin(sp.θ)*sin(sp.ϕ)
    end
    return Sy
end

function cal_Sz_stag(config::Lattice)
    Sz = 0
    for (ind,sp) in enumerate(config.spin_conf)
        coord = coordinate(ind; L1=config.L1, L2=config.L2)
        Sz += (-1)^(coord[1]+coord[2])*sp.norm*cos(sp.θ)
    end
    return Sz
end
@debug println()
@debug println(cal_Sx_stag(Lattice(2,2,spin_conf=[Spin(norm=1.0,ϕ=0),Spin(norm=1.0,ϕ=π),Spin(norm=1.0,ϕ=π),Spin(norm=1.0,ϕ=0)])))
@debug println(cal_Sy_stag(Lattice(2,2)))
@debug println(cal_Sz_stag(Lattice(2,2)))

function cal_n(config::Lattice)
    n = 0
    for (ind,sp) in enumerate(config.spin_conf)
        n += sp.norm
    end
    n = n / (config.L1*config.L2)
    return n
end

@debug println(cal_n(Lattice(2,2)))

"""
Represent a Microstate.

fields:  T, L, steps, dimension, spin
"""
mutable struct Microstate
    J::Float64
    Δ::Float64
    Dz::Float64
    config::Lattice
    E::Float64
    S::Array{Float64,1}
    S_stag::Array{Float64, 1}
    n::Float64

    function Microstate(;
        J = 1.0,
        Δ::Float64 = 1.0,
        Dz::Float64 = 0.1,
        config::Lattice = Lattice(2,2),
    )
        E = Energy(;J=J, Δ=Δ, Dz=Dz, config=config)
        S = [cal_Sx(config),cal_Sy(config),cal_Sz(config)]
        S_stag = [cal_Sx_stag(config),cal_Sy_stag(config),cal_Sz_stag(config)]
        n = cal_n(config)
        new(J, Δ, Dz, config, E, S, S_stag, n)
    end
end
@debug micstate = Microstate()


"""
Calculate the energy difference after rotating one spin.

Parameter:
micstate::Microstate : Input the
numbering::Int64 : the numbering of the rotated spin
"""
function Energy_Diff(micstate::Microstate, numbering::Int64, sp_new::Spin)
    J = micstate.J
    Δ = micstate.Δ
    Dz = micstate.Dz
    L1 = micstate.config.L1
    L2 = micstate.config.L2
    nearest_neib_list = micstate.config.nearest_neib_list
    spin_conf = micstate.config.spin_conf
    sp_old = spin_conf[numbering]
    @assert(
        0 < numbering <= L1*L2,
        "The numbering of rotated spin should be within range"
    )
    coord = coordinate(numbering,L1=L1, L2=L2)
    ######Start calculating energy difference before and after
    ΔE = Δ * (sp_new.norm - sp_old.norm)
    for neib_num in nearest_neib_list[numbering]
        sp_neib = spin_conf[neib_num]
        if isodd(coord[1]+coord[2]) #selected site is odd site
            ΔE += J*sp_new.norm*sp_neib.norm*((sp_new*sp_neib) + Dz*(sp_new ⨱ sp_neib))
            ΔE -= J*sp_old.norm*sp_neib.norm*((sp_old*sp_neib) + Dz*(sp_old ⨱ sp_neib))
        else #selected site is even site
            ΔE += J*sp_neib.norm*sp_new.norm*((sp_neib*sp_new) + Dz*(sp_neib ⨱ sp_new))
            ΔE -= J*sp_neib.norm*sp_old.norm*((sp_neib*sp_old) + Dz*(sp_neib ⨱ sp_old))
        end
    end
    return ΔE
end
@debug ΔE = Energy_Diff(micstate, 1, Spin(ϕ=π/2))

"""
Calculate the total magnetization S difference of a Microstate after rotating one spin.

"""
function S_Diff(micstate::Microstate, numbering::Int64, sp_new::Spin)
    L1 = micstate.config.L1
    L2 = micstate.config.L2
    spin_conf = micstate.config.spin_conf
    sp_old = spin_conf[numbering]
    @assert(
        0 < numbering <= L1*L2,
        "The numbering of rotated spin should be within range"
    )
    coord = coordinate(numbering, L1=L1, L2=L2)
    ######Start calculating S difference before and after
    ΔS = [sp_new.norm*sin(sp_new.θ)*cos(sp_new.ϕ) - sp_old.norm*sin(sp_old.θ)*cos(sp_old.ϕ),
    sp_new.norm*sin(sp_new.θ)*sin(sp_new.ϕ) - sp_old.norm*sin(sp_old.θ)*sin(sp_old.ϕ),
    sp_new.norm*cos(sp_new.θ) - sp_old.norm*cos(sp_old.θ)]
    return ΔS
end
@debug ΔS=S_Diff(Microstate(), 1, Spin(ϕ=π/2))

"""
Calculate the average norm n difference of a Microstate after rotating one spin.
"""
function n_Diff(micstate::Microstate, numbering::Int64, sp_new::Spin)
    L1 = micstate.config.L1
    L2 = micstate.config.L2
    spin_conf = micstate.config.spin_conf
    sp_old = spin_conf[numbering]
    @assert(
        0 < numbering <= L1*L2,
        "The numbering of rotated spin should be within range"
    )
    ######Start calculating n difference before and after
    Δn = (sp_new.norm - sp_old.norm) / (L1*L2)
    return Δn
end
@debug println(n_Diff(Microstate(), 1, Spin(norm=0)))

"""
Calculate the staggered total magnetization S difference of a Microstate after rotating one spin.
"""
function S_Diff_stag(micstate::Microstate, numbering::Int64, sp_new::Spin)
    L1 = micstate.config.L1
    L2 = micstate.config.L2
    spin_conf = micstate.config.spin_conf
    sp_old = spin_conf[numbering]
    @assert(
        0 < numbering <= L1*L2,
        "The numbering of rotated spin should be within range"
    )
    coord = coordinate(numbering, L1=L1, L2=L2)
    ######Start calculating S difference before and after
    ΔS = (-1)^(coord[1]+coord[2]).*[sp_new.norm*sin(sp_new.θ)*cos(sp_new.ϕ) - sp_old.norm*sin(sp_old.θ)*cos(sp_old.ϕ),
    sp_new.norm*sin(sp_new.θ)*sin(sp_new.ϕ) - sp_old.norm*sin(sp_old.θ)*sin(sp_old.ϕ),
    sp_new.norm*cos(sp_new.θ) - sp_old.norm*cos(sp_old.θ)]
    return ΔS
end
@debug println(S_Diff_stag(Microstate(), 1, Spin(ϕ=π/2)))

"""
Updating the current Microstate to a new Microstate. Also change the energy and S.

Parameter:
micstate::Microstate : current configuration.
numbering::Int64 : numbering of the changed spin.
sp_new::Spin : updated spin value for the selected spin.
"""
function update(micstate::Microstate, numbering::Int64, sp_new::Spin)
    micstate.E += Energy_Diff(micstate, numbering, sp_new)
    micstate.S += S_Diff(micstate, numbering, sp_new)
    micstate.S_stag += S_Diff_stag(micstate, numbering, sp_new)
    micstate.n += n_Diff(micstate, numbering, sp_new)
    micstate.config.spin_conf[numbering] = sp_new

    nothing
end
micstate = Microstate()
update(micstate, 1, Spin(ϕ=π/2))
