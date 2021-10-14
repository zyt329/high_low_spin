import Base.*

"""
Represent a spin on a site. Angles are in radian([0,2π))

fields:
norm::Float64 : Total spin of the spin
θ::Float64 : θ angle the spin makes with z axis, takes value [0,π)
ϕ::Float64 : ϕ angle the spin makes with x axis, takes value [0,2π)
"""
mutable struct Spin
    norm::Float64
    θ::Float64
    ϕ::Float64
    function Spin(;norm=1.0, θ=π/2, ϕ=0.0)
        norm = Float64(norm)
        #make sure θ is within 0 to π
        θ = Float64(θ)
        θ1 = mod(θ, 2π)
        θ = min(θ1, 2π-θ1)
        #make sure ϕ is within 0 to 2π
        ϕ = Float64(mod(ϕ, 2*π))
        new(norm, θ, ϕ)
    end
end

"""
Define spin dot product. Returns the dot product of two spins.
"""
function *(sp1::Spin, sp2::Spin)
    norm_part = sp1.norm * sp2.norm
    if norm_part == 0
        return 0.0
    else
        dir_part = sin(sp1.θ)*sin(sp2.θ)*cos(sp1.ϕ - sp2.ϕ) + cos(sp1.θ)*cos(sp2.θ)
        return norm_part * dir_part
    end
end

"""
    Computing Z component of dot product of 2 spins
"""
function ⨱(sp1::Spin, sp2::Spin)
    norm_part = sp1.norm * sp2.norm
    if norm_part == 0
        return 0.0
    else
        dir_part = - sin(sp1.θ)*sin(sp2.θ)*sin(sp1.ϕ - sp2.ϕ)
        return norm_part * dir_part
    end
end


"""
    Rotate the spin by Δθ and Δϕ.
    Subtlety:
    1. when θ is rotated to be larger than π or smaller than 0, it should also add to ϕ an additional π.
"""
function rotate(sp::Spin; Δθ=0, Δϕ=0)
    Δθ=Float64(Δθ); Δϕ=Float64(Δϕ)
    # update θ
    sp.θ += Δθ
    θ1 = mod(sp.θ, 2π)
    if !(0<=θ1<=π) # we need to add π to ϕ
        sp.θ = 2π-θ1
        sp.ϕ += π # Put π to proper range later(in the update-ϕ step)
    end
    # update ϕ
    sp.ϕ += Δϕ
    sp.ϕ = mod(sp.ϕ, 2π)
    println("θ=$(sp.θ), ϕ=$(sp.ϕ)")
    return sp
end

function test()
    sp=Spin(θ=0,ϕ=0)
    println(sp)
    for i = 1:7
        sp=rotate(sp,Δθ=π/3)
    end
end
#test()
#rotate(Spin(θ=5π/6,ϕ=0), Δθ=2π/6)
