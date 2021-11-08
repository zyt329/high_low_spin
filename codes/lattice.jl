include("spin.jl")

function coordinate(n::Int64; L1::Int64=2, L2::Int64=2)
    @assert((n ≤ L1 * L2) && (1 ≤ n),"The numbering (bit position) of a site shouldn't exceed the total number of sites $(L1 * L2), and should be bigger than 0.")
    i::Int64 = Int(ceil(n/L1))
    j::Int64 = mod1(n,L1)  #site i is at i-th row, j-th column
    return (i,j)
end

function numbering(coordinate::Tuple{Int64,Int64};L1::Int64=2,L2::Int64=2)
    @assert((coordinate[1] ≤ L2) && (coordinate[2] ≤ L1),"The cooridnate should be within the range of the lattice size $L1 by $L2")
    n = (coordinate[1]-1)*L1 + coordinate[2]
    return n
end

function neib(n::Int64;L1::Int64=2, L2::Int64=2)
    coord = coordinate(n,L1=L1, L2=L2)
    neibs = Tuple{Int64,Int64}[]
    push!(neibs, (mod1(coord[1]+1,L2), coord[2]))
    push!(neibs, (mod1(coord[1]-1,L2), coord[2]))
    push!(neibs, (coord[1], mod1(coord[2]+1,L1)))
    push!(neibs, (coord[1], mod1(coord[2]-1,L1)))
    if iseven(coord[1]+coord[2])
        push!(neibs, (mod1(coord[1]+1,L2), mod1(coord[2]-1,L1)))
        push!(neibs, (mod1(coord[1]-1,L2), mod1(coord[2]+1,L1)))
    else
        push!(neibs, (mod1(coord[1]+1,L2), mod1(coord[2]+1,L1)))
        push!(neibs, (mod1(coord[1]-1,L2), mod1(coord[2]-1,L1)))
    end
    #=convert coordinations to positions in bits=#
    neibs_numbering = Set{Int64}()
    for neib in neibs
        push!(neibs_numbering, numbering(neib, L1=L1, L2=L2))
    end
    return neibs_numbering
end

function nearest_neib(n::Int64;L1::Int64=2, L2::Int64=2)
    coord = coordinate(n,L1=L1, L2=L2)
    neibs = Tuple{Int64,Int64}[]
    push!(neibs, (mod1(coord[1]+1,L2), coord[2]))
    push!(neibs, (mod1(coord[1]-1,L2), coord[2]))
    push!(neibs, (coord[1], mod1(coord[2]+1,L1)))
    push!(neibs, (coord[1], mod1(coord[2]-1,L1)))
    #=convert coordinations to positions in bits=#
    neibs_numbering = Set{Int64}()
    for neib in neibs
        push!(neibs_numbering, numbering(neib, L1=L1, L2=L2))
    end
    return neibs_numbering
end

function second_neib(n::Int64;L1::Int64=2, L2::Int64=2)
    coord = coordinate(n,L1=L1, L2=L2)
    neibs = Tuple{Int64,Int64}[]
    if iseven(coord[1]+coord[2])
        push!(neibs, (mod1(coord[1]+1,L2), mod1(coord[2]-1,L1)))
        push!(neibs, (mod1(coord[1]-1,L2), mod1(coord[2]+1,L1)))
    else
        push!(neibs, (mod1(coord[1]+1,L2), mod1(coord[2]+1,L1)))
        push!(neibs, (mod1(coord[1]-1,L2), mod1(coord[2]-1,L1)))
    end
    #=convert coordinations to positions in bits=#
    neibs_numbering = Set{Int64}()
    for neib in neibs
        push!(neibs_numbering, numbering(neib, L1=L1, L2=L2))
    end
    return neibs_numbering
end

function neib_list_gen(;L1::Int64=2, L2::Int64=2)
    neib_list = Set{Int64}[]
    for n in 1:L1*L2
        push!(neib_list, neib(n, L1=L1, L2=L2))
    end
    return neib_list
end

function nearest_neib_list_gen(;L1::Int64=2, L2::Int64=2)
    neib_list = Set{Int64}[]
    for n in 1:L1*L2
        push!(neib_list, nearest_neib(n, L1=L1, L2=L2))
    end
    return neib_list
end

function second_neib_list_gen(;L1::Int64=2, L2::Int64=2)
    neib_list = Set{Int64}[]
    for n in 1:L1*L2
        push!(neib_list, second_neib(n, L1=L1, L2=L2))
    end
    return neib_list
end

function odd_even_site_list_gen(;L1::Int64=2, L2::Int64=2)
    odd_even_site_list = Dict()
    odd_even_site_list[:odd] = Int64[]
    odd_even_site_list[:even] = Int64[]
    for i in 1:L1*L2
        if isodd(sum(coordinate(i; L1=L1, L2=L2)))
            push!(odd_even_site_list[:odd], i)
        else
            push!(odd_even_site_list[:even], i)
        end
    end
    return odd_even_site_list
end

function init_lattice_gen(;L1,L2)
    lattice = []
    ϕ_init = rand(Uniform(0,2π))
    for i in 1:L1*L2
        if isodd(sum(coordinate(i;L1=L1,L2=L2)))
            push!(lattice, Spin(ϕ=ϕ_init))
        else
            push!(lattice, Spin(ϕ=ϕ_init+π))
        end
    end
    return lattice
end

mutable struct Lattice
    L1::Int64
    L2::Int64
    odd_sites::Array{Int64,1}
    even_sites::Array{Int64,1}
    nearest_neib_list::Array{Set{Int64},1}
    second_neib_list::Array{Set{Int64},1}
    all_neib_list::Array{Set{Int64},1}
    spin_conf::Array{Spin, 1}
    function Lattice(
        L1,L2;
        nearest_neib_list=nearest_neib_list_gen(L1=L1, L2=L2),
        second_neib_list=second_neib_list_gen(L1=L1, L2=L2),
        all_neib_list=neib_list_gen(L1=L1, L2=L2),
        spin_conf = init_lattice_gen(L1=L1,L2=L2)
        )
        odd_sites = odd_even_site_list_gen(L1=L1, L2=L2)[:odd]
        even_sites = odd_even_site_list_gen(L1=L1, L2=L2)[:even]

        new(L1,L2,
        odd_sites,
        even_sites,
        nearest_neib_list,
        second_neib_list,
        all_neib_list,
        spin_conf)
    end
end

L1=2;L2=2
init_config = Lattice(2,2)
