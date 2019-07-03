include("./helper.jl")

localModulePath = "/home/fireofearth/Research/mitchell-ian/2019s-verification/verification/module"
if(!(localModulePath in LOAD_PATH))
    push!(LOAD_PATH, localModulePath)
end

using Test, Random, IntervalArithmetic
using ModalIntervalArithmetic
using AffineArithmetic

@testset "modal interval getters" begin
    X = ModalInterval(1,2) # proper
    Y = ModalInterval(2,1) # improper
    Z = ModalInterval(2,2) # real_number
    dX = dual(X)
    dY = dual(Y)
    dZ = dual(Z)
    @test inf(X) == 1 && sup(X) == 2
    @test inf(dX) == 2 && sup(dX) == 1
    @test inf(Y) == 2 && sup(Y) == 1
    @test inf(dY) == 1 && sup(dY) == 2
    @test inf(Z) == sup(Z) == 2
    @test inf(dZ) == sup(dZ) == 2
    @test dX == Y && dY == X
end

@testset "modal interval comparison" begin
    X = ModalInterval(2,3)
    Y = ModalInterval(1,4)
    dX = ModalInterval(3,2)
    dY = ModalInterval(4,1)
    # ([2,3]',∃) = [2,3] ⊆ [1,4] = ([1,4]',∃)
    # => [2,3]' ⊆ [1,4]'
    @test  X ⊆ Y  &&  Y ⊇ X
    # ([2,3]',∀) = [3,2] ⊇ [4,1] = ([1,4]',∀)
    # => [2,3]' ⊆ [1,4]'
    @test dX ⊇ dY && dY ⊆ dX 
end

@testset "modal interval additive group" begin
    A = ModalInterval(1,2)
    B = ModalInterval(3,4)
    @test A + B == ModalInterval(1+3,2+4)
    @test -A == ModalInterval(-2,-1)
    @test A - B == ModalInterval(1-4,2-3)
    @test A * B == ModalInterval(1 * 3, 2 * 4)
end

@testset "modal interval multiplicative group" begin
    #A = ModalInterval()
end

function getGroupSolCoeffs(
    f::Function, ind1::Vector, ind2::Vector ,dev1::Vector, dev2::Vector
)
    mm   = max(last(ind1), last(ind2))
    rg   = Array(1:mm)
    addSol = Array{Union{Nothing,AbstractFloat}}(nothing, mm)
    for ii in rg
        idx1 = findfirst(x -> x == ii, ind1)
        idx2 = findfirst(x -> x == ii, ind2)
        if(idx1 != nothing && idx2 != nothing)
            addSol[ii] = f(dev1[idx1], dev2[idx2])
        elseif(idx1 != nothing)
            addSol[ii] = f(dev1[idx1], 0)
        elseif(idx2 != nothing)
            addSol[ii] = f(0, dev2[idx2])
        end
    end
    return sort(
        filter(x -> x[2] != nothing, tuple.(rg,addSol)), 
        by = y -> y[1]
    )
end

#=
function getLeftModuleSolCoeffs(
    f::Function, ind::Vector, dev::Vector
)
    rg = last(ind)
    #for ii in rg
end
=#

@testset "affine arithmetic" begin
    center1 = rand(Float64) + rand(0:99)
    n1 = rand(6:10)
    dev1 = rand(Float64,n1) .+ rand(-9:9, n1)
    ind1 = sort(shuffle(Array(1:20))[1:n1])
    a1 = AAF(center1,dev1,ind1)

    center2 = rand(Float64) + rand(0:99)
    n2 = rand(6:10)
    dev2 = rand(Float64,n2) .+ rand(-9:9, n2)
    ind2 = sort(shuffle(Array(1:20))[1:n2])
    a2 = AAF(center2,dev2,ind2)

    @testset "affine-affine addition" begin
        aAddn = a1 + a2
        @test aAddn[0] ≈ center1 + center2
        addnSol = getGroupSolCoeffs(+,ind1,ind2,dev1,dev2)
        @test length(aAddn) == length(aAddn.indexes) == length(addnSol)
        addnComp = [(aAddn.indexes[ii], dev) for (ii, dev) in enumerate(aAddn.deviations)]
        for ((compIdx, compDev), (solIdx, solDev)) in tuple.(addnComp,addnSol)
            @test compIdx == solIdx
            @test compDev ≈ solDev
        end
    end

    @testset "affine-affine substraction" begin
        aSubt = a1 - a2
        @test aSubt[0] ≈ center1 - center2
        subtSol = getGroupSolCoeffs(-,ind1,ind2,dev1,dev2)
        @test length(aSubt) == length(aSubt.indexes) == length(subtSol)
        subtComp = [(aSubt.indexes[ii], dev) for (ii, dev) in enumerate(aSubt.deviations)]
        for ((compIdx, compDev), (solIdx, solDev)) in tuple.(subtComp,subtSol)
            @test compIdx == solIdx
            @test compDev ≈ solDev
        end
    end

    @testset "affine negative" begin
        aNeg = -a2
        @test aNeg[0] == -center2
        negSol = tuple.(ind2, -dev2)
        @test length(aNeg) == length(aNeg.indexes) == length(negSol)
        negComp = [(aNeg.indexes[ii], dev) for (ii, dev) in enumerate(aNeg.deviations)]
        for ((compIdx, compDev), (solIdx, solDev)) in tuple.(negComp,negSol)
            @test compIdx == solIdx
            @test compDev ≈ solDev
        end
    end

    @testset "affine add constant" begin
    end
end

