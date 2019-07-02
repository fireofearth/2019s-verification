using Test, Random

using IntervalArithmetic

include("helper.jl")
include("ModalInterval.jl")
include("Affine.jl")

@testset "modal interval getters" begin
    X = ModalInterval(1,2) # proper
    Y = ModalInterval(2,1) # improper
    dX = dual(X)
    @test inf(X) == 1
    @test sup(X) == 2
    @test inf(dX) == 2
    @test sup(dX) == 1
    @test inf(Y) == 2
    @test sup(Y) == 1
    @test dX == Y
end

@testset "modal interval arithmetic" begin
    A = ModalInterval(1,2)
    B = ModalInterval(3,4)
    @test A + B == ModalInterval(1+3,2+4)
    @test -A == ModalInterval(-2,-1)
    @test A - B == ModalInterval(1-4,2-3)
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
        filter(x -> x[2] != nothing, pair.(rg,addSol)), 
        by = y -> y[1]
    )
end

function getLeftModuleSolCoeffs(
    f::Function, ind::Vector, dev::Vector
)
    rg = last(ind)
    #for ii in rg
end

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
        for ((compIdx, compDev), (solIdx, solDev)) in pair.(addnComp,addnSol)
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
        for ((compIdx, compDev), (solIdx, solDev)) in pair.(subtComp,subtSol)
            @test compIdx == solIdx
            @test compDev ≈ solDev
        end
    end

    @testset "affine negative" begin
        aNeg = -a2
        @test aNeg[0] == -center2
        negSol = pair.(ind2, -dev2)
        @test length(aNeg) == length(aNeg.indexes) == length(negSol)
        negComp = [(aNeg.indexes[ii], dev) for (ii, dev) in enumerate(aNeg.deviations)]
        for ((compIdx, compDev), (solIdx, solDev)) in pair.(negComp,negSol)
            @test compIdx == solIdx
            @test compDev ≈ solDev
        end
    end

    @testset "affine add constant" begin
    end
end

