module AffineArithmetic 
#=
 # Affine arithmetic module
 #
 # This is a implementation of the An Affine Arithmetic C++
 # Library ( http://aaflib.sourceforge.net ) in Julia
 #
 # Specification:
 # - we use Float64 to represent double, and Vector{Float64}
 # to represent arrays in heap memory.
 # - Julia operations rely on the use of the `last` index to
 # keep track of indeterminate coefficients
 #
 # TODO: rename AAF => Affine
 # TODO: documentation?
 # TODO: aaflib arithmetic functions +,-,*,/
 # TODO: aaflib trig. functions sin, cos
 # TODO: aaflib power functions pow, ^
 # TODO: figure out what changes Goubault/Putot made for
 # aaflib
 # TODO: Goubault/Putot functions mult_eps, hull
=#

localModulePath = "/home/fireofearth/Research/mitchell-ian/2019s-verification/verification/module"
if(!(localModulePath in LOAD_PATH))
    push!(LOAD_PATH, localModulePath)
end

import IntervalArithmetic: Interval

# since we will likely use AffineArithmetic along with IntervalArithmetic, we want to avoid namespace conflicts so we will import inf, sup here
import IntervalArithmetic: inf, sup

import Base:
    zero, one, iszero, isone, convert,
    getindex, length, repr, size, firstindex, lastindex,
    <, <=, >, >=, ==, +, -, *, /, inv

export
    AAFCoeff, AAFInd, AAF,
    convert,
    getindex, length, repr, firstindex, lastindex,
    <, <=, >, >=,
    Interval,
    rad, getMax, getMin, getAbsMax, getAbsMin, inv,
    ==, +, -, *, /

# TODO: testing only
export getLastAAFIndex, resetLastAAFIndex, ApproximationType

 #=
 # Module-wide constants
 #
 # TODO: are we using the right constants?
=#
@enum ApproximationType MINRANGE CHEBYSHEV SECANT
MINRAD = 1E-10
TOL = 1E-15
EPSILON = 1E-20

 #=
 # Type declarations
=#
AAFCoeff = Float64
AAFInd = Int64
approximationType = CHEBYSHEV

 #=
 # last keeps record of last coefficient index of affine forms accoss all AAF instances.
 #
 # Specification:
 # - force setLastAFFIndex to call whenever a new AAF instance is created.
 #
 # TODO: turn this into a decorator/macro and force calls
 # TODO: which functions is `last` for exactly?
=#
let lastAAFIndex::Int = 0

    global function resetLastAAFIndex()
        lastAAFIndex = 0
    end

     #=
     # Similar to getDefault() in aaflib
    =#
    global function getLastAAFIndex()
        return lastAAFIndex
    end

     #=
     # Assign new coefficient index i for noise symbol μᵢ.
     # To be used when assigning computation noise.
     # Similar to inclast() in aaflib
    =#
    global function addAAFIndex()
        lastAAFIndex += 1
        return [lastAAFIndex]
    end

    global function addAAFIndex(indexes::Vector{AAFInd})
        lastAAFIndex += 1
        return vcat(indexes, lastAAFIndex)
    end

     #=
     # TODO: finish and make sure this is correct
    =#
    global function setLastAAFIndex(indexes::Vector{AAFInd})
        m = last(indexes)
        if(m > lastAAFIndex)
            lastAAFIndex = m
        end
        return indexes
    end
end

 #=
 # AAF represents an affine form
 #
 # Specification:
 # - AAF is a collection with deviations as its elements.
 #
 # Invariants:
 # - AAF indexes are always in sorted order from lowest to highest
 # - elts in AAF indexes are unique
 #
 # TODO: it's not clear whether `length` and `size` are redundant
 # TODO: refactor to get rid of length and size, can query vector length
 # TODO: we don't need to write most getters, since AAF is immutable
 # TODO: make this inherit AbstractArrays
 # TODO: enable iterator, indexing
=#
struct AAF
    cvalue::AAFCoeff  # central value 
    #length::Int # length of indexes 
    #size::Int   # array size of indexes and deviations
    deviations::Vector{AAFCoeff}
    indexes::Vector{AAFInd}

     #=
     # Creates an AAF without deviations.
    =#
    AAF(v0::AAFCoeff = 0.0) = new(v0, Vector{AAFCoeff}(), Vector{AAFInd}())

     #=
     # Creates an AAF with deviations
     # TODO: check viability of removint 't' from constructor
     # TODO: manage assertions in dev/production versioning
    =#
    function AAF(v0::AAFCoeff, dev::Vector{AAFCoeff}, ind::Vector{AAFInd})
        @assert length(ind) == length(dev)
        for ii in 1:(length(ind) - 1)
            @assert ind[ii] <= ind[ii + 1]
        end
        new(v0, dev, setLastAAFIndex(ind))
    end
    
     #=
     # Constructor from intervals
    =#
    AAF(iv::Interval) = new(mid(iv), [radius(iv)], addAAFIndex())

     #=
     # Constructor that assigns new center to AAF
    =#
    AAF(a::AAF, cst::AAFCoeff) = new(cst, a.deviations, a.indexes)

     #=
     # Constructor that assigns new center, and diff to AAF. We assume indexes unchanged
    =#
    function AAF(a::AAF, cst::AAFCoeff, dev::Vector{AAFCoeff})
        @assert length(dev) == length(a.deviations)
        new(cst, dev, a.indexes)
    end
end

 #=
 # Copy constructor are implicitly supported
 # similar to `function AAF(p::AAF)`
 #
 # Julia has no explicit type assigning, so we skip all assignment overloading.
=#

# current approximation type: <CHEBYSHEV> (default), <MINRANGE> or <SECANT> 
# approximationType::tApproximationType

# highest deviation symbol in use
# last::Int

function getindex(a::AAF, ind::Int)::AAFCoeff
    if(ind < 0 || ind > length(a.deviations))
        return 0.0
    elseif(ind == 0)
        return a.cvalue
    else
        return a.deviations[ind]
    end
end

function length(a::AAF)
    return length(a.deviations)
end

function repr(a::AAF)
    s = "$(a[0])"
    if(length(a) > 0)
        for i in 1:length(a)
            s *= " + $(a[i])ϵ$(a.indexes[i])"
        end
    end
    return s
end

convert(::Type{AAF}, x::Number) = AAF(x)
one(::Type{AAF}) = convert(AAF, 1.)
one(x::AAF) = convert(AAF, 1.)
zero(::Type{AAF}) = convert(AAF, 0.)
zero(x::AAF) = convert(AAF, 0.)

 #=
 # Get the total deviation of an AAF (i.e. the sum of all deviations (their abs value))
=#
rad(a::AAF)::AAFCoeff = sum(abs.(a.deviations))

Interval(a::AAF) = Interval(a[0] - rad(a), a[0] + rad(a))

 #=
 # Goubault+Putot methods
 # Specification: convert_int, reduce_aaf, rad
 # TODO: implement
=#

  #=
  # Get maximum / minimum of an AAF
 =#
getCenter(a::AAF) = a[0]
getMax(a::AAF)::AAFCoeff = a[0] + rad(a)
getMin(a::AAF)::AAFCoeff = a[0] - rad(a)
getAbsMax()::AAFCoeff = max(abs(a[0] - rad(a)), abs(a[0] + rad(a)))
getAbsMin()::AAFCoef  = min(abs(a[0] - rad(a)), abs(a[0] + rad(a)))

firstindex(a::AAF) = length(a) > 0 ? a.indexes[1] : 0
lastindex(a::AAF)  = length(a) > 0 ? last(a.indexes) : 0 

#zero(x::DataType) = 
#one
#iszero
#isone

 #=
 # Unknown methods
 # Specification: compact, sumup
 # TODO: implement
=#

  #=
  # Conditionals
 =#
<(a::AAF,  p::AAF) = (a[0] + rad(a)) <  (p[0] - rad(p))
<=(a::AAF, p::AAF) = (a[0] + rad(a)) <= (p[0] - rad(p))
>(a::AAF,  p::AAF) = (a[0] - rad(a)) >  (p[0] + rad(p))
>=(a::AAF, p::AAF) = (a[0] - rad(a)) >= (p[0] + rad(p))

 #=
 # Equality
 # Specification: affine equality compares cvalues, and deviations TOL=1E-15
 # TODO: is it better to use ≈ instead?
 # TODO: refactor length out
=#
function ==(a::AAF, p::AAF)
    if(length(a) != length(p))
        return false
    end

    if(abs(a[0]) < 1 && abs(p[0]) < 1)
        if(abs(a[0] - p[0]) > TOL)
            return false
        end
    else
        if(abs((a[0] - p[0]) / (a[0] + p[0])) > TOL)
            return false
        end
    end

    for i in 1:length(a)
        if(a.indexes[i] != p.indexes[i])
            return false
        end
        if(abs(a[i]) < 1 && abs(p[i]) < 1)
            if(abs(a[i] - p[i]) > TOL)
                return false
            end
        else
            if(abs((a[i] - p[i]) / 
                   (a[i] + p[i])) > TOL)
                return false
            end
        end
    end
    
    return true
end

+(a::AAF, cst::AAFCoeff)::AAF = AAF(a, a[0] + cst)
-(a::AAF, cst::AAFCoeff)::AAF = AAF(a, a[0] - cst)
+(cst::AAFCoeff, a::AAF)::AAF = AAF(a, cst + a[0])
-(cst::AAFCoeff, a::AAF)::AAF = AAF(a, cst - a[0])

*(a::AAF, cst::AAFCoeff)::AAF = AAF(a, a[0] * cst, cst * a.deviations)
*(cst::AAFCoeff, a::AAF)::AAF = AAF(a, cst * a[0], cst * a.deviations)

function /(a::AAF, cst::AAFCoeff)::AAF
    @assert cst != 0.0
    AAF(a, a.cvalue * (1.0 / cst), (1.0 / cst) * deviations)
end

 #=
 # a + p, where a, p are AAF
=#
function +(a::AAF, p::AAF)::AAF
    if(length(p) == 0)
        return a + p[0]
    elseif(length(a) == 0)
        return a[0] + p
    end
    indt = [ii for ii in union(a.indexes, p.indexes)]
    sort!(indt)
    devt = fill(0.0, length(indt)) #  Vector(undef,length(indt))
    pcomp = tuple.(indexin(indt, a.indexes), indexin(indt, p.indexes))
    for (ii, (ia,ip)) in enumerate(pcomp)
        @assert ia != nothing || ip != nothing
        devt[ii] = (ia == nothing ? p[ip] : 
         (ip == nothing ? a[ia] : a[ia] + p[ip]))
    end
    return AAF(a[0] + p[0], devt, indt)
end

 #=
 # a - p where a, p are AAF
=#
function -(a::AAF, p::AAF)::AAF
    if(length(p) == 0)
        return a - p[0]
    elseif(length(a) == 0)
        return a[0] - p
    end
    indt = [ii for ii in union(a.indexes, p.indexes)]
    sort!(indt)
    devt = fill(0.0, length(indt)) #  Vector(undef,length(indt))
    pcomp = tuple.(indexin(indt, a.indexes), indexin(indt, p.indexes))
    for (ii, (ia,ip)) in enumerate(pcomp)
        @assert ia != nothing || ip != nothing
        devt[ii] = (ia == nothing ? -p[ip] : 
         (ip == nothing ? a[ia] : a[ia] - p[ip]))
    end
    return AAF(a[0] - p[0], devt, indt)
end

function -(a::AAF)::AAF
    return AAF(a, -a[0], -1 * a.deviations)
end

 #=
 # Approximates a * p where a, p are AAF
 # There are three affine products based on the appoximation of the coefficient for μₖ
 #   xy = x₀ŷ₀ + ∑ᴺᵢ(xᵢy₀+yᵢx₀)ϵᵢ + ½∑[over 1⩽i,j⩽n] |xᵢyⱼ+yᵢxⱼ|μₖ
 #   xy = x₀ŷ₀ + ∑ᴺᵢ(xᵢy₀+yᵢx₀)ϵᵢ + (∑ᴺᵢ|xᵢ|)(∑ᴺᵢ|yᵢ|)μₖ
 #   xy = x₀ŷ₀ + ½∑ᴺᵢxᵢyᵢ + ∑ᴺᵢ(xᵢy₀+yᵢx₀)ϵᵢ + [(∑ᴺᵢ|xᵢ|)(∑ᴺᵢ|yᵢ|) - ½∑ᴺᵢ|xᵢyᵢ|]μₖ
 # The last approximation is obtainable by observing that for products of like terms 
 # we have squares of noise symbols: xᵢyᵢϵᵢ² 
 # Since ϵᵢ² ∈ [0,1] the term has a center at ½xᵢyᵢ with magnitude of deviation ½|xᵢyᵢ|
 #
 # Specification:
 #   xy = x₀ŷ₀ + ½∑ᴺᵢxᵢyᵢ + ∑ᴺᵢ(xᵢy₀+yᵢx₀)ϵᵢ + [(∑ᴺᵢ|xᵢ|)(∑ᴺᵢ|yᵢ|) - ½∑ᴺᵢ|xᵢyᵢ|]μₖ
=#
function *(a::AAF, p::AAF)::AAF
    if(length(p) == 0)
        return a * p[0]
    elseif(length(a) == 0)
        return a[0] * p
    end
    # create new index with length = length(a.indexes) + length(p.indexes) + 1
    adjDeviation2 = 0.0
    adjCenter2    = 0.0
    indt  = [ii for ii in union(a.indexes, p.indexes)]
    sort!(indt)
    indt  = addAAFIndex(indt)
    lindt = length(indt)
    devt  = fill(0.0, lindt)
    pcomp = tuple.(indexin(indt, a.indexes), indexin(indt, p.indexes))
    for (ii, (ia,ip)) in enumerate(pcomp)
        if(ia == nothing && ip == nothing) # happens at the end of indt
            # do nothing
        elseif(ia == nothing)
            devt[ii] = a[0] * p[ip]
        elseif(ip == nothing)
            devt[ii] = a[ia] * p[0]
        else
            devt[ii] = a[ia] * p[0] + a[0] * p[ip]
            adjDeviation2 = abs(a[ia] * p[ip])
            adjCenter2    = a[ia] * p[ip]
        end
    end
    devt[lindt] = rad(a)*rad(p) - 0.5*adjDeviation2
    return AAF(a[0]*p[0] + 0.5*adjCenter2, devt, indt)
end

 #=
 # 1/a where a is AAF
=#
function inv(p::AAF)::AAF
    if(length(p) == 0)
        return AAF(1.0 / p[0])
    end

    a = p[0] - r;
    b = p[0] + r;
    if(a*b < EPSILON)
        throw(DomainError(p, "trying to invert zero"))
    end

    inva = 1. / b
    invb = 1. / b

    if(approximationType = CHEBYSHEV)
        if(r > MINRAD)
            alpha = (invb - inva) / (b - a)
        else
            alpha = inva
        end
        u = log(alpha)
        delta = 0.5*(inva + (u - a - 1.0)*alpha)
        dzeta = inva - a*alpha - delta
    elseif(approximationType = MINRANGE)
        error("incomplete")
    else # if(approximationType = SECANT)
        error("incomplete")
    end

    indt = addAAFIndex(p.indexes)
    devt = alpha * p.deviations
    devt = vcat(devt, delta)
    return AAF(alpha*p[0] + dzeta, devt, indt)
end

/(cst::AAFCoeff, a::AAF) = cst * inv(a)
/(a::AAF, cst::AAFCoeff) = AAF(a, a[0] / cst, a.deviations / cst)
/(a::AAF, p::AAF)::AAF   = a * inv(p)

function ^(p::AAF, n::Int)
    if(length(p) == 0)
        return AAF(p[0]^n)
    end

    if(n == 0)
        return one(AAF)
    elseif(n == 1)
        return p
    elseif(n == -1)
        return inv(p)
    end
    
    r = rad(p)
    a = p[0] - r
    b = p[0] + r
    pa = a^n
    pb = b^n
    if(a*b < EPSILON && n < 0)
        throw(DomainError(p, "trying to invert zero"))
    end

    if(approximationType = CHEBYSHEV)
        if(r > MINRAD)
            alpha = (pb - pa) / (b - a)
        else
            alpha = n * pa / (a + EPSILON)
        end

        xₐ = - abs(alpha / n)^(1 / (n - 1))
        xᵦ = -xₐ

        if(xₐ > a)
            pxₐ = xₐ^n
        else
            xₐ  = a
            pxₐ = pa
        end

        if(xᵦ < b)
            pxᵦ = xᵦ^n
        else
            xᵦ = b
            pxᵦ = pb
        end

        yₐ= pxₐ - alpha*xₐ
        yᵦ= pxᵦ - alpha*xᵦ
        delta = 0.5*(yₐ - yᵦ)
        dzeta = 0.5*(yₐ + yᵦ)

    elseif(approximationType = MINRANGE)
        error("incomplete")
    else # if(approximationType = SECANT)
        error("incomplete")
    end

    indt = addAAFIndex(p.indexes)
    devt = alpha * p.deviations
    devt = vcat(devt, delta)
    return AAF(alpha*p[0] + dzeta, devt, indt)
end

#function ^(a::AAF) const;

end # module AffineArithmetic
