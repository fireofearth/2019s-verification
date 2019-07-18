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
 # TODO: finish documentation
 # TODO: expand arithmetic functions for all subtypes of Real
 # TODO: aaflib trig. functions sin, cos; requires me to know
 # what implementation changes Goubault/Putot made for sin, cos
 # TODO: aaflib pow is supported?
 # TODO: figure out what changes Goubault/Putot made for aaflib
 # TODO: Goubault/Putot functions mult_eps, hull
 # TODO: complete + test support for ForwardDiff
=#

import IntervalArithmetic: Interval

# since we will likely use AffineArithmetic along with IntervalArithmetic, we want to avoid namespace conflicts so we will import inf, sup here
import IntervalArithmetic: inf, sup

import Base:
    zero, one, iszero, isone, convert, isapprox,
    getindex, length, repr, size, firstindex, lastindex,
    <, <=, >, >=, ==, +, -, *, /, inv, ^

export
    AffineCoeff, AffineInd, AffineInt, Affine,
    zero, one, convert, isapprox,
    getindex, length, repr, firstindex, lastindex,
    <, <=, >, >=,
    Interval, inf, sup,
    rad, getMax, getMin, getAbsMax, getAbsMin, inv,
    ==, +, -, *, /, ^

# TODO: testing only
export getLastAffineIndex, resetLastAffineIndex, ApproximationType

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
AffineCoeff = Float64 # type for coefficients and constants
AffineInt   = Int64 # type for integers
AffineInd   = Int64 # type for indexes
approximationType = CHEBYSHEV # default rounding mode

disp(msg) = print("$(msg)\n")
debug() = print("DEBUG\n")

 #=
 # last keeps record of last coefficient index of affine forms accoss all Affine instances.
 #
 # Specification:
 # - force setLastAffineIndex to call whenever a new Affine instance is created.
 #
 # TODO: turn this into a decorator/macro and force calls
 # TODO: which functions is `last` for exactly?
=#
let lastAffineIndex::Int = 0

    global function resetLastAffineIndex()
        lastAffineIndex = 0
    end

     #=
     # Similar to getDefault() in aaflib
    =#
    global function getLastAffineIndex()
        return lastAffineIndex
    end

     #=
     # Assign new coefficient index i for noise symbol μᵢ.
     # To be used when assigning computation noise.
     # Similar to inclast() in aaflib
    =#
    global function addAffineIndex()
        lastAffineIndex += 1
        return [lastAffineIndex]
    end

    global function addAffineIndex(indexes::Vector{AffineInd})
        lastAffineIndex += 1
        return vcat(indexes, lastAffineIndex)
    end

     #=
     # TODO: test
    =#
    global function setLastAffineIndex(indexes::Vector{AffineInd})
        m = last(indexes)
        if(m > lastAffineIndex)
            lastAffineIndex = m
        end
        return indexes
    end
end

 #=
 # Affine represents an affine form
 #
 # Specification:
 # - Affine is a collection with deviations as its elements.
 #
 # Invariants:
 # - Affine indexes are always in sorted order from lowest to highest
 # - elts in Affine indexes are unique
 #
 # TODO: it's not clear whether `length` and `size` are redundant
 # TODO: refactor to get rid of length and size, can query vector length
 # TODO: we don't need to write most getters, since Affine is immutable
 # TODO: make this inherit AbstractArrays
 # TODO: enable iterator, indexing
=#
struct Affine <: Real
    cvalue::AffineCoeff  # central value 
    #length::Int # length of indexes 
    #size::Int   # array size of indexes and deviations
    deviations::Vector{AffineCoeff}
    indexes::Vector{AffineInd}

     #=
     # Creates an Affine without deviations.
    =#
    Affine(v0::AffineCoeff = 0.0) = new(v0, Vector{AffineCoeff}(), Vector{AffineInd}())
    Affine(v0::AffineInt) = new(AffineCoeff(v0), Vector{AffineCoeff}(), Vector{AffineInd}())

     #=
     # Creates an Affine with deviations
     # TODO: check viability of removing 't' from constructor
     # TODO: manage assertions in dev/production versioning
    =#
    function Affine(v0::AffineCoeff, dev::Vector{AffineCoeff}, ind::Vector{AffineInd})
        @assert length(ind) == length(dev)
        for ii in 1:(length(ind) - 1)
            @assert ind[ii] <= ind[ii + 1]
        end
        new(v0, dev, setLastAffineIndex(ind))
    end
    
     #=
     # Constructor from intervals
    =#
    Affine(iv::Interval) = new(mid(iv), [radius(iv)], addAffineIndex())

     #=
     # Constructor that assigns new center to Affine
    =#
    Affine(a::Affine, cst::AffineCoeff) = new(cst, a.deviations, a.indexes)

     #=
     # Constructor that assigns new center, and diff to Affine. We assume indexes unchanged
    =#
    function Affine(a::Affine, cst::AffineCoeff, dev::Vector{AffineCoeff})
        @assert length(dev) == length(a.deviations)
        new(cst, dev, a.indexes)
    end
end

 #=
 # Copy constructor are implicitly supported
 # similar to `function Affine(p::Affine)`
 #
 # Julia has no explicit type assigning, so we skip all assignment overloading.
=#

function getindex(a::Affine, ind::Int)::AffineCoeff
    if(ind < 0 || ind > length(a.deviations))
        return 0.0
    elseif(ind == 0)
        return a.cvalue
    else
        return a.deviations[ind]
    end
end

function length(a::Affine)
    return length(a.deviations)
end

function repr(a::Affine)
    s = "$(a[0])"
    if(length(a) > 0)
        for i in 1:length(a)
            s *= " + $(a[i])ϵ$(a.indexes[i])"
        end
    end
    return s
end

convert(::Type{Affine}, x::Number) = Affine(x)
one(::Type{Affine}) = convert(Affine, 1.)
one(x::Affine) = convert(Affine, 1.)
zero(::Type{Affine}) = convert(Affine, 0.)
zero(x::Affine) = convert(Affine, 0.)

 #=
 # Get the total deviation of an Affine (i.e. the sum of all deviations (their abs value))
=#
rad(a::Affine)::AffineCoeff = sum(abs.(a.deviations))

Interval(a::Affine) = Interval(a[0] - rad(a), a[0] + rad(a))

 #=
 # Goubault+Putot methods
 # Specification: convert_int, reduce_aaf, rad
 # TODO: implement
=#

  #=
  # Get maximum / minimum of an Affine
 =#
getCenter(a::Affine) = a[0]
getMax(a::Affine)::AffineCoeff = a[0] + rad(a)
getMin(a::Affine)::AffineCoeff = a[0] - rad(a)
getAbsMax()::AffineCoeff = max(abs(a[0] - rad(a)), abs(a[0] + rad(a)))
getAbsMin()::AffineCoef  = min(abs(a[0] - rad(a)), abs(a[0] + rad(a)))

firstindex(a::Affine) = length(a) > 0 ? a.indexes[1] : 0
lastindex(a::Affine)  = length(a) > 0 ? last(a.indexes) : 0 

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
<(a::Affine,  p::Affine) = (a[0] + rad(a)) <  (p[0] - rad(p))
<=(a::Affine, p::Affine) = (a[0] + rad(a)) <= (p[0] - rad(p))
>(a::Affine,  p::Affine) = (a[0] - rad(a)) >  (p[0] + rad(p))
>=(a::Affine, p::Affine) = (a[0] - rad(a)) >= (p[0] + rad(p))

 #=
 # Equality
 # Specification: affine equality compares cvalues, and deviations up to
 # some tolerance which defaults to TOL
=#
function equalityInterval(a::Affine, p::Affine; tol::Float64=TOL)
    if(length(a) != length(p))
        return false
    end

    if(abs(a[0]) < 1 && abs(p[0]) < 1)
        if(abs(a[0] - p[0]) > tol)
            return false
        end
    else
        if(abs((a[0] - p[0]) / (a[0] + p[0])) > tol)
            return false
        end
    end

    for i in 1:length(a)
        if(a.indexes[i] != p.indexes[i])
            return false
        end
        if(abs(a[i]) < 1 && abs(p[i]) < 1)
            if(abs(a[i] - p[i]) > tol)
                return false
            end
        else
            if(abs(a[i] - p[i]) / 
                   (abs(a[i]) + abs(p[i])) > tol)
                return false
            end
        end
    end
    
    return true
end

==(a::Affine, p::Affine) = equalityInterval(a, p)

function isapprox(a::Affine, p::Affine; tol::Float64=TOL)
    return equalityInterval(a, p; tol=tol)
end

+(a::Affine, cst::AffineCoeff)::Affine = Affine(a, a[0] + cst)
-(a::Affine, cst::AffineCoeff)::Affine = Affine(a, a[0] - cst)
+(cst::AffineCoeff, a::Affine)::Affine = Affine(a, cst + a[0])
-(cst::AffineCoeff, a::Affine)::Affine = Affine(a, cst - a[0])

+(a::Affine, cst::AffineInt)::Affine = Affine(a, AffineCoeff(a[0] + cst))
-(a::Affine, cst::AffineInt)::Affine = Affine(a, AffineCoeff(a[0] - cst))
+(cst::AffineInt, a::Affine)::Affine = Affine(a, AffineCoeff(cst + a[0]))
-(cst::AffineInt, a::Affine)::Affine = Affine(a, AffineCoeff(cst - a[0]))

*(a::Affine, cst::AffineCoeff)::Affine = Affine(a, a[0] * cst, cst * a.deviations)
*(cst::AffineCoeff, a::Affine)::Affine = Affine(a, cst * a[0], cst * a.deviations)

*(a::Affine, cst::AffineInt)::Affine = Affine(a, AffineCoeff(a[0] * cst), 
                                              convert(Vector{AffineCoeff}, cst * a.deviations))
*(cst::AffineInt, a::Affine)::Affine = Affine(a, AffineCoeff(cst * a[0]), 
                                              convert(Vector{AffineCoeff}, cst * a.deviations))

function /(a::Affine, cst::Union{AffineCoeff, AffineInt})::Affine
    if(cst == zero(cst))
        throw(DomainError(a, "trying to divide by zero"))
    end
    Affine(a, a.cvalue * (1.0 / cst), (1.0 / cst) * a.deviations)
end

 #=
 # a + p, where a, p are Affine
=#
function +(a::Affine, p::Affine)::Affine
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
    return Affine(a[0] + p[0], devt, indt)
end

 #=
 # a - p where a, p are Affine
=#
function -(a::Affine, p::Affine)::Affine
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
    return Affine(a[0] - p[0], devt, indt)
end

function -(a::Affine)::Affine
    return Affine(a, -a[0], -1 * a.deviations)
end

 #=
 # Approximates a * p where a, p are Affine
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
function *(a::Affine, p::Affine)::Affine
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
    indt  = addAffineIndex(indt)
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
    return Affine(a[0]*p[0] + 0.5*adjCenter2, devt, indt)
end

 #=
 # Obtain 1/a where a is Affine
 # TODO: explain how this works?
=#
function inv(p::Affine)::Affine
    if(length(p) == 0)
        return Affine(1.0 / p[0])
    end

    r = rad(p)
    a = p[0] - r;
    b = p[0] + r;
    if(a*b < EPSILON)
        throw(DomainError(p, "trying to invert zero"))
    end

    inva = 1. / a
    invb = 1. / b

    if(approximationType == CHEBYSHEV)
        alpha = -inva * invb
        u = sqrt(a * b)

        if(a > 0)
            delta = 0.5*(inva + invb - 2.0/u)
            dzeta = inva + invb - delta
        else
            delta = -0.5*(inva + invb + 2.0/u)
            dzeta = inva + invb + delta
        end
        
    elseif(approximationType == MINRANGE)
        error("incomplete")
    else # if(approximationType == SECANT)
        error("incomplete")
    end

    indt = addAffineIndex(p.indexes)
    devt = alpha * p.deviations
    devt = vcat(devt, delta)
    return Affine(alpha*p[0] + dzeta, devt, indt)
end

/(cst::AffineCoeff, a::Affine) = cst * inv(a)
/(a::Affine, p::Affine)::Affine   = a * inv(p)

 #=
 # Obtain a^n where a is Affine and n is an integer
 # TODO: explain how this works?
=#
function ^(p::Affine, n::Int)
    if(length(p) == 0)
        return Affine(p[0]^n)
    end

    if(n == 0)
        return one(Affine)
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

    if(approximationType == CHEBYSHEV)
        if(r > MINRAD)
            alpha = (pb - pa) / (b - a)
        else
            alpha = n * pa / (a + EPSILON)
        end

        xₐ = -abs(alpha / n)^(1.0 / (n - 1.0))
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
            xᵦ  = b
            pxᵦ = pb
        end

        yₐ= pxₐ - alpha*xₐ
        yᵦ= pxᵦ - alpha*xᵦ
        delta = 0.5*(yₐ - yᵦ)
        dzeta = 0.5*(yₐ + yᵦ)

    elseif(approximationType == MINRANGE)
        error("incomplete")
    else # if(approximationType == SECANT)
        error("incomplete")
    end

    indt = addAffineIndex(p.indexes)
    devt = alpha * p.deviations
    devt = vcat(devt, delta)
    return Affine(alpha*p[0] + dzeta, devt, indt)
end

# TODO
function sin(p::Affine)::Affine
    if(length(p) == 0)
        return Affine(sin(p[0]))
    end

end

#function ^(a::Affine) const;

end # module AffineArithmetic
