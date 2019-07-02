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
 # TODO: aaflib arithmetic functions +,-,*,/
 # TODO: aaflib trig. functions sin, cos
 # TODO: aaflib power functions pow, ^
 # TODO: figure out what changes Goubault/Putot made for
 # aaflib
 # TODO: Goubault/Putot functions mult_eps, hull
=#

using IntervalArithmetic, StaticArrays

include("ModalInterval.jl")

import IntervalArithmetic

import Base

import Base:
    ==

 #=
 # Type declarations
=#
AAFCoeff = Float64
AAFInd = Int64

@enum tApproximationType MINRANGE CHEBYSHEV SECANT

 #=
 # last keeps record of last coefficient index of affine forms accoss all AAF instances.
 #
 # Specification:
 # - force setLastAFFIndex to call whenever a new AAF instance is created.
 #
 # TODO: turn this into a decorator/macro and force calls
 # TODO: which functions is `last` for exactly?
=#
let last::Int = 0

    global function resetLastAAFIndex()
        last = 0
    end

     #=
     # Similar to getDefault() in aaflib
    =#
    global function getLastAAFIndex()
        return last
    end

     #=
     # Assign new coefficient index i for noise symbol μᵢ.
     # To be used when assigning computation noise.
     # Similar to inclast() in aaflib
    =#
    global function addAAFIndex()
        last += 1
        return [last]
    end

    global function addAAFIndex(indexes::Vector{AAFInd})
        last += 1
        return vcat(indexes,last)
    end

     #=
     # TODO: finish and make sure this is correct
    =#
    global function setLastAAFIndex(indexes::Vector{AAFInd})
        m = last(indexes)
        if(m > last)
            last = m
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
     # TODO: not sure if dev, ind that are longer than t in aaflib is applicable
     # TODO: what if t > length(dev) ?
     # TODO: get rid of l; clean up comments
    =#
    #AAF(v0::AAFCoeff, dev::Vector{AAFCoeff}, ind::Vector{AAFCoeff}, t:Unsigned) = new(
    #    v0, t, t, dev[1:t], setLastAFFIndex(ind[1:t])
    #)
    #function AAF(v0::AAFCoeff, dev::Vector{AAFCoeff}, ind::Vector{AAFInd}, l::Int)
    #    @assert length(ind) == l
    #    @assert length(dev) == l
    #    new(v0, l, l, dev, setLastAAFIndex(ind))
    #end
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

function Base.getindex(a::AAF, ind::Int)::AAFCoeff
    if(ind < 0 || ind > length(a.deviations))
        return 0.0
    elseif(ind == 0)
        return a.cvalue
    else
        return a.deviations[ind]
    end
end

function Base.length(a::AAF)
    return length(a.deviations)
end

function Base.repr(a::AAF)
    s = "$(a[0])"
    if(length(a) > 0)
        for i in 1:length(a)
            s *= " + $(a[i])ϵ$(a.indexes[i])"
        end
    end
    return s
end

Interval(a::AAF) = Interval(a[0] - rad(a), a[0] + rad(a))

 #=
 # Goubault+Putot methods
 # Specification: convert_int, reduce_aaf, rad
 # TODO: implement
=#

 #=
 # Get the total deviation of an AAF (i.e. the sum of all deviations (their abs value))
=#
rad(a::AAF)::AAFCoeff = sum(abs.(a.deviations))

  #=
  # Get maximum / minimum of an AAF
 =#
getMax(a::AAF)::AAFCoeff = a[0] + rad(a)
getMin(a::AAF)::AAFCoeff = a[0] - rad(a)
getAbsMax()::AAFCoeff = max(abs(a[0] - rad(a)), abs(a[0] + rad(a)))
getAbsMin()::AAFCoef  = min(abs(a[0] - rad(a)), abs(a[0] + rad(a)))

Base.firstindex(a::AAF) = length(a) > 0 ? a.indexes[1] : 0
Base.lastindex(a::AAF)  = length(a) > 0 ? last(a.indexes) : 0 

 #=
 # Unknown methods
 # Specification: compact, sumup
 # TODO: implement
=#

  #=
  # Conditionals
 =#
 Base.:<(a::AAF,  p::AAF) = (a[0] + rad(a)) <  (p[0] - rad(p))
 Base.:<=(a::AAF, p::AAF) = (a[0] + rad(a)) <= (p[0] - rad(p))
 Base.:>(a::AAF,  p::AAF) = (a[0] - rad(a)) >  (p[0] + rad(p))
 Base.:>=(a::AAF, p::AAF) = (a[0] - rad(a)) >= (p[0] + rad(p))

 #=
 # Equality
 # Specification: affine equality compares cvalues, and deviations TOL=1E-15
 # TODO: is it better to use ≈ instead?
 # TODO: refactor length out
=#
TOL = 1E-15
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

*(a::AAF, cst::AAFCoeff)::AAF = AAF(a, a[0] * cst, cst * deviations)
*(cst::AAFCoeff, a::AAF)::AAF = AAF(a, cst * a[0], cst * deviations)

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
    pcomp = pair.(indexin(indt, a.indexes), indexin(indt, p.indexes))
    for (ii, (ia,ip)) in enumerate(pcomp)
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
    pcomp = pair.(indexin(indt, a.indexes), indexin(indt, p.indexes))
    for (ii, (ia,ip)) in enumerate(pcomp)
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
 #
 # Specification:
 #
 # TODO: there are various multiplication forms
 #   xy = x₀ŷ₀ + ∑ᴺᵢ(xᵢy₀+yᵢx₀)ϵᵢ + ½∑[over 1⩽i,j⩽n] |xᵢyⱼ+yᵢxⱼ|μₖ
 #   xy = x₀ŷ₀ + ∑ᴺᵢ(xᵢy₀+yᵢx₀)ϵᵢ + (∑ᴺᵢ|xᵢ|)(∑ᴺᵢ|yᵢ|)μₖ (using this one)
 #
 # TODO: aaflib uses an unknown approximation method for μₖ
=#
function *(a::AAF, P::AAF)::AAF
    if(length(p) == 0)
        return a * p[0]
    elseif(length(a) == 0)
        return a[0] * p
    end
    # create new index with length = length(a.indexes) + length(p.indexes) + 1
    indt = sort([ii for ii in union(a.indexes, p.indexes)])
    indt = addAAFIndex(indt)
    lindt = length(indt)
    devt  = fill(0.0, lindt)
    pcomp = pair.(indexin(indt, a.indexes), indexin(indt, p.indexes))
    for (ii, (ia,ip)) in enumerate(pcomp)
        devt[ii] = (ia == nothing ? a[0] * p[ip] : 
                    (ip == nothing ? a[ia] * p[0] : 
                     a[ia] * p[0] + a[0] * p[ip]))
        #devt[lindt] += (ia == nothing || ip == nothing) ? 0 :
        #                    a[ia] * p[ip]
    end
    devt[lindt] = rad(a) * rad(p)
    return AAF(a[0] * p[0], devt, indt)
end

function inv(a::AAF)::AAF
    if(length(a) == 0)
        return AAF(1.0 / a[0])
    end
end

function /(a::AAF, P::AAF)::AAF
    #
end

#function ^(a::AAF) const;
#function ^(const int) const;
#function -(a::AAF) const;
#function *(double);
#function /(double) const;

#AAF operator * (double, const AAF);
#AAF operator / (double, const AAF);
#AAF operator + (double, const AAF);
#AAF operator - (double, const AAF);






