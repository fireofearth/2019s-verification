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

import Base

 #=
 # Type declarations
 #
 #
=#
AAFCoeff = Float64

@enum tApproximationType MINRANGE CHEBYSHEV SECANT

 #=
 # last keeps record of last coefficient index of affine forms accoss all AAF instances.
 #
 # Specification:
 # - force setLastAFFIndex to call whenever a new AAF instance is created.
 #
 # TODO: turn this into a decorator/macro and force calls
=#
let last::Unsigned = 0

    global function resetLastAFFIndex()
        last = 0
    end

    global function getLastAFFIndex()
        return last
    end

     #=
     # Assign new coefficient index i for noise symbol μᵢ.
     # To be used when assigning computation noise.
    =#
    global function setLastAAFIndex()::Bool
        last += 1
        return [last]
    end

     #=
     # TODO:
    =#
    global function setLastAAFIndex(indexes::Vector{Unsigned})::Bool
        if(pop(indexes) > last)
            last = pop(indexes)
            return true
        else
            return false
        end
    end
end

 #=
 # AAF represents an affine form
 #
 # Specification:
 # - AAF is a collection with deviations as its elements.
 #
 # TODO: it's not clear whether `length` and `size` are redundant
 # TODO: should be able to get rid of length, can query vector length
=#
struct AAF
    cvalue::AAFCoeff  # central value 
    length::Unsigned # length of indexes 
    size::Unsigned   # array size of indexes and deviations
    deviations::Vector{AAFCoeff}
    indexes::Vector{Unsigned}

     #=
     # Creates an AAF without deviations.
    =#
    function AAF(v0::AAFCoeff = 0.0)
        new(v0, 0, 0, Vector{AAFCoeff}(), Vector{AAFCoeff}())
    end
    
     #=
     # Creates an AAF with deviations
     # TODO: not sure if dev, ind that are longer than t in aaflib is applicable
    =#
    function AAF(v0::AAFCoeff, dev::Vector{AAFCoeff}, ind::Vector{AAFCoeff}, t:Unsigned )
        setLastAFFIndex(ind[1:t])
        new(v0, t, t, dev[1:t], ind[1:t])
    end
    
     #=
     # Constructor from intervals
    =#
    function AAF(iv::Interval)
        cvalue = mid(iv)
        new(cvalue, 1, 1, [radius(iv)], setLastAAFIndex())
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
# last::Unsigned

function getClassSize(a::AAF)::Unsigned
    #
end


function getClassSizeIndex(a::AAF)::Unsigned
    #
end

function getClassSizeDeviations(a::AAF)::Unsigned
    #
end

function Base.getindex(a::AAF, i::Unsigned)::AAFCoeff
    #
end

function Base.length(a::AAF)
    return a.length
end

function getCenter(a::AAF)
    #
end

function convert(a::AAF)
    #
end

 #=
 # Goubault+Putot methods
 #
 # Specification: convert_int, reduce_aaf, rad
 # 
 # TODO: implement
=#

function rad(a::AAF)::AAFCoeff
    #
end

function getMax(a::AA

function <(a::AAF, P::AAF)::Bool
    #
end

function <=(a::AAF, P::AAF)::Bool
    #
end

function  >(a::AAF, P::AAF)::Bool
    #
end

function >=(a::AAF, P::AAF)::Bool
    #
end

function ==(a::AAF, P::AAF)::Bool
    #
end

# function =(a::Float64)

#function =(a::AAF, P::AAF)::AAF
#end

function +(a::AAF, P::AAF)::AAF
    #
end

function -(a::AAF, P::AAF)::AAF
    #
end

function *(a::AAF, P::AAF)::AAF
    #
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






