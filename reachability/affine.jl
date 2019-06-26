
 #=
 # Affine arithmetic module
 #
 # This is a implementation of the An Affine Arithmetic C++ Library ( http://aaflib.sourceforge.net ) in Julia
 #
 # TODO: aaflib contains 
=#

using IntervalArithmetic, StaticArrays

import Base:
    getindex

@enum tApproximationType MINRANGE CHEBYSHEV SECANT

struct AAF
    cvalue::Float64  # central value 
    length::Unsigned # length of indexes 
    size::Unsigned   # array size of indexes and deviations

    deviations::Array{Float64}
    indexes::Array{Unsigned}
    
    function AAF(v0::double = 0.0)
    end
end

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

function getindex(a::AAF, i::Unsigned)::Float64
    #
end


function getindex(a::AAF, i::Unsigned)::Float64
    #
end

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

function =(a::AAF, P::AAF)::AAF
    #
end

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

AAF operator * (double, const AAF);
AAF operator / (double, const AAF);
AAF operator + (double, const AAF);
AAF operator - (double, const AAF);






