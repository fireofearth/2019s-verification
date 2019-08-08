module ModalIntervalArithmetic
#=
Functions are obtained from Modal interval analysis book.

Sainz, M.A., Armengol, J., Calm, R., Herrero, P., Jorba, L. and Vehi, J., 2014. Modal interval analysis. Lecture Notes in Mathematics, 2091.

=#

import Base:
    +, -, *, /, ^,
    <, >, ==, !=, <=,
    issubset,
#    ∩, ∪, ⊆, ⊇, ∈,
    min, max,
	show, convert, promote_rule, isreal,
    mod, abs, abs2, real,
    sqrt, exp, log, sin, cos, tan, inv
#    exp2in, exp10, log2, log10,
#    asin, acos, atan,
#    sinh, cosh, tanh, asinh, acosh, atanh, sinpi, cospi,
#    union, intersect, isempty,
#    eltype, size,
#    BigFloat, float, widen, big,
#    eps

using IntervalArithmetic

import IntervalArithmetic:
    Interval, inf, sup, ⊂, ⊃

export
    ModalInterval, Interval,
    Predicate, PROPER, IMPROPER, REAL_NUMBER,
    mod, inf, sup, show, fields, dual, prop,
	min, max,
	one, zero,
    ==, issubset, ⊂, ⊃
    #+, -, *, /, inv

# In Modal Intervals, proper ≡ ∃, improper ≡ ∀, and real are both ∀,∃
@enum Predicate PROPER=1 REAL_NUMBER=0 IMPROPER=-1 

# Modal Interval of dimension 1
struct ModalInterval{T<:Real} <: Number
    inf::T
    sup::T
end

ModalInterval(a::Real, b::Real) = ModalInterval(promote(a, b)...)
ModalInterval{T}(a::Real) where {T<:Real} = ModalInterval{T}(a, a)
ModalInterval(a::Real) = ModalInterval(a, a)
ModalInterval(X::Interval) = ModalInterval(X.lo, X.hi)
ModalInterval{T}(X::Interval) where {T<:Real} = ModalInterval{T}(X.lo, X.hi)
ModalInterval(X::ModalInterval) = X

# Getter functions
inf(X::ModalInterval) = X.inf
sup(X::ModalInterval) = X.sup
fields(X::ModalInterval) = inf(X), sup(X)
min(X::ModalInterval)    = min(fields(X)...)
max(X::ModalInterval)    = max(fields(X)...)
mod(X::ModalInterval)    = inf(X) == sup(X) ? REAL_NUMBER : 
    (inf(X) < sup(X) ? PROPER : IMPROPER)
prop(X::ModalInterval)   = ModalInterval(min(X), max(X))
dual(X::ModalInterval)   = ModalInterval(sup(X), inf(X))
isreal(X::ModalInterval) = mod(X) == REAL_NUMBER

# Modal inclusion / equality
# In Modal Intervals, proper ≡ ∃, improper ≡ ∀, and real are both ∀,∃
==(A::ModalInterval, B::ModalInterval) = inf(A) == inf(B) && sup(B) == sup(A)
⊂(A::ModalInterval, B::ModalInterval) = inf(A) > inf(B) && sup(A) < sup(B)
⊃(A::ModalInterval, B::ModalInterval) = B ⊂ A
# ⊆ = issubset
issubset(A::ModalInterval, B::ModalInterval) = A ⊂ B || A == B

 #=
 # Obtain the string representation of modal interval
 # String representation is [x₁, x₂] =: X
=# 
show(io::IO, X::ModalInterval) = print(io, "[$(inf(X)), $(sup(X))]")

function Interval(X::ModalInterval)
    if(mod(X) == IMPROPER)
        throw(DomainError(X, "improper modal integrals cannot be converted to proper integrals; use Interval(Dual(X)) instead."))
    end
    Interval(inf(X), sup(X))
end

# Promotion rules
promote_rule(::Type{ModalInterval{T}}, ::Type{ModalInterval{S}}) where {T<:Real, S<:Real} =
    ModalInterval{promote_type(T, S)}
promote_rule(::Type{ModalInterval{T}}, ::Type{S}) where {T<:Real, S<:Real} =
    ModalInterval{promote_type(T, S)}

# Conversion
convert(::Type{ModalInterval{T}}, x::S) where {T<:Real, S<:Real} = 
    ModalInterval{T}(x)
convert(::Type{ModalInterval}, x::S) where {S<:Real} = 
    ModalInterval{S}(x)
convert(::Type{ModalInterval{T}}, X::Interval{S}) where {T<:Real, S<:Real} = 
    ModalInterval{T}(X)
convert(::Type{ModalInterval}, X::Interval{S}) where {S<:Real} = 
    ModalInterval{S}(X)

# Arithmetic
+(A::ModalInterval, B::ModalInterval) = ModalInterval(inf(A)+inf(B),sup(A)+sup(B))
-(X::ModalInterval) = ModalInterval(-sup(X),-inf(X))
-(A::ModalInterval, B::ModalInterval) = A + (-B)

# Kaucher multiplication
function *(A::ModalInterval, B::ModalInterval)
    (a₁, a₂) = fields(A)
    (b₁, b₂) = fields(B)
    if(a₁ ≥ 0 && a₂ ≥ 0 && b₁ ≥ 0 && b₂ ≥ 0)
        return ModalInterval(a₁*b₁, a₂*b₂)
    elseif(a₁ ≥ 0 && a₂ ≥ 0 && b₁ ≥ 0 && b₂ < 0)
        return ModalInterval(a₁*b₁, a₁*b₂)
    elseif(a₁ ≥ 0 && a₂ ≥ 0 && b₁ < 0 && b₂ ≥ 0)
        return ModalInterval(a₂*b₁, a₂*b₂)
    elseif(a₁ ≥ 0 && a₂ ≥ 0 && b₁ < 0 && b₂ < 0)
        return ModalInterval(a₂*b₁, a₁*b₂)
    elseif(a₁ ≥ 0 && a₂ < 0 && b₁ ≥ 0 && b₂ ≥ 0)
        return ModalInterval(a₁*b₁, a₂*b₁)
    elseif(a₁ ≥ 0 && a₂ < 0 && b₁ ≥ 0 && b₂ < 0)
        return ModalInterval(max(a₂*b₂, a₁*b₁), min(a₂*b₁, a₁*b₂))
    elseif(a₁ ≥ 0 && a₂ < 0 && b₁ < 0 && b₂ ≥ 0)
        return ModalInterval(0, 0)
    elseif(a₁ ≥ 0 && a₂ < 0 && b₁ < 0 && b₂ < 0)
        return ModalInterval(a₂*b₂, a₁*b₂)
    elseif(a₁ < 0 && a₂ ≥ 0 && b₁ ≥ 0 && b₂ ≥ 0)
        return ModalInterval(a₁*b₂, a₂*b₂)
    elseif(a₁ < 0 && a₂ ≥ 0 && b₁ ≥ 0 && b₂ < 0)
        return ModalInterval(0, 0)
    elseif(a₁ < 0 && a₂ ≥ 0 && b₁ < 0 && b₂ ≥ 0)
        return ModalInterval(min(a₁*b₂, a₂*b₁), max(a₁*b₁, a₂*b₂))
    elseif(a₁ < 0 && a₂ ≥ 0 && b₁ < 0 && b₂ < 0)
        return ModalInterval(a₂*b₁, a₁*b₁)
    elseif(a₁ < 0 && a₂ < 0 && b₁ ≥ 0 && b₂ ≥ 0)
        return ModalInterval(a₁*b₂, a₂*b₁)
    elseif(a₁ < 0 && a₂ < 0 && b₁ ≥ 0 && b₂ < 0)
        return ModalInterval(a₂*b₂, a₂*b₁)
    elseif(a₁ < 0 && a₂ < 0 && b₁ < 0 && b₂ ≥ 0)
        return ModalInterval(a₁*b₂, a₁*b₁)
    else # if(a₁ < 0 && a₂ < 0 && b₁ < 0 && b₂ < 0)
        return ModalInterval(a₂*b₂, a₁*b₁)
    end
end

function ^(A::ModalInterval, n::Int)
    (a₁, a₂) = fields(A)
    if(isodd(n))
        return ModalInterval(a₁^n, a₂^n)
    else
        if(a₁ ≥ 0, x₂ ≥ 0)
            return ModalInterval(a₁^n, a₂^n)
        elseif(a₁ < 0, x₂ < 0)
            return ModalInterval(a₂^n, a₁^n)
        elseif(a₁ < 0, x₂ ≥ 0)
            return ModalInterval(0, max(abs(a₂)^n, abs(a₁)^n))
        else # if(a₁ ≥ 0, x₂ < 0)
            return ModalInterval(max(abs(a₂)^n, abs(a₁)^n), 0)
        end
    end
end

function inv(B::ModalInterval)
    if(0 ∈ prop(B))
        throw(DomainError(B, "proper interval must not contain zero"))
    end
    return ModalInterval(1/sup(B), 1/inf(B))
end

 #=
 # Kaucher division
=#
function /(A::ModalInterval, B::ModalInterval)
    if(0 ∈ prop(B))
        throw(DomainError(B, "proper interval must not contain zero"))
    end
    (a₁, a₂) = fields(A)
    (b₁, b₂) = fields(B)
    if(a₁ ≥ 0 && a₂ ≥ 0 && b₁ > 0 && b₂ > 0)
        return ModalInterval(a₁/b₂, a₂/b₁)
    elseif(a₁ ≥ 0 && a₂ ≥ 0 && b₁ < 0 && b₂ < 0)
        return ModalInterval(a₂/b₂, a₁/b₁)
    elseif(a₁ ≥ 0 && a₂ < 0 && b₁ > 0 && b₂ > 0)
        return ModalInterval(a₁/b₂, a₂/b₂)
    elseif(a₁ ≥ 0 && a₂ < 0 && b₁ < 0 && b₂ < 0)
        return ModalInterval(a₂/b₁, a₁/b₁)
    elseif(a₁ < 0 && a₂ ≥ 0 && b₁ > 0 && b₂ > 0)
        return ModalInterval(a₁/b₁, a₂/b₁)
    elseif(a₁ < 0 && a₂ ≥ 0 && b₁ < 0 && b₂ < 0)
        return ModalInterval(a₂/b₂, a₁/b₂)
    elseif(a₁ < 0 && a₂ < 0 && b₁ > 0 && b₂ > 0)
        return ModalInterval(a₁/b₁, a₂/b₂)
    else # if(a₁ < 0 && a₂ < 0 && b₁ < 0 && b₂ < 0)
        return ModalInterval(a₂/b₁, a₁/b₂)
    end
end

#===================#
# temporary code
import AffineArithmetic: Affine

ModalInterval(a::Affine) = ModalInterval(min(a), max(a))
#===================#

end # module ModalIntervalArithmetic
