#__precompile__(true)
#
#module IntervalArithmetic

#=

Functions are obtained from Modal interval analysis book.

Sainz, M.A., Armengol, J., Calm, R., Herrero, P., Jorba, L. and Vehi, J., 2014. Modal interval analysis. Lecture Notes in Mathematics, 2091.
=#

# TODO: turn this file into a module, optimize compilation\runtime
# TODO: complete tests for ModalInterval
# TODO: create IntervalBox

using IntervalArithmetic

import Base:
    +, -, *, /, //, fma,
    <, >, ==, !=, ⊆, ^, <=,
    in, zero, one, eps, typemin, typemax, abs, abs2, real, min, max,
    sqrt, exp, log, sin, cos, tan, inv,
    exp2, exp10, log2, log10,
    asin, acos, atan,
    sinh, cosh, tanh, asinh, acosh, atanh, sinpi, cospi,
    union, intersect, isempty,
    convert, promote_rule, eltype, size,
    BigFloat, float, widen, big,
    ∩, ∪, ⊆, ⊇, ∈, eps

function debug(); print("DEBUG\n"); end

# In Modal Intervals, proper ≡ ∃, improper ≡ ∀, and real are both ∀,∃
@enum Predicate proper=1 real_number=0 improper=-1 

# Modal Interval of dimension 1
struct ModalInterval
    prime::Interval
    pred::Predicate

    function ModalInterval(infX::Real,supX::Real)
        pred = infX == supX ? real_number : 
            (infX < supX ? proper : improper)
        new(interval(min(infX,supX), max(infX,supX)), pred)
    end
    
    function ModalInterval(X::Interval)
        pred = X.lo == X.hi ? real_number : proper
        new(interval(X), pred)
    end
end

# Getter functions
mod(X::ModalInterval) = X.pred
prop(X::ModalInterval) = interval(X.prime)

 #=
 # Returns sup(X) = x₂ where X = [x₁, x₂]
=# 
sup(X::ModalInterval) = mod(X) == improper ? X.prime.lo : X.prime.hi

 #=
 # Returns inf(X) = x₁ where X = [x₁, x₂]
=# 
inf(X::ModalInterval) = mod(X) == improper ? X.prime.hi : X.prime.lo

 #=
 # Prints [x₁, x₂] =: X
=# 
display(X::ModalInterval) = print("[$(inf(X)), $(sup(X))]")

function ret(X::ModalInterval)
    return inf(X), sup(X)
end
function dual(X::ModalInterval)
    return ModalInterval(sup(X),inf(X))
end

# Modal inclusion / equality
# In Modal Intervals, proper ≡ ∃, improper ≡ ∀, and real are both ∀,∃
⊆(A::ModalInterval, B::ModalInterval) = inf(A) ≥ inf(B) && sup(A) ≤ sup(B)
⊇(A::ModalInterval, B::ModalInterval) = A ⊆ B
==(A::ModalInterval, B::ModalInterval) = inf(A) == inf(B) && sup(B) == sup(A)

# Arithmetic
+(A::ModalInterval, B::ModalInterval) = ModalInterval(inf(A)+inf(B),sup(A)+sup(B))
-(X::ModalInterval) = ModalInterval(-sup(X),-inf(X))
-(A::ModalInterval, B::ModalInterval) = A + (-B)

# Kaucher multiplication
function *(A::ModalInterval, B::ModalInterval)
    (a₁, a₂) = ret(A)
    (b₁, b₂) = ret(B)
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

function Base.:^(A::ModalInterval, n::Int)
    (a₁, a₂) = ret(A)
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
    (a₁, a₂) = ret(A)
    (b₁, b₂) = ret(B)
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









