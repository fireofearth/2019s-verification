#__precompile__(true)
#
#module IntervalArithmetic

#=

Functions are obtained from Modal interval analysis book.

Sainz, M.A., Armengol, J., Calm, R., Herrero, P., Jorba, L. and Vehi, J., 2014. Modal interval analysis. Lecture Notes in Mathematics, 2091.
=#

# TODO: turn this file into a module, optimize compilation\runtime
# TODO: finish ModalInterval
# TODO: create IntervalBox

using IntervalArithmetic, StaticArrays

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
sup(X::ModalInterval) = mod(X) == improper ? X.prime.lo : X.prime.hi
inf(X::ModalInterval) = mod(X) == improper ? X.prime.hi : X.prime.lo
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
# Obtained
function *(A::ModalInterval, B::ModalInterval)
    (a₁, a₂) = ret(A)
    (b₁, b₂) = ret(B)
    if(a₁ ≥ 0 && a₂ ≥ 0 && b₁ ≥ 0 && b₂ ≥ 0)
        return ModalInterval(a₁*b₁, a₂*b₂)
    elseif(a₁ ≥ 0 && a₂ ≥ 0 && b₁ ≥ 0 && b₂ < 0)
        return ModalInterval(a₁*b₁, a₁*b₂)
    elseif(a₁ ≥ 0 && a₂ ≥ 0 && b₁ < 0 && b₂ ≥ 0)
        return ModalInterval(a₂*b₁, a₂*b₂)
    else
        error("function * is incomplete")
    end
end


