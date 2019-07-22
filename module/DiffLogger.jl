module DiffLogger

import Base:
    zero, one, iszero, isone, convert, isapprox, promote_rule, repr,
    ==, +, -, *, /, inv, ^

using Logging

export
    zero, one, iszero, isone, convert, isapprox, promote_rule, repr,
    ==, +, -, *, /, inv, ^,
    TestReal

struct TestReal <: Real
    x::Real

    TestReal(x::Real) = new(x)
    TestReal(a::TestReal) = new(a.x)
end

convert(::Type{TestReal}, x::TestReal) = TestReal(x)
convert(::Type{TestReal}, x::Real) = TestReal(x)

one(::Type{TestReal})  = convert(TestReal, 1)
one(x::TestReal)       = convert(TestReal, 1)
zero(::Type{TestReal}) = convert(TestReal, 0)
zero(x::TestReal)      = convert(TestReal, 0)

isone(x::TestReal)     = x == one(Affine)
iszero(x::TestReal)    = x == zero(Affine)

promote_rule(::Type{TestReal}, ::Type{Real}) = TestReal
promote_rule(::Type{Real}, ::Type{TestReal}) = TestReal

==(a::TestReal, p::TestReal) = a.x == p.x

function +(a::TestReal, b::TestReal)
    @info "in +(a::TestReal, b::TestReal)\n $(a.x) + $(b.x) = $(a.x + b.x)"
    TestReal(a.x + b.x)
end

function +(a::TestReal, x::Real)
    @info "in +(a::TestReal, x::Real)\n $(a.x) + $(x) = $(a.x + x)"
    TestReal(a.x + x)
end

function +(x::Real, b::TestReal)
    @info "in +(x::Real, b::TestReal)\n $(x) + $(b.x) = $(x + b.x)"
    TestReal(x + b.x)
end

function -(a::TestReal, b::TestReal)
    @info "in -(a::TestReal, b::TestReal)\n $(a.x) - $(b.x) = $(a.x - b.x)"
    TestReal(a.x - b.x)
end

function -(a::TestReal, x::Real)
    @info "in -(a::TestReal, x::Real)\n $(a.x) - $(x) = $(a.x - x)"
    TestReal(a.x - x)
end

function -(x::Real, b::TestReal)
    @info "in -(a::TestReal, b::TestReal)\n $(x) - $(b.x) = $(x - b.x)"
    TestReal(x - b.x)
end

function *(a::TestReal, b::TestReal)
    @info "in *(a::TestReal, b::TestReal)\n $(a.x) * $(b.x) = $(a.x * b.x)"
    TestReal(a.x * b.x)
end

function *(a::TestReal, x::Real)
    @info "in *(a::TestReal, x::Real)\n $(a.x) * $(x) = $(a.x * x)"
    TestReal(a.x * x)
end

function *(x::Real, b::TestReal)
    @info "in *(x::Real, b::TestReal)\n $(x) * $(b.x) = $(x * b.x)"
    TestReal(x * b.x)
end

function /(a::TestReal, b::TestReal)
    @info "in /(a::TestReal, b::TestReal)\n $(a.x) / $(b.x) = $(a.x / b.x)"
    TestReal(a.x / b.x)
end

function /(a::TestReal, x::Real)
    @info "in /(a::TestReal, x::Real)\n $(a.x) / $(x) = $(a.x / x)"
    TestReal(a.x / x)
end

function /(x::Real, b::TestReal)
    @info "in /(x::Real, b::TestReal)\n $(x) / $(b.x) = $(x / b.x)"
    TestReal(x / b.x)
end

function +(a::TestReal)
    @info "in +(a::TestReal)\n +($(a.x)) = $(+a.x)"
    TestReal(+a.x)
end

function -(a::TestReal)
    @info "in -(a::TestReal)\n -($(a.x)) = $(-a.x)"
    TestReal(-a.x)
end

function inv(a::TestReal)
    @info "in inv(a::TestReal)\n inv($(a.x)) = $(inv(a.x))"
    TestReal(inv(a.x))
end

function ^(a::TestReal, n::Int)
    @info "in ^(a::TestReal, n::Int)\n $(a.x)^$(n) = $(a.x^n)"
    TestReal(a.x^n)
end

end
