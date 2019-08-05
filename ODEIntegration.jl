using LinearAlgebra
using ForwardDiff
using AffineArithmetic

 #=
 # Construct Taylor Model
 #
 # Specification:
 # - Given problem ̇z' = f(z), we construct the function
 # [z](t, tⱼ, [zⱼ]) = [zⱼ] + ∑{i=1,…,k-1} (t-tⱼ)ⁱ/i! f⁽ⁱ⁾([zⱼ]) + (t-tⱼ)ᵏ/k! f⁽ᵏ⁾([rⱼ₊₁])
 # used by the HSCC'17 article by Goubault+Putot and returns it
 #
 # - Assumes that f: Rᴺ → Rᴺ
 #
 # - We do the calculation with affine forms as initial conditions, this implies 
 # [zₒ] ≡ Affine(center([zₒ])), etc... However we permit evaluating for any type T for zⱼ
 #
 # TODO:
 # - ForwardDiff fails when applying Jacobian >7 times, so we require order ≤ 7
 # and we may need to use TaylorSeries after all
 # - May want to collapse constructTM{Affine, Real} to one function
=#

#=
@generated function constructTM(f::Function; order::Int=5, T::Type=Affine)
    # store the lie derivatives f⁽ⁱ⁾ (i = 1,…,order) in vf
    vf = [f]
    for i in 2:order
        fi = (z::Vector -> ForwardDiff.jacobian(vf[i - 1], z) * f(z))
        vf = vcat(vf, fi)
    end
end
=#

function constructTMAffine(f::Function, order::Int)
    # store the lie derivatives f⁽ⁱ⁾ (i = 1,…,order) in vf
    vf = [f]
    for i in 2:order
        fi = (z::Vector -> ForwardDiff.jacobian(vf[i - 1], z) * f(z))
        vf = vcat(vf, fi)
    end

    # constructed TM
    # cTM(t,tⱼ,[zⱼ],[rⱼ₊₁]) = [zⱼ] + ∑ᵢ(t-tⱼ)ⁱ/i! f⁽ⁱ⁾([zⱼ]) + (t-tⱼ)ᵏ/k! f⁽ᵏ⁾([rⱼ₊₁])
    function cTM(t::Real, tj::Real, zj::Vector{Affine}, r::Vector{Affine})
        acc = zj
        for i in 1:(order - 1)
            term = 1.0
            for l in 1:i
                term *= (t - tj) / l
            end
            acc += term * vf[i](zj)
        end
        term = 1.0
        for l in 1:order
            term *= (t - tj) / l
        end
        acc += term * vf[order](r)
        return acc
    end

    return cTM
end

function constructTMReal(f::Function, order::Int)
    # store the lie derivatives f⁽ⁱ⁾ (i = 1,…,order) in vf
    vf = [f]
    for i in 2:order
        fi = (z::Vector{<:Real} -> ForwardDiff.jacobian(vf[i - 1], z) * f(z))
        vf = vcat(vf, fi)
    end

    # constructed TM
    # cTM(t,tⱼ,[zⱼ],[rⱼ₊₁]) = [zⱼ] + ∑ᵢ(t-tⱼ)ⁱ/i! f⁽ⁱ⁾([zⱼ]) + (t-tⱼ)ᵏ/k! f⁽ᵏ⁾([rⱼ₊₁])
    function cTM(t::Real, tj::Real, zj::Vector{<:Real}, r::Vector{<:Real})
        acc = zj
        for i in 1:(order - 1)
            term = 1.0
            for l in 1:i
                term *= (t - tj) / l
            end
            acc += term * vf[i](zj)
        end
        term = 1.0
        for l in 1:order
            term *= (t - tj) / l
        end
        acc += term * vf[order](r)
        return acc
    end

    return cTM
end

function constructTM(f::Function; order::Int=5, T::Type=Affine)
    if(T == Affine)
        return constructTMAffine(f, order)
    elseif(T == Real)
        return constructTMReal(f, order)
    else
        error("constructTM: T is not real or an affine form")
    end
end

 #=
 # Construct Taylor Model of the Jacobian
 #
 # Specification:
 # - Given f, we construct the function
 # [J](t, tⱼ, [zⱼ]) 
 #  = [Jⱼ] + ∑{i=1,…,k-1} (t-tⱼ)ⁱ/i! Jac f⁽ⁱ⁾([zⱼ]) [Jⱼ] + (t-tⱼ)ᵏ/k! Jac f⁽ᵏ⁾([rⱼ₊₁]) [Rⱼ₊₁]
 #
 # - Assumes that f: Rᴺ → Rᴺ
 #
 # - We do the calculation with affine forms as initial conditions, this implies 
 # [zₒ] ≡ Affine(center([zₒ])), etc... However we permit evaluating for any type T for zⱼ
 #
 # TODO:
 # - possible overflow when calculating fartorials, orders
=#
function constructJacTMAffine(f::Function, order::Int)

    # store the lie derivatives f⁽ⁱ⁾ (i = 1,…,order) in vf
    vf = [f]
    for i in 2:order
        fi = (z::Vector{Affine} -> ForwardDiff.jacobian(vf[i - 1], z) * f(z))
        vf = vcat(vf, fi)
    end

    # store the jacobians of the lie derivatives Jac(f⁽ⁱ⁾) in vJacf
    vJacf = [ ]
    for i in 1:order
        Jacfi = (z::Vector{Affine} -> ForwardDiff.jacobian(vf[i], z))
        vJacf = vcat(vJacf, Jacfi)
    end

    # cJacTM(t,tⱼ,[zⱼ],[Jⱼ],[rⱼ₊₁],[Rⱼ₊₁])
    # = [Jⱼ] + ∑{i=1,…,k-1} (t-tⱼ)ⁱ/i! Jac f⁽ⁱ⁾([zⱼ]) [Jⱼ] + (t-tⱼ)ᵏ/k! Jac f⁽ᵏ⁾([rⱼ₊₁]) [Rⱼ₊₁]
    function cJacTM(t::Real, tj::Real, zj::Vector{Affine}, 
                    Jj::Matrix{Affine}, r::Vector{Affine}, R::Matrix{Affine})
        acc = Jj
        for i in 1:(order - 1)
            acc += ((t - tj)^i / factorial(i)) * vJacf[i](zj) * Jj
        end
        acc += ((t - tj)^order / factorial(order)) * vJacf[order](r) * R
        return acc
    end

    return cJacTM
end

function constructJacTMReal(f::Function, order::Int)

    # store the lie derivatives f⁽ⁱ⁾ (i = 1,…,order) in vf
    vf = [f]
    for i in 2:order
        fi = (z::Vector{<:Real} -> ForwardDiff.jacobian(vf[i - 1], z) * f(z))
        vf = vcat(vf, fi)
    end

    # store the jacobians of the lie derivatives Jac(f⁽ⁱ⁾) in vJacf
    vJacf = [ ]
    for i in 1:order
        Jacfi = (z::Vector{<:Real} -> ForwardDiff.jacobian(vf[i], z))
        vJacf = vcat(vJacf, Jacfi)
    end

    # cJacTM(t,tⱼ,[zⱼ],[Jⱼ],[rⱼ₊₁],[Rⱼ₊₁])
    # = [Jⱼ] + ∑{i=1,…,k-1} (t-tⱼ)ⁱ/i! Jac f⁽ⁱ⁾([zⱼ]) [Jⱼ] + (t-tⱼ)ᵏ/k! Jac f⁽ᵏ⁾([rⱼ₊₁]) [Rⱼ₊₁]
    function cJacTM(t::Real, tj::Real, zj::Vector{<:Real}, 
                    Jj::Matrix{<:Real}, r::Vector{<:Real}, R::Matrix{<:Real})
        acc = Jj
        for i in 1:(order - 1)
            term = 1.0
            for l in 1:i
                term *= (t - tj) / l
            end
            acc += term * vJacf[i](zj) * Jj
        end
        term = 1.0
        for l in 1:order
            term *=  (t - tj) / l
        end
        acc += term * vJacf[order](r) * R
        return acc
    end

    return cJacTM
end

function constructJacTM(f::Function; order::Int=5, T::Type=Affine)
    if(T == Affine)
        return constructJacTMAffine(f, order)
    elseif(T == Real)
        return constructJacTMReal(f, order)
    else
        error("constructTM: T is not real or an affine form")
    end
end

 #=
 # Compute a priori enclosure for [rⱼ₊₁]
 #
 # Specification:
 # - Compute a priori enclosure [rⱼ₊₁] of f using the interval extension of the 
 # Picard-Lindelof operator.
 # The Picard Lindelof is expressed as F[z(t)] = z₀ + ∫{tⱼ to t} f(z(s)) ds for t ∈ [tⱼ, tⱼ₊₁]
 # With the interval extension F[z] = z₀ + (tⱼ₊₁ - tⱼ)[f]([z]) with [f] being the natural
 # interval extension of f. Both versions of F admits a unique fixed point if f is Lipschitz.
=#

function fixedPoint(f::Function, z0::Vector{<:Interval}, τ::Real)
    iter = 1
    t    = Interval(0, τ)
    x    = Interval(-1,  1)
    zi   = z0
    Fzi  = z0 + t*f(z0)

    while(iter ≤ 1 || reduce(&, Fzi .⊈ zi))
        @assert iter < 50
        disp("$(iter) $(repr(Interval(zi[1]))) $(repr(Interval(zi[2])))")

        if(iter > 25)
            β = 1.0
        elseif(iter > 20)
            β = 0.1
        elseif(iter > 15)
            β = 0.01
        elseif(iter > 10)
            β = 0.001
        elseif(iter > 5)
            β = 0.0001
        else
            β = 0.00001
        end

        if(iter > 2)
            zi = Fzi + β*x*Fzi
        else
            zi = Fzi
        end
        Fzi = z0 + t*f(zi)

        iter += 1
    end

    return Affine.(zi)
end

fixedPoint(f::Function, z0::Vector{Affine}, τ::Real) = fixedPoint(f, Interval.(z0), τ)

 #=
 # Computes the inner and outer approximations of the flowpipes formed by ODE z' = f(z)
 #
 # Specifications:
 # - initializes variables
 #
 # TODO: finish procedure
=#
#function solveODE(f::Function, tspan::NTuple{2,<:Real}, τ<:Real,
#                  inputs::Vector{<:Interval}; order::Int=4)
#
#    # initialize variables
#    Jⱼ      = Matrix{Affine}(I, 2, 2)
#    zⱼ      = Affine.(inputs)
#    ̃z₀ⱼ     = Affine.(mid.(inputs))
#    tⱼ      = tspan[0]
#    tₙ      = tspan[1]
#    T       = constructTM(f; order=order)
#    JacT    = constructJacTM(f; order=order)
#    # preallocate array to store tⱼ, zⱼ, iiⱼ
#
#    while(tⱼ < tₙ)
#        r     = fixpoint(f, zⱼ, τ)
#        zⱼ    = T(tⱼ + τ, tⱼ, zⱼ, r)
#        r₀    = fixpoint(f, z₀ⱼ, τ)
#        z₀ⱼ   = T(tⱼ + τ, tⱼ, z₀ⱼ, r₀)
#        # fixpoint
#        R     = NaN
#        Jⱼ    = JacT(tⱼ + τ, tⱼ, zⱼ, Jⱼ, r, R)
#        # compute inner approximation
#        iiⱼ   = NaN
#        # save the outer approx. zⱼ and inner approx. iiⱼ
#        tⱼ =+ τ
#    end
#    
#    # get outputs
#end
