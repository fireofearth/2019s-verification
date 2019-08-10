using LinearAlgebra
using IntervalArithmetic
using AffineArithmetic
using ModalIntervalArithmetic
using ForwardDiff

#=
Flowpipes module

Currently contains solver for inner and outer reachability approximations of flowpipes.

TODO:
- turn this into a module.
=#

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
 # - No longer need T
=#

function constructTM(f::Function; order::Int=5, T::Type=Affine)
    # store the lie derivatives f⁽ⁱ⁾ (i = 1,…,order) in vf
    vf = [f]
    for i in 2:order
        fi = (z::Vector -> ForwardDiff.jacobian(vf[i - 1], z) * f(z))
        vf = vcat(vf, fi)
    end

    # construct TM and methods for (Affine, Real, etc)
    # cTM(t,tⱼ,[zⱼ],[rⱼ₊₁]) = [zⱼ] + ∑ᵢ(t-tⱼ)ⁱ/i! f⁽ⁱ⁾([zⱼ]) + (t-tⱼ)ᵏ/k! f⁽ᵏ⁾([rⱼ₊₁])
    for T in (:Affine, :(<:Real))
        @eval function cTM(t::Real, tj::Real, zj::Vector{$T}, r::Vector{$T})
            acc = zj
            for i in 1:($order - 1)
                term = 1.0
                for l in 1:i
                    term *= (t - tj) / l
                end
                acc += term * $vf[i](zj)
            end
            term = 1.0
            for l in 1:$order
                term *= (t - tj) / l
            end
            acc += term * $vf[$order](r)
            return acc
        end
    end

    return cTM
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
=#
function constructJacTM(f::Function; order::Int=5)

    # store the lie derivatives f⁽ⁱ⁾ (i = 1,…,order) in vf
    vf = [f]
    for i in 2:order
        fi = (z::Vector -> ForwardDiff.jacobian(vf[i - 1], z) * f(z))
        vf = vcat(vf, fi)
    end

    # store the jacobians of the lie derivatives Jac(f⁽ⁱ⁾) in vJacf
    vJacf = [ ]
    for i in 1:order
        Jacfi = (z::Vector -> ForwardDiff.jacobian(vf[i], z))
        vJacf = vcat(vJacf, Jacfi)
    end

    # construct jacobian TM and methods for (Affine, Real, etc)
    # cJacTM(t,tⱼ,[zⱼ],[Jⱼ],[rⱼ₊₁],[Rⱼ₊₁])
    # = [Jⱼ] + ∑{i=1,…,k-1} (t-tⱼ)ⁱ/i! Jac f⁽ⁱ⁾([zⱼ]) [Jⱼ] + (t-tⱼ)ᵏ/k! Jac f⁽ᵏ⁾([rⱼ₊₁]) [Rⱼ₊₁]
    for T in (:Affine, :(<:Real))
        @eval function cJacTM(t::Real, tj::Real, zj::Vector{$T}, 
                        Jj::Matrix{$T}, r::Vector{$T}, R::Matrix{$T})
            acc = Jj
            for i in 1:($order - 1)
                term = 1.0
                for l in 1:i
                    term *= (t - tj) / l
                end
                acc += term * $vJacf[i](zj) * Jj
            end
            term = 1.0
            for l in 1:$order
                term *=  (t - tj) / l
            end
            acc += term * $vJacf[$order](r) * R
            return acc
        end
    end

    return cJacTM
end

 #=
 # Compute a priori enclosure for [rⱼ₊₁]
 #
 # Specification:
 # - Compute a priori enclosure [rⱼ₊₁] of f using the interval extension of the 
 # Picard-Lindelof operator.
 # The Picard Lindelof is expressed as F[z(t)] = z₀ + ∫{tⱼ to t} f(z(s)) ds for t ∈ [tⱼ, tⱼ₊₁]
 # With the interval extension F[z] = z₀ + [0, τ] [f]([z]) with [f] being the natural
 # interval extension of f. Both versions of F admits a unique fixed point if f is Lipschitz.
 #
 # TODO: convert intervals to affines and see whether there is a difference
=#

function fixedPoint(f::Function, z0::Vector{Affine}, τ::Real)
    iter = 1
    t    = Interval(0, τ)
    x    = Interval(-1,  1)
    zi   = z0
    Fzi  = z0 + t*f(z0)

    #while(iter ≤ 1 || reduce(&, Fzi .⊈ zi))
    while(iter ≤ 1 || reduce(|, Fzi .⊈ zi))
        @assert iter < 50
        #disp("$(iter) $(repr(Interval(zi[1]))) $(repr(Interval(zi[2])))")

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

fixedPoint(f::Function, z0::Vector{<:Interval}, τ::Real) = fixedPoint(f, Affine.(z0), τ)

 #=
 # Compute a priori enclosure for [Rⱼ₊₁]
 #
 # TODO: not clear about implementation

fixedPoint uses
L[z] = J₀ + [0, τ] [J∘f]([z])

whereas fixedJacPoint uses
[Jᵣ] = Jf([r])
F[J] = J₀ + [0, τ] [Jᵣ] [J]
=#
function fixedJacPoint(f::Function, J₀::Matrix{<:Interval}, r::Vector{Affine}, τ::Real)
    # RINO:
    #Jac1_g_rough[j][k] = odeVAR_g.x[j][1].d(k);
    #fixpoint(J_rough, Jac1_g_rough, J, tau); output = J_rough
    #void fixpoint(vector<vector<AAF>> &y0, vector<vector<AAF>> &Jac1_g_rough, vector<vector<AAF>> &J0, double tau)
    #multMiMi(fJ0, Jac1_g_rough, y0); // fJ0 = Jac1_g_rough * y0
    #J1[i][j] = J0[i][j] + interval(0, tau) * fJ0[i][j].convert_int();
    iter = 1
    t    = Interval(0, τ)
    x    = Interval(-1,  1)
    Jᵢ   = J₀
    Jᵣ = ForwardDiff.jacobian(f, r)
    Jᵣ = Interval.(Jᵣ)
    FJᵢ = J₀ + t*Jᵣ*J₀

    while(iter ≤ 1 || reduce(|, FJᵢ .⊈ Jᵢ))
        @assert iter < 50
        #disp("$(iter) $(repr(Interval(Jᵢ[1,1]))) $(repr(Interval(Jᵢ[1,2])))")

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
            Jᵢ = FJᵢ + β*x*FJᵢ
        else
            Jᵢ = FJᵢ
        end
        FJᵢ = J₀ + t*Jᵣ*Jᵢ

        iter += 1
    end

    return Affine.(Jᵢ)
end

fixedJacPoint(f::Function, J₀::Matrix{Affine}, 
              r::Vector{Affine}, τ::Real) = fixedJacPoint(f, Interval.(J₀), r, τ)

function innerApproximate(z₀ⱼ::Vector{T}, Jⱼ::Matrix{T}, z₀::Vector{T}) where 
        T <: ModalInterval
    zt₀ = mid.(z₀)
    ia = z₀ⱼ + Jⱼ*(dual.(z₀) - zt₀)
    if(isimproper(ia))
        return Interval.(prop.(ia))
    else
        return fill(NaN, length(z₀))
    end
end

innerApproximate(z₀ⱼ::Vector{Affine}, Jⱼ::Matrix{Affine}, 
                 z₀::Vector{T}) where T <: Interval = 
    innerApproximate(ModalInterval.(z₀ⱼ), ModalInterval.(Jⱼ), ModalInterval.(z₀))

 #=
 # Computes the inner and outer approximations of the flowpipes formed by ODE z' = f(z)
 #
 # Arguments:
 # - f:       function to ODE z' = f(z)
 # - tspan:   tuple containing time interval (t₀, tₙ) (t₀ start time; tₙ end time)
 # - τ:       the amount of time for each time step
 #            fixedPoint, fixedJacPoint *usually* converges when τ ≤ 0.05
 # - z₀:      interval bound on initial conditions
 # - [order]: order of taylor approximation
 #
 # Output: a tuple containing in order
 # - st:  set of 
 # - sz:  s
 # - sia: set intervals containing of inner approximationss, or NaN
 #
 # Specifications:
 # - initializes variables
 #
 # TODO: preallocate arrays st, sz, sia
=#
function approximate(f::Function, tspan::NTuple{2,<:Real}, τ::Real,
                  z₀::Vector{<:Interval}; order::Int=4)

    # initialize variables
    iter = 1
    Jⱼ      = Matrix{Affine}(I, 2, 2)
    tⱼ      = tspan[1]
    # flowpipe outer appoximating ODE with inputs as initial conditions
    # at each point tⱼ
    zⱼ      = Affine.(z₀)
    # forms approximating the center of inputs at time tⱼ
    z₀ⱼ     = Affine.(mid.(z₀))
    tₙ      = tspan[2]
    T       = constructTM(f; order=order)
    JacT    = constructJacTM(f; order=order)
    # preallocate array to store tⱼ, zⱼ, iaⱼ
    st  = [tⱼ]
    sz  = z₀
    sia = fill(NaN, length(z₀))

    while(tⱼ < tₙ)
        disp("iteration = $(iter); tⱼ = $(tⱼ)")
        r     = fixedPoint(f, zⱼ, τ)
        zⱼ    = T(tⱼ + τ, tⱼ, zⱼ, r)
        r₀    = fixedPoint(f, z₀ⱼ, τ)
        z₀ⱼ   = T(tⱼ + τ, tⱼ, z₀ⱼ, r₀)
        R     = fixedJacPoint(f, Jⱼ, r, τ)
        Jⱼ    = JacT(tⱼ + τ, tⱼ, zⱼ, Jⱼ, r, R)
        # compute inner approximation
        iaⱼ   = innerApproximate(z₀ⱼ, Jⱼ, z₀)
        # save the outer approx. zⱼ and inner approx. iaⱼ
        tⱼ += τ
        st    = vcat(st,  tⱼ)
        sz    = hcat(sz,  Interval.(zⱼ))
        sia   = hcat(sia, iaⱼ)
        iter += 1
    end
    
    return st, sz, sia
end
