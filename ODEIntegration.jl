
using ForwardDiff
using AffineArithmetic

 #=
 # Construct Taylor Model
 #
 # Specification:
 # - Given problem ̇z' = f(z), we construct the function
 # [z](t, tⱼ, [zⱼ]) = [zⱼ] = ∑{i=1,…,k-1} (t-tⱼ)ⁱ/i! f⁽ⁱ⁾([zⱼ]) + (t-tⱼ)ᵏ/k! f⁽ᵏ⁾([rⱼ₊₁])
 # used by the HSCC'17 article by Goubault+Putot and returns it
 #
 # - Assumes that f: Rᴺ → Rᴺ
 #
 # - We do the calculation with affine forms as initial conditions, this implies 
 # [zₒ] ≡ Affine(center([zₒ])), etc... However we permit evaluating for any type T for zⱼ
 #
 # TODO:
 # - possible overflow when calculating fartorials, orders
=#
function constructTMAffine(f::Function, order::Int)
    # store the lie derivatives f⁽ⁱ⁾ (i = 1,…,order) in vf
    vf = [f]
    for i in 2:order
        fi = (z::Vector{Affine} -> ForwardDiff.jacobian(vf[i - 1], z) * f(z))
        vf = vcat(vf, fi)
    end

    # constructed TM
    # cTM(t,tⱼ,[zⱼ],[rⱼ₊₁]) = [zⱼ] = ∑ᵢ(t-tⱼ)ⁱ/i! f⁽ⁱ⁾([zⱼ]) + (t-tⱼ)ᵏ/k! f⁽ᵏ⁾([rⱼ₊₁])
    function cTM(t::Real, tj::Real, zj::Vector{Affine}, r::Vector{Affine})
        acc = zj
        for i in 1:(order - 1)
            acc += ((t - tj)^i / factorial(i)) * vf[i](zj)
        end
        acc += ((t - tj)^order / factorial(order)) * vf[order](r)
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
    # cTM(t,tⱼ,[zⱼ],[rⱼ₊₁]) = [zⱼ] = ∑ᵢ(t-tⱼ)ⁱ/i! f⁽ⁱ⁾([zⱼ]) + (t-tⱼ)ᵏ/k! f⁽ᵏ⁾([rⱼ₊₁])
    function cTM(t::Real, tj::Real, zj::Vector{<:Real}, r::Vector{<:Real})
        acc = zj
        for i in 1:(order - 1)
            acc += ((t - tj)^i / factorial(i)) * vf[i](zj)
        end
        acc += ((t - tj)^order / factorial(order)) * vf[order](r)
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
=#
function constructJacTMAffine(f::Function; order::Int=5)

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
        vJacf = vJacf(vJacf, Jacfi)
    end

    # cJacTM(t,tⱼ,[zⱼ],[Jⱼ],[rⱼ₊₁],[Rⱼ₊₁])
    # = [Jⱼ] = ∑{i=1,…,k-1} (t-tⱼ)ⁱ/i! Jac f⁽ⁱ⁾([zⱼ]) [Jⱼ] + (t-tⱼ)ᵏ/k! Jac f⁽ᵏ⁾([rⱼ₊₁]) [Rⱼ₊₁]
    function cJacTM(t::Real, tj::Real, zj::Vector{Affine}, 
                    Jj::Matrix{Affine}, r::Vector{Affine}, R::Matrix{Affine})
        acc = Jj
        for i in 1:(order - 1)
            acc += ((t - tj)^i / factorial(i)) * vJacf[i](zj) * Jj
        end
        acc += ((t - tj)^order / factorial(order)) * vf[order](r) * R
        return acc
    end

    return cJacTM
end

function constructJacTMReal(f::Function; order::Int=5)

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
        vJacf = vJacf(vJacf, Jacfi)
    end

    # cJacTM(t,tⱼ,[zⱼ],[Jⱼ],[rⱼ₊₁],[Rⱼ₊₁])
    # = [Jⱼ] = ∑{i=1,…,k-1} (t-tⱼ)ⁱ/i! Jac f⁽ⁱ⁾([zⱼ]) [Jⱼ] + (t-tⱼ)ᵏ/k! Jac f⁽ᵏ⁾([rⱼ₊₁]) [Rⱼ₊₁]
    function cJacTM(t::Real, tj::Real, zj::Vector{<:Real}, 
                    Jj::Matrix{<:Real}, r::Vector{<:Real}, R::Matrix{<:Real})
        acc = Jj
        for i in 1:(order - 1)
            acc += ((t - tj)^i / factorial(i)) * vJacf[i](zj) * Jj
        end
        acc += ((t - tj)^order / factorial(order)) * vf[order](r) * R
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
 # The Picard Lindelof is expressed as F[z(t)] = z₀ + ∫{tⱼ to tⱼ₊₁} f(z(s)) ds 
 # With the interval extension F[z] = z₀ + (tⱼ₊₁ - tⱼ)[f]([z]) with [f] being the natural
 # interval extension of f. Both versions of F admits a unique fixed point if f is Lipschitz.
 #
 # TODO: apply sumup to affine zi so we don't have too accumulate coefficients
=#
function fixedPoint(f::Function, z0::Vector{Affine}, n::Int, τ::Real)
    iiter = 1
    zi = z0
    T = fill(Affine(interval(0, τ)), n)
    X = fill(Affine(interval(-1,  1)),   n)

    while true
        Fzi = z0 + (T .* f(zi))
        i > 1 && Interval(Fzi) ⊆ Interval(zi) && break
        zi = Fzi
        
        # reasoning of these choices?
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
        elseif(iter > 2)
            β = 0.00001
        end
        if(iter > 2) # what is this guard for?
            zi = zi + β*(X .* zi)
        end
        iter += 1
    end

    return zi
end

 #=
 # Entrypoint to solver
 #
 # Specifications:
 # - initializes variables
 #
 # TODO: finish procedure
=#
#=
function solveODE(odef::ODEFunc)

     #=
     # Initializes arrays
     # - J::Matrix{AAF} the Jacobian
     # - x::Vector{AAF}
     # - xCenter::Vector{AAF}
     # - eps::Vector{Interval}
     #
     # Merging to set_initialconditions() second half of init_system() in RINO,
     # and discarding init_subdiv()
    =#
    J = fill(AAF(0.0), odev.sysdim, odev.sysdim)
    for ii in 1:jacdim
        J[ii, ii] = AAF(1.0)
    end
    x = [AAF(odev.inputs[ii]) for ii in 1:odev.sysdim]
    xCenter = [AAF(getCenter(odev.inputs[ii])) for ii in 1:odev.sysdim]
#    eps = [
#        Interval(odev.inputs[ii] - getCenter(odev.inputs[ii]))
#            for ii in 1:odev.sysdim
#    ]

     #=
     # Print output, generate plots, etc
    =#
end
=#
