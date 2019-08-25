using AffineArithmetic
using ForwardDiff
using Profile

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


f(z::Vector) = [1.0 - 2.5*z[1] + z[2]*z[1]*z[1]; 1.5*z[1] - z[2]*z[1]*z[1]]
z₀ = [0.9, 0.1]

tm = constructTM(f; order=4, T=Real)
Profile.init(n = 10^8, delay = 0.0001)
@profile a = tm(0.1, 1, z₀, z₀)
#a = tm(0.1, 1, z₀, z₀)
Profile.print()













