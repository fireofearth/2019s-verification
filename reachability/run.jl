using IntervalArithmetic

include("interval.jl")

# In Modal Intervals, proper ≡ ∃, improper ≡ ∀, and real are both ∀,∃
#=
function oldsubseteq(A::ModalInterval, B::ModalInterval)
    if(mod(A) == mod(B) == proper)
        return prop(A) ⊆ prop(B)
    elseif(mod(A) == mod(B) == improper)
        return prop(A) ⊇ prop(B)
    elseif(mod(A) == improper && mod(B) == proper)
        return prop(A) ∩ prop(B) = ∅
    else
        return prop(A) == prop(B) && inf(A) == sup(A)
    end
end
=#

X = interval(1,2)
print("Inf(X) = $(X.lo)\n")
print("Sup(X) = $(X.hi)\n")
print("mid(X) = $(mid(X))\n")
print("mid(X,α=0.7) = $(mid(X,0.7))\n")

X = ModalInterval(1,2) # proper
print("X' = [$(X.prime.lo), $(X.prime.hi)]\n")
print("X = [$(inf(X)), $(sup(X))]\n")

Y = ModalInterval(3,0) # improper
print("Y = [$(inf(Y)), $(sup(Y))]\n")
display(X + Y)
display(X - Y)


