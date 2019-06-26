using IntervalArithmetic, StaticArrays

include("interval.jl")
include("affine.jl")

@enum Problem NN BRUSSELATOR VANDERPOLL CIRCLE BOX LORENZATTR BALLISTIC LOTKAVOLTERRA FITHUGHNAGUMO

 #=
 # Initial variables of ODE problem
 #
 # TODO: it is not clear whether jacdim = sysdim in RINO
=#
struct ODEInitVar
    const sysdim::Unsigned
    const jacdim::Unsigned
    const numSubdiv::Unsigned
    const t_begin::AbstractFloat
    const t_end::AbstractFloat
    const tau::AbstractFloat
    const order::Unsigned

    const inputs::NTuple{<:AbstractInterval}
    const isControlled::NTuple{<:Bool}
    const isInitialConditoin::NTuple{<:Bool}
    const isVariable::NTuple{<:Bool}
    const numControlled::Unsigned
    const numUncontrolled::Unsigned

     #=
     # Initialize variables where inputs types are not specified
    =#
    function ODEInitVar(sysdim::Unsigned, jacdim::Unsigned, nb_subdiv_init::Unsigned, t_begin::AbstractFloat, t_end::AbstractFloat, tau::AbstractFloat, order::Unsigned, inputs::NTuple{<:AbstractInterval})
        numInputs = length(inputs)
        new(sysdim, jacdim, nb_subdiv_init, t_begin, t_end, tau, order, inputs,
            Tuple(false for i=1:numInputs),
            Tuple(false for i=1:numInputs),
            Tuple(false for i=1:numInputs),
            numInputs, numInputs)
    end
end

struct ODEFunc
    p::Problem

    function ODEFunc(p::Problem)
        new(p)
    end
end

 #=
 #
 #
 # TODO: finish ODEFunc, add more problems
=#
function getInitialVar(o::ODEFunc)::ODEInitVar
    if(o.p == NN)
        error("ODEFunc contains no problem")
    elseif(o.p == BRUSSELATOR)
        return ODEInitVar( 2, 2, 1, 0.05, 0, 10., 4,
                      ( Interval(0.9, 1), Interval(0, 0.1) ))
    else
        error("ODEFunc no problem specified")
    end
end

 #=
 #
 #
 # TODO: finish ODEFunc, add more problems
=#
function (o::ODEFunc)(x::Tuple{AbstractFloat,Vararg{<:Unsigned}})::Tuple{AbstractFloat,Vararg{<:Unsigned}}
    if(o.p == NN)
        error("ODEFunc contains no problem")
    elseif(o.p == BRUSSELATOR)
        return (1. - (2.5 * x[1]) + (x[1] * x[1] * x[2]),
                (1.5 * x[1]) - (x[1] * x[1] * x[2]))
    else
        error("ODEFunc no problem specified")
    end
end

struct IOR
    J::Array{Union{Missing,AAF}}
    x::Array{Union{Missing,AAF}}
    xCenter::Array{Union{Missing,AAF}}
    
    function IOR(sysdim::Unsigned, jacdim::Unsigned)
Vector{Union{Missing,AAF}}(missing, sysdim)
        new()
end

 #=
 #
 #
 # TODO: finish function
=#
function initSubdiv(odev::ODEInitVar)
end

 #=
 #
 #
 # TODO: finish procedure
=#
function solveODE(odef::ODEFunc)
    odev = getInitialVar(odef)

     #=
     # Initialize x, J
    =#

     #=
     # Solver loop for subdivisions
     #
     # TODO: what is a subdivision exactly?
    =#
    for currentSubdiv in 1:numSubdiv
        if(numSubdiv > 1) initSubdiv(odev)
    end

     #=
     # Print output, generate plots, etc
    =#
end

# example using Brusselator
f = ODEFunc(BRUSSELATOR)
print(f((0.5, 0.5)))
