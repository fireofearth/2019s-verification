
include("ModalInterval.jl")
include("Affine.jl")

@enum Problem NN BRUSSELATOR VANDERPOLL CIRCLE BOX LORENZATTR BALLISTIC LOTKAVOLTERRA FITHUGHNAGUMO

 #=
 # Initial variables of ODE problem
 #
 # TODO: it is not clear whether jacdim = sysdim in RINO
=#
struct ODEInitVar
    sysdim::Unsigned
    jacdim::Unsigned
    numSubdiv::Unsigned
    tBegin::AbstractFloat
    tEnd::AbstractFloat
    tau::AbstractFloat
    order::Unsigned

    inputs::NTuple{<:AbstractInterval}
    isControlled::NTuple{<:Bool}
    isInitialConditoin::NTuple{<:Bool}
    isVariable::NTuple{<:Bool}
    numControlled::Unsigned
    numUncontrolled::Unsigned

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

 #=
 #
 #
 #
=#
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
 # TODO: do we want to use tuple?
=#
#function (o::ODEFunc)(x::Tuple{AbstractFloat,Vararg{<:Unsigned}})::Tuple{AbstractFloat,Vararg{<:Unsigned}}
function (o::ODEFunc)(x::AbstractArray)::AbstractArray
    if(o.p == NN)
        error("ODEFunc contains no problem")
    elseif(o.p == BRUSSELATOR)
        return SVector(1. - (2.5 * x[1]) + (x[1] * x[1] * x[2]),
                       (1.5 * x[1]) - (x[1] * x[1] * x[2]))
    else
        error("ODEFunc no problem specified")
    end
end
