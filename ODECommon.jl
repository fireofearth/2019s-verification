
@enum Problem NN BRUSSELATOR VANDERPOLL CIRCLE BOX LORENZATTR BALLISTIC LOTKAVOLTERRA FITHUGHNAGUMO

 #=
 # Initial variables of ODE problem
=#
struct ODEInitVar
    sysdim::Unsigned
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
    function ODEInitVar(sysdim::Unsigned, nb_subdiv_init::Unsigned, t_begin::AbstractFloat, t_end::AbstractFloat, tau::AbstractFloat, order::Unsigned, inputs::NTuple{<:AbstractInterval})
        numInputs = length(inputs)
        new(sysdim, jacdim, nb_subdiv_init, t_begin, t_end, tau, order, inputs,
            Tuple(false for i=1:numInputs),
            Tuple(false for i=1:numInputs),
            Tuple(false for i=1:numInputs),
            numInputs, numInputs)
    end
end

 #=
 # Struct specifies which ODE we choose
=#
struct ODEFunc
    p::Problem

    function ODEFunc(p::Problem)
        new(p)
    end
end

 #=
 # (Initial) values for problem instance
 #
 #
 # 
 #
 # TODO: finish ODEFunc, add more problems
=#
function getInitialVar(o::ODEFunc)::ODEInitVar
     #=
     # Assigns values to the following:
     #
     # ODEInitVar.sysdim::Unsigned
     # ODEInitVar.jacdim::Unsigned
     # ODEInitVar.nb_subdiv_init::Unsigned
     # ODEInitVar.t_begin::AbstractFloat
     # ODEInitVar.t_end::AbstractFloat
     # ODEInitVar.tau::AbstractFloat
     # ODEInitVar.order::Unsigned
     # ODEInitVar.inputs::NTuple{<:AbstractInterval})
     #
    =#
    if(o.p == NN)
        error("ODEFunc contains no problem")
    elseif(o.p == BRUSSELATOR)
        return ODEInitVar( 2, 1, 0.05, 0, 10., 4,
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
