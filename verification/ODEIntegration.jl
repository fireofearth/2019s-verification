include("ODECommon.jl")

 #=
 # params:
 # J - Jacobian of ODE function
 # x - vector of affine forms representing ODE parameters
 # t - the current time
 #
 # Specification:
 # - Based on RINO by Goubault/Putot
 # - Use Immutable (functional) types. No global variables or in place object member mutations.
 # - 
 #
 # TODO: finish IOR
 # TODO: check J dimensions (row, col of J(sysdim, vector<AAF>(jacdim))?); RINO uses which onesfor cols, rows?
 # TODO: is it the case sysdim = jacdim? Redundent variable?
=#
struct IOR
    J::Matrix{<:Union{Missing,AAF}}
    x::Vector{Union{Missing,AAF}}
    xCenter::Vector{Union{Missing,AAF}}
    centerInputs::Vector{Union{Missing,AAF}}
    eps::Vector{Union{Missing,AAF}}
    t::AbstractFloat
    
    function IOR(sysdim::Unsigned, jacdim::Unsigned)
        new(Matrix{T}(missing, jacdim, sysdim),
            Vector{Union{Missing,AAF}}(missing, sysdim),
            Vector{Union{Missing,AAF}}(missing, sysdim),
            Vector{Union{Missing,AAF}}(missing, sysdim),
            Vector{Union{Missing,AAF}}(missing, sysdim),
            0)
    end
end

 #=
 #
 #
 # TODO: finish function
=#
function initSubdiv(odev::ODEInitVar, ior::IOR)
end

struct HybridStepODE
    # ODEFunc is declared here in RINO
    
    function HybridStepODE()
    end
end

 #=
 # Entrypoint to solver
 #
 # Specifications:
 # - initializes variables
 #
 # TODO: finish procedure
=#
function solveODE(odef::ODEFunc)

     #=
     # Get initial values of ODE problem
     # Specification: members in odev are never mutated
     # Similar to first half of init_system in RINO.
    =#
    odev = getInitialVar(odef)

     #=
     # Preallocates the array 
     # - J
     # - x
     # - xCenter
     #
     # Similar to initialize J, x, x_center in RINO
    =#
    J = fill(AAF(0.0), odev.sysdim, odev.sysdim)
    x = fill(AAF(0.0), odev.sysdim)
    xCenter = fill(AAF(0.0), odev.sysdim)

     #=
     # Initializes array
     # - centerInputs::Vector{AAF}
     # - eps::Vector{Interval}
     #
     # Similar to second half of init_system in RINO.
    =#
    centerInputs = [
        AAF(getCenter(odev.inputs[ii])) for ii in 1:odev.sysdim
    ]
    eps = [
        Interval(odev.inputs[ii] - getCenter(odev.inputs[ii]))
            for ii in 1:odev.sysdim
    ]

     #=
     # Solver loop for subdivisions
    =#
    for currentSubdiv in 1:numSubdiv

        inputs = odev.inputs
         #=
         # Similar to init_subdiv()
         # TODO: improve this code
         # TODO: what is a subdivision exactly?
         # TODO: why is param_to_subdivide = 0 always in RINO?
        =#
        if(odev.numSubdiv > 1)
            delta = ( sup(odev.inputs[1]) - inf(odev.inputs[1]) ) / odev.numSubdiv
            inputs[1] = AAF(Interval(
                 getMin(inputs[1]) + delta*(currentSubdivid - 1),
                 getMin(inputs[1]) + delta*currentSubdivid
            ))
            centerInputs[1] = getCenter(inputs[1])
            eps[1] = Interval(inputs[1] - getCenter(inputs[1]))
        end

         #=
         # Initialize J = Id containing AAF, x = inputs, xCenter to centerInputs
         # Similar to set_initialconditions() in RINO.
        =#
        for ii in 1:jacdim
            J[ii, ii] = AAF(1.0)
            x[ii] = inputs[ii]
            xCenter
        end

        ior.t = odev.tBegin
        # RINO prints stats here
        step = HybridStepODE()

         #=
         # Loop iteration procedure.
         # For each step from tBegin to tEnd do
         # (2) TM build
         # (3) TM evaluate (and print intermediates)
         # (3.1)
        =#
        while()
    end

     #=
     # Print output, generate plots, etc
    =#
end
