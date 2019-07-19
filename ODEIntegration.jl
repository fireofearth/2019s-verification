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
 # TODO: 
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
     #
     # Similar to first half of init_system in RINO.
    =#
    odev = getInitialVar(odef)

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
    eps = [
        Interval(odev.inputs[ii] - getCenter(odev.inputs[ii]))
            for ii in 1:odev.sysdim
    ]

    # RINO prints stats here
    
    # 
     #=
     #
     # init_ode()
     #
     # TODO: passing in odef, xcenter, x, J, tn, tau, order
     # TODO: set up these variables:
     # odeVAR_x, odeVAR_g of form 
    =#
    step = HybridStepODE()

    # TM_val

         #=
         # Loop iteration procedure.
         # For each step from tBegin to tEnd do
         # (2) TM build
         # (3) TM evaluate (and print intermediates)
         # (3.1)
        =#
        #while

     #=
     # Print output, generate plots, etc
    =#
end
