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
 # TODO: check J dimensions (row, col of J(sysdim, vector<AAF>(jacdim))?)
 # TODO: do I want to initialize arrays
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
 #
 #
 # TODO: finish procedure
=#
function solveODE(odef::ODEFunc)

     #=
     # Similar to initSystem in RINO
    =#
    odev = getInitialVar(odef)

     #=
     # Similar to initialize J, x, x_center in RINO
    =#
    ior = IOR(sysdim, jacdim)

     #=
     # Solver loop for subdivisions
     #
     # TODO: what is a subdivision exactly?
    =#
    for currentSubdiv in 1:numSubdiv
        if(numSubdiv > 1)
            initSubdiv(odev, ior)
        end
        setInitialConditions(ior)
        ior.t = odev.tBegin
        # RINO prints stats here
        step = HybridStepODE()
    end

     #=
     # Print output, generate plots, etc
    =#
end
