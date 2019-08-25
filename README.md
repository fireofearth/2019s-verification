# Flowpipes.jl

In the conference HSCC'17-19, Goubault and Putot presented a series of papers listed in the [related works section](#related-works) and an accompanying demonstration software called [RINO](#related-software).
This project uses RINO to create a numerical solver that computes inner and outer approximations of flowpipes for ODEs. Currently it only implements the algorithm present in HSCC'17.

## Main Dependencies

A modified version of JuliaDiff/[ForwardDiff.jl](https://github.com/fireofearth/ForwardDiff.jl)  
[AffineArithmetic.jl](https://github.com/fireofearth/AffineArithmetic.jl)  
[ModalIntervalArithmetic.jl](https://github.com/fireofearth/ModalIntervalArithmetic.jl)  
[IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl)  
[TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl) used in some example code  

## Installation

1. Assuming Julia is already installed, install `IntervalArithmetic` and `TaylorSeries` from official repositories. 
2. Install ForwardDiff dependencies `DiffResults`, `DiffRules`, `StaticArrays`, `SpecialFunctions`, `NaNMath` , `CommonSubexpressions` from official repositories.
3. Clone `AffineArithmetic.jl`, `ModalIntervalArithmetic.jl` and `ForwardDiff.jl` from the Github repositories.
4. So the Julia REPL knows where to find the modules you've cloned, add the following to `~/.julia/config/startup.jl` (Julia startup file):

```julia
verificationModulePaths = [
    "/path/to/AffineArithmetic.jl/src"
    "/path/to/ModalIntervalArithmetic.jl/src"
    "/path/to/ForwardDiff.jl/src"
]

for path in verificationModulePaths
    if(!(path in LOAD_PATH))
        push!(LOAD_PATH, path)
    end
end
```

5. Clone `Flowpipes.jl` from this repository. To test that it works do `include("test/runtests.jl")` in the REPL.

## Usage

Flowpipes comes with a single function `approximate(f::Function, tspan::NTuple{2,<:Real}, τ::Real, z₀::Vector{<:Interval}; order::Int=4)`

Arguments:

- `f` is the function representing the system of autonomous first order ODEs `\dot{x} = f(x)`  
- `tspan=(tstart, tend)` contains the start and end time the solution of the system is defined on. We build the inner and outer approximation of flowpipe over this time period.  
- `τ` is the time step  
- `z₀` is an interval representing the bounds of the initial values that forms the flowpipe.  
- `order` (optional) is the order of taylor approximation we use to build the outer approximation at each step.

Returns:

 - `st` vector of evenly spaced points of time from `t₀` to `tn`
 - `sz` vector of intervals representing outer approximations at each time point. sz[i] is the outer approximation at st[i]
 - `sia` vector of intervals representing of inner approximationss, or NaN if does not exist at each time point. sia[i] is the inner approimxation/NaN at st[i]

To use simply add `include("~/path/to/Flowpipes.jl")` in your file.
For a runnable example take a look at `./examples/approx_2.jl`.

## TODOs

- ForwardDiff takes up much CPU cycles and has a long runtime. See `./benchmarks/profiler.txt` for Profiler results. Most of the cycles are spent in `jacobian.jl` and `partials.jl`.
This is most likely because ForwardDiff does not use source code transformation, rather it computes derivatives (and higher-order derivatives) with operator overloading using dual numbers. This means for each higher-order derivative, hyper-dual numbers must be instantiated and calculated instead of one pass with pre-constructed functions `f'`, `f''`, `f'''`, etc.

- ForwardDiff does not natively support obtaining derivatives `f'(x)` of functions evaluated with types other than real numbers (i.e. Interval and Affine). Instead ForwardDiff must be altered by replacing instances of `Real` with `ResolvableType = Union{Affine, Taylor1, Real}`.
This may mean that 

- Flowpipes gives woefully poor inner/outer approximations compared to RINO. I'm not sure exactly why and it's something to investigate.

## Related Works

Eric Goubault and Sylvie Putot. 2017. Forward inner-approximated reachability of non-linear continuous systems. In *Proceedings of the 20th. ACM International Conference on Hybrid Systems: Computation and Control (HSCC '17)*. ACM Press, New York, NY, 1-10. DOI:<https://doi.org/10.1145/3049797.3049811> PDF:<http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/hscc17.pdf>

Alexandre Goldsztejn1, David Daney1, Michel Rueher1, and Patrick Taillibert. 2005. Modal Intervals Revisited: a mean-value extension to generalized intervals. In *Proceedings of the First International Workshop on Quantification in Constraint Programming (QCP '05)*. ??? PDF:<http://goldsztejn.com/publications/QCP2005.Goldsztejn-Daney-Rueher-Taillibert.pdf>

Miguel Sainz et al. 2014. *Modal Interval Analysis: New Tools for Numerical Information*. Springer Lecture Notes in Mathematics, New York, NY.
PDF:<https://www.springer.com/gp/book/9783319017204>

Jorge Stolfi and Luiz Henrique de Figueiredo. 1997. Self-Validated Numberical Methods and Applications. In *Proceedings of the 21st Brazilian Mathematics Colloquium*. ??? PDF:<http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.36.8089&rep=rep1&type=pdf>

## Related Software

Robust INner and Outer Approximated Reachability (RINO)
<https://github.com/cosynus-lix/RINO>

Affine Arithmetic C++ Library (aaflib) 
<http://aaflib.sourceforge.net>

TaylorModels: Rigorous function approximation using Taylor models in Julia. This package produces tight outer approximations to flowpipes similar to Flowpipes.
<https://github.com/JuliaIntervals/TaylorModels.jl >

