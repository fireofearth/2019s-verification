
## Original implementation

Good news! We have the software that was used for Goubault+Putot's last 3 years of research: Robust INner and Outer Approximated Reachability (RINO) <https://github.com/cosynus-lix/RINO> as well as the dependencies FILIB++ <http://www2.math.uni-wuppertal.de/wrswt/software/filib.html> for interval computations, and FADBAD <http://www.fadbad.com/fadbad.html> for automatic diff.

### Details and specifications

#### Unused for loop

The below for loop and function `init_subdiv()` is never used so we can remove it entirely. To simplify we may ignore the variables `inputs_save`, `center_inputs`.

```
for (current_subdiv=1; current_subdiv <= nb_subdiv_init; current_subdiv++) {
    if (nb_subdiv_init > 1) init_subdiv(current_subdiv, inputs_save, 0);
    // ...
}
```

#### Redundant variables

In all simulations, `jacdim` and `sysdim` contain identical values. I am choosing `sysdim`.

#### Modal Arithmetic

Modal arithmetic is used at the last step of the procedure run at every time step via the `InnerOuter()` function to obtain inner approximations from the taylor models. All other calculations are done using taylor models using affine arithmetic.

#### Affine Arithmetic

Affine arithmetic is done using the `aaflib` library

#### Automatic Differentiation

Automatic differentiation is called on the jacobians of the ODE function at each time step.

### Classes

#### OdeVar

Is used to represent the taylor coefficients of the taylor model in `TM_Jac`. Automatic differentiation using FADBAD is called on the member variables of `OdeVar`.

x (vector<T<AAF>>): independent variables  
xp (vector<T<AAF>>): dependent variables  

#### Ode

Is used to represent the taylor coefficients of the taylor model in `TM_val`.

`Ode` members:

x (vector<T<AAF>>): independent variables  
xp (vector<T<AAF>>): dependent variables  

`Ode` functions:

`Ode::Ode()` constructor  
`Ode::Ode(OdeFunc f)` constructor that sets `x`, `xp` by evaluating them using the function `f` passed to it.  

#### TM_val

Is a vector Taylor Model

TODO: is it necessary that all points are AAF?

`TM_val` members:

ode_x (Ode): taylor coefficients from 0 to k-1 of TM  
ode_g (Ode): the k-th taylor coefficient of TM  

`TM_val` functions:

`TM_val::build()`  


`TM_val::eval()`  

#### HybridStep_ode

`HybridStep_ode` members:

bf (OdeFunc):  
TMcenter (TM_val):  
TMJac (TM_Jac):  
order (int): order of the taylor model  
tn (double): the time at the n-th iteration  
tau (double): the time step  

`HybridStep_ode` functions:

`HybridStep_ode::HybridStep_ode()` constructor  
`HybridStep_ode::TM_build()` build Taylor Model using ODE function `bf` by calling `TMcenter.build()` and `TMJack.build()`.  
`HybridStep_ode::TM_eval()` calls `TMcenter.eval()` and `TMJac.eval()`.  
`HybridStep_ode::init_nextstep()`  

### Steps

Step 4: create an ODE object `HybridStep_ode` using `HybridStep_ode::init_ode()`. 

Step 5.1: build Taylor Model using `HybridStep_ode::TM_build()`. 

## Julia RINO implementation

Required components to reproduce RINO for Julia:

- automatic diff <https://github.com/JuliaDiff>
    + I'm concerned this library is unusable for functions that have inputs and outputs that are not type `Real`. 

- intervals <https://github.com/JuliaIntervals/IntervalArithmetic.jl>
    + I'm implementing modal intervals on top of IntervalArithmetic library.
    + I'm considering the reimplementation of IntervalArithmetic itself because the rounding of improper integrals uses inner rounding.

- affine arithmetic <https://github.com/JuliaIntervals/AffineArithmetic.jl>
    + This library is a thin implementation so I am implementing affine arithmetic from the beginning.

- taylor models <https://github.com/JuliaIntervals/TaylorModels.jl>


### Components

ModalInterval  
AffineArithmetic  
AffineTaylorModel  
ODECommon  
ODEIntegration  

#### Automatic differentiation

### Procedure

I am only implementing RINO for ODE functions, and am skipping DDE for now.

TODO: figure out implementation

Step 2: preallocates the J, x, xCenter

Step 3: loop in solveODE

Step 3.2: obtain the center of the inputs, and eps

Step 3.3: set initial values for J, x, xCenter

Step 3.4 set initial ODE system
