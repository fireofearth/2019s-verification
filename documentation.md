
## Verification Study

Two remote repositories:

```
git remote set-url origin https://github.com/fireofearth/2019s-verification
git remote set-url origin https://bitbucket.org/ian_mitchell/julia-intervals
git push origin master
```

*20190617*

We have the software that was used for Goubault+Putot's last 3 years of research: Robust INner and Outer Approximated Reachability (RINO) [here](https://github.com/cosynus-lix/RINO) as well as the dependencies[FILIB++](http://www2.math.uni-wuppertal.de/wrswt/software/filib.html) for interval computations, and [FADBAD](http://www.fadbad.com/fadbad.html > for automatic diff)

Additionally for Julia (all required components to reproduce RINO for Julia):
- automatic diff <https://github.com/JuliaDiff/>
- intervals <https://github.com/JuliaIntervals/IntervalArithmetic.jl>
- affine arithmetic <https://github.com/JuliaIntervals/AffineArithmetic.jl>
- taylor models <https://github.com/JuliaIntervals/TaylorModels.jl>

It was used in the latest paper for HSCC 2019 "Inner and Outer Reachability for the Verification of Control Systems"; and the paper we're reading that is HSCC 2017 "Forward inner-approximated reachability of non-linear continuous systems"; and the paper CAV 2018 "Inner and Outer Approximating Flowpipes for Delay Differential Equations".

*20190626*

Affine arithmetic does not seem to be supported (well) in Julia. Interval arithmetic does not have modal interval extensions in Julia. Ian Mitchell (IM) asked me to check with JuliaIntervals to ask whether developers (Dr David P. Sanders, Chris Rackauckas) would be open to contributions from for libraries.

*20190702*

Step 2: for each subdivision

call set_initialconditions (ode_def.cpp) to set x to inputs which are intervals, xCenter which are midpoints of intervals (as floats), and J = Id. Recall x, xCenter and J are all vectors and matrices containing AAF.

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

### Testing

I'm writing a test suite in RINO to test their modification of the aaflib affine arithmetic library

```
g++ -ggdb -frounding-math -DMAXORDER=40 -I. -I${HOME}/lib/filib-3.0.2/include \
    -I/usr/include -I$(pwd)/aaflib-0.1 -fpermissive -std=c++11 -c testsuite.cpp

g++ -L/usr/lib -L$(pwd)/aaflib-0.1 -L${HOME}/lib/filib-3.0.2/lib -o testsuite \
    testsuite.o -laaf -lprim -lgsl -llapack -lblas -lcblas -lstdc++ \
    -lboost_unit_test_framework

./testsuite --log_level=test_suite

# in the case that the testsuite does not detect libaaf.so we need to add the library path for this shared resource.
LD_LIBRARY_PATH=/usr/lib:/home/fireofearth/res/mitchell-ian/programs/rino-hscc17/aaflib-0.1/
export LD_LIBRARY_PATH
```

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

[Cassette](https://jrevels.github.io/Cassette.jl/latest/whycassette.html) was originally designed for better language support for automatic differentiation.
[AutoGrad](https://github.com/denizyuret/AutoGrad.jl) is a port of Python [autograd](https://github.com/HIPS/autograd). I'm avoiding this package due to lack of documentation

Both [AutoDiffSource](https://github.com/gaika/AutoDiffSource.jl) and ReverseDiffSource are deprecated and limited to Julia version 0.5 (and hence in JuliaAttic)

By the process of elimination [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) is the one I will be using. I have to modify the source code so that we support Affine types. In addition Affine must support the following:

`one` - return a multiplicative identity for x: a value such that one(x)\*x == x\*one(x) == x. Alternatively one(T) can take a type T, in which case one returns a multiplicative identity for any x of type T.

`isone` - return true if x == one(x) if x is an array, this checks whether x is an identity matrix.

### Procedure

I am only implementing RINO for ODE functions, and am skipping DDE for now.

TODO: figure out implementation

Step 2: preallocates the J, x, xCenter

Step 3: loop in solveODE

Step 3.2: obtain the center of the inputs, and eps
    
Step 3.3: set initial values for J, x, xCenter

Step 3.4 set initial ODE system
