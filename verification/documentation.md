
## Original implementation

Step 3: `inputs` are saved to `inputs_save`. The variable `inputs` is reassigned to be subdivisions of the original inputs now stored as `inputs_save`.
`init_subdiv()` sets up inputs

```
center_inputs = vector<AAF>(jacdim);
eps = vector<interval>(jacdim);
    
interval save = inputs_save[param_to_subdivide].convert_int();
double delta = (save.sup()-save.inf())/nb_subdiv_init;

if ((current_subdiv > 1) && (current_subdiv < nb_subdiv_init))
    inputs[param_to_subdivide] = interval(save.inf()+delta*(current_subdiv-1-recovering),save.inf()+delta*(current_subdiv+recovering));
else if (current_subdiv == 1)
    inputs[param_to_subdivide] = interval(save.inf()+delta*(current_subdiv-1),save.inf()+delta*(current_subdiv+recovering));
else if (current_subdiv == nb_subdiv_init)
    inputs[param_to_subdivide] = interval(save.inf()+delta*(current_subdiv-1-recovering),save.inf()+delta*(current_subdiv));

cout << "inputs[param_to_subdivide] " << inputs[param_to_subdivide] << endl;
    
interval temp = inputs[param_to_subdivide].convert_int();
center_inputs[param_to_subdivide] = mid(temp);
eps[param_to_subdivide] = temp-mid(temp);
```

## Julia RINO implementation

Step 2: preallocates the J, x, xCenter

Step 3: loop in solveODE

Step 3.2: obtain the center of the inputs, and eps

Step 3.3: set initial values for J, x, xCenter

Step 3.4 set initial ODE system
