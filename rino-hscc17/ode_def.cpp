
/* ============================================================================
 File   : ode_def.cpp
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
The place to define initial conditions and parameters the systems of ODEs or DDEs on which to perform reachability
============================================================================ */

#include <assert.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>

#include "filib_interval.h"
#include "tadiff.h" 
#include "fadiff.h"
#include "fadbad_aa.h"
#include "ode_def.h"
#include "matrix.h"

using namespace std;

int sysdim; // dimension of system of ODE/DDE
int jacdim;  //  Jacobian will be dimension sysdim * jacdim
int sysdim_params;
//int sysdim_jacparams;

// parameters of the system of the ODEs
vector<AAF> params;  // params of the ODE (nondeterministic disturbances)
vector<AAF> inputs; // uncertain inputs and parameters : some will be used in initial condition, some as uncertain parameters
vector<AAF> center_inputs;
vector<interval> eps;

// for subdivisions of the initial domain to refine precision
int nb_subdiv_init; // number of subdivisiions
double recovering; // percentage of recovering between subdivisions
vector<vector<vector<interval>>> Xouter_print, Xouter_robust_print, Xouter_minimal_print, Xinner_print, Xinner_robust_print, Xinner_minimal_print, Xexact_print; // store results of subdivision
vector<double> t_print; // times where results are stored
int current_subdiv;
int current_iteration;

// for robust inner-approximations
int uncontrolled; // number of uncontrolled parameters (forall params)
int controlled; // number of controlled parameters (forall params)
vector<bool> is_uncontrolled; // for each input, uncontrolled or controlled (for robust inner-approx)
vector<bool> is_initialcondition; // for each input, initial condition or parameter (for robust outer-approx)
int variable; // number of non constant parameters
vector<bool> is_variable;  // for each parameter, constant or variable

// for ODEs : initialize the state variable (and center for inner-approximation)
void set_initialconditions(vector<AAF> &x, vector<AAF> &xcenter, vector<vector<AAF>> &J)
{
    // par défaut
    for (int i=0 ; i<sysdim ; i++) {
        x[i] = inputs[i];
        xcenter[i] = center_inputs[i];
    }
    setId(J);
}

// for ODEs and DDEs: define bounds for parameters and inputs, value of delay d0 if any, and parameters of integration (timestep, order of TM)
void init_system(double &t_begin, double &t_end, double &tau, double &d0, int &nb_subdiv, int &order /*, vector<interval> &ix*/) {

    /**
     * Colin: problem dimensions
     */
    sysdim = 2;
    jacdim = 2;
    sysdim_params = 2;
    // nb of initial subdivisions of the input range
    nb_subdiv_init = 1;
    /***/

    interval temp;
    int nb_points;
    
    inputs = vector<AAF>(jacdim);
    if (sysdim_params > 0)
        params = vector<AAF>(sysdim_params);
    
    uncontrolled = 0;
    controlled = 0;
    is_uncontrolled = vector<bool>(jacdim);
    variable = 0;
    is_variable = vector<bool>(jacdim);
    is_initialcondition = vector<bool>(jacdim);
    for (int i=0 ; i<jacdim; i++) {
        is_uncontrolled[i] = false;
        is_variable[i] = false;
        is_initialcondition[i] = false; // by definition, initial conditions are controlled
    }

    /**
     * Colin: fixed runtime arguments
     * Runs the HSCC 2014 Brusselator
     */
    t_begin = 0;
    tau = 0.05;
    t_end = 10.;
    order = 4;
    /***/

    /**
     * Colin: for initial inputs
     * x1(0), x2(0)
     */
    inputs[0] = interval(0.9,1);
    inputs[1] = interval(0,0.1);

    /**
     *
     *
     */
    params[0] = 1;
    params[1] = 1.5;
    nb_points = (t_end-t_begin)/tau+1;
    
    // common to EDO and DDE
    center_inputs = vector<AAF>(jacdim);
    eps = vector<interval>(jacdim);
    for (int i=0 ; i<jacdim ; i++)
    {
        if (is_uncontrolled[i])
            uncontrolled ++;
        if (!is_uncontrolled[i] && !is_initialcondition[i])
            controlled++;
        temp = inputs[i].convert_int();
        center_inputs[i] = mid(temp);
        eps[i] = temp-mid(temp);
    }
    
    cout << "controlled=" << controlled  << " uncontrolled=" << uncontrolled << endl;
    
    // for saving results
  //  cout << "(t_end-t_begin)*nb_subdiv/d0+1=" << ((t_end-t_begin)/d0+1)*(nb_subdiv+1) << endl;
    Xouter_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    Xouter_robust_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    Xouter_minimal_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    Xinner_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    Xinner_robust_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    Xinner_minimal_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    Xexact_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    t_print = vector<double>(nb_points);
        
}

void init_subdiv(int current_subdiv, vector<AAF> inputs_save, int param_to_subdivide)
{
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
    
   
     interval   temp = inputs[param_to_subdivide].convert_int();
        center_inputs[param_to_subdivide] = mid(temp);
        eps[param_to_subdivide] = temp-mid(temp);
    
}


// initial condition on [t0,t0+d0] given as x = g(t)
vector <T<AAF>> Initfunc(const  T<AAF> &t, vector<AAF> &beta)
{
    vector<T<AAF>> res(sysdim);

    // Colin: DDE only, should have no impact
//    // by default
//    for (int i=0 ; i<sysdim; i++)
//    res[i] = beta[i];
//
//    if (syschoice == 1)  // running example
//        res[0] = (1+beta[0]*t)*(1+beta[0]*t);  // ix[0] = beta
//    else if (syschoice == 2) //  example 5.15
//    {
//        res[0] = beta[0]*exp(t);
//        res[1] = beta[1]*(1-exp(-1));
//    }
//    else if (syschoice == 4) // Szczelina_1 2014
//    {
//        res[0] = beta[0]*sin(M_PI*t/2.0);
//    }
//    else if (syschoice == 5) // Szczelina_2 2014
//    {
//        res[0] = beta[0]*sin(M_PI*t/2.0);
//    }
    
    return res;
}



// initial condition on [-d0,0] given as x = g(t)
vector <T<F<AAF>>> Initfunc(const  T<F<AAF>> &t, vector<T<F<AAF>>> &beta)
{
    vector<T<F<AAF>>> res(sysdim);

    // Colin: DDE only, should have no impact
//    // by default
//    for (int i=0 ; i<sysdim; i++)
//    res[i] = beta[i];
//
//    if (syschoice == 1)  // running example
//        res[0] = (1+beta[0]*t)*(1+beta[0]*t);  // ix[0] = beta
//    else if (syschoice == 2) //  example 5.15
//    {
//        res[0] = beta[0]*exp(t);
//        res[1] = beta[1]*(1-exp(-1));
//    }
//    else if (syschoice == 4) // Szczelina_1 2014
//    {
//        res[0] = beta[0]*sin(M_PI*t/2.0);
//    }
//    else if (syschoice == 5) // Szczelina_2 2014
//    {
//        res[0] = beta[0]*sin(M_PI*t/2.0);
//    }
    return res;
}

// analytical solution if any (for comparison purpose)
vector <interval> AnalyticalSol(double t, vector<AAF> &beta, double d0)
{
    vector<interval> res(sysdim);
    vector<interval> Xouter_min(sysdim), Xouter_max(sysdim);
    
    double beta_inf = beta[0].convert_int().inf();
    double beta_sup = beta[0].convert_int().sup();
    
    // running example
    //res[0] = (1+beta[0]*t)*(1+beta[0]*t);  // ix[0] = beta
    
    if ((syschoice == 1) && beta_sup <= 1)  // running example
    {
        if (t <= 0)   // on [-d0,0], solution is defined by init function
        {
            Xouter_min[0] = ((1.+beta_inf*t)*(1.+beta_inf*t));
            Xouter_max[0] = ((1.+beta_sup*t)*(1.+beta_sup*t));
            res[0] = hull(Xouter_min[0],Xouter_max[0]);
        }
        else if (t >= 0 && t <= d0)
        {
            Xouter_min[0] = exp((-1./(3.*beta_inf)*(pow(1.+(t-1.)*beta_inf,3)-pow(1.-beta_inf,3))));
            Xouter_max[0] = exp((-1./(3.*beta_sup)*(pow(1.+(t-1.)*beta_sup,3)-pow(1.-beta_sup,3))));
            res[0] = hull(Xouter_min[0],Xouter_max[0]);
            
        }
        else if (t >= d0 && t <= 2*d0)
        {
            double aux, temp1, temp2;
            aux = exp((-1./(3.*beta_inf)*(pow(1.+(d0-1.)*beta_inf,3)-pow(1.-beta_inf,3))));
            temp1 = pow(1+(t-2)*beta_inf,3)/(3*beta_inf);
            temp2 = pow(1-beta_inf,3)/(3*beta_inf);
            Xouter_min[0] =  aux * exp(-exp(temp2)*pow(3*beta_inf,-2/3.0)*2.6789385347077476337*(gsl_sf_gamma_inc_P(1/3.0,temp1)-gsl_sf_gamma_inc_P(1/3.0,temp2)));
            aux = exp((-1./(3.*beta_sup)*(pow(1.+(d0-1.)*beta_inf,3)-pow(1.-beta_sup,3))));
            temp1 = pow(1+(t-2)*beta_sup,3)/(3*beta_sup);
            temp2 = pow(1-beta_sup,3)/(3*beta_sup);
            Xouter_max[0] =  aux * exp(-exp(temp2)*pow(3*beta_sup,-2/3.0)*2.6789385347077476337*(gsl_sf_gamma_inc_P(1/3.0,temp1)-gsl_sf_gamma_inc_P(1/3.0,temp2)));
            res[0] = hull(Xouter_min[0],Xouter_max[0]);
            // res[0] = interval(min(inf(Xouter_min[0]),inf(Xouter_max[0])),max(sup(Xouter_min[0]),sup(Xouter_max[0])));
          //  cout << "beta_inf = " << beta_inf << "beta_sup=" << beta_sup << " res[0] = " << res[0] << endl;
        }
        else
        {
            res[0] = interval(1,-1); // no analytical solution : bot
        }
        
    }
    else if (syschoice == 2) //  example 5.15
    {
        // example 5.15
        if (t <0)
        {
            res[0] = beta[0].convert_int()*exp(t);
            res[1] = beta[1].convert_int()*(1-exp(-1.));
        }
        else
        {
            res[0] = beta[0].convert_int()*exp(t);
            res[1] = beta[1].convert_int()*(exp(t)-exp(t-1.));
        }
    }
    else if (syschoice == 3) // Xue et al. 2017 ex. 3
    {
        for (int i=0 ; i<sysdim ; i++)
            res[i] = interval(1,-1); // no analytical solution : bot
    }
    else if (syschoice == 4 || syschoice == 5) // Szczelina_1 2014 : no analytical solution : bot
    {
        res[0] = interval(1,-1);
    }
    else if (syschoice == 6)  // self-driving car
    {
        for (int i=0 ; i<sysdim ; i++)
            res[i] = interval(1,-1); // no analytical solution : bot
    }
    else if (syschoice == 7)  // self-driving car
    {
        for (int i=0 ; i<sysdim ; i++)
            res[i] = interval(1,-1); // no analytical solution : bot
    }
    else if (syschoice == 8)  // self-driving car
    {
        for (int i=0 ; i<sysdim ; i++)
            res[i] = interval(1,-1); // no analytical solution : bot
    }
    else if (syschoice == 9) // Zou CAV 2015
    {
        res[0] = interval(1,-1);
    }
    else if (syschoice == 10 || syschoice == 11) // platoon
    {
        res[0] = interval(1,-1);
    }
    return res;
}




