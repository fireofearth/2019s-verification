/* ============================================================================
 File   : ode_def.h
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
 The place to declare the systems of ODEs or DDEs on which he/she wants to perform reachability
 Setup of dimension of the system (sysdim), initial conditions and parameters is in ode_def.cpp
 ============================================================================ */
#ifndef ODE_DEF_H
#define ODE_DEF_H

#include <exception>
#include "filib_interval.h"
#include "tadiff.h" 
#include "fadiff.h"
#include "fadbad_aa.h"

using namespace std;

#define innerapprox 1

enum Problem {BRUSSELATOR, VANDERPOLL};

// dimension of the system of ODE/DDE to analyze:
extern int sysdim;
// dimension of uncertain input
extern int jacdim; // Jacobian will be dimension sysdim * jacdim
extern vector<AAF> inputs;   // uncertain inputs and parameters : some will be used in initial condition, some as uncertain parameters
extern vector<AAF> center_inputs;
extern vector<interval> eps;

// for subdivisions of the initial domain to refine precision
extern int nb_subdiv_init; // number of subdivisiions
extern double recovering; // percentage of recovering between subdivisions
extern vector<vector<vector<interval>>> Xouter_print, Xouter_robust_print, Xouter_minimal_print, Xinner_print,
    Xinner_robust_print, Xinner_minimal_print, Xexact_print; // store results of subdivision
extern vector<double> t_print; // times where results are stored
extern int current_subdiv;
extern int current_iteration;

// for robust inner-approximations
extern int uncontrolled;  // number of uncontrolled parameters (forall params)
extern int controlled;  // number of controlled parameters (forall params)
extern vector<bool> is_uncontrolled; // for each input, uncontrolled or controlled (for robust inner-approx)
extern vector<bool> is_initialcondition; // for each input, initial condition or parameter (for robust inner-approx)
extern int variable;  // number of non constant parameters
extern vector<bool> is_variable; // for each parameter, constant or variable

struct OdeInit {
    int sysdim;
    int jacdim;
    int nb_subdiv_init;
    double tau;
    double t_begin;
    double t_end;
    int order;
    vector<AAF> inputs;
};

class OdeFunc {
private:
    Problem p;
public:
    OdeFunc(Problem p) {
        this-> p = p;
    }

    struct OdeInit get_initial_variables() {
        switch(this->p) {
            case BRUSSELATOR:
                return { 2, 2, 1, 0.05, 0, 10., 4,
                         { interval(0.9,1), interval(0,0.1) }
                };
            case VANDERPOLL:
                throw logic_error("incomplete functionality");
                break;
            default:
                throw logic_error("problem is not specified");
        }
    }

    template <class C>
    void operator()(vector<C> &yp, vector<C> y) {
        /**
         * Colin: This is hardcoded to be the Brusselator
         */
        switch(this->p) {
            case BRUSSELATOR:
                yp[0] = 1 - (1.5 + 1) * y[0] + y[0] * y[0] * y[1];
                yp[1] = 1.5 * y[0] - y[0] * y[0] * y[1];
                break;
            case VANDERPOLL:
                throw logic_error("incomplete functionality");
                break;
            default:
                throw logic_error("problem is not specified");
        }
    }
};

// for ODEs : initialize the state variable (and center for inner-approximation)
void set_initialconditions(vector<AAF> &x, vector<AAF> &xcenter, vector<vector<AAF>> &J);

// for ODEs and DDEs: define bounds for parameters and inputs, value of delay d0 if any, and parameters of integration (timestep, order of TM)
void init_system(OdeFunc &odef, double &t_begin, double &t_end, double &tau, double &d0, int &nb_subdiv, int &order);

// specific to subdivisions
void init_subdiv(int current_subdiv, vector<AAF> inputs_save, int param_to_subdivide);

#endif