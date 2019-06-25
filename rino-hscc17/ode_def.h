/* ============================================================================
 File   : ode_def.h
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
 The place to declare the systems of ODEs or DDEs on which he/she wants to perform reachability
 Setup of dimension of the system (sysdim), initial conditions and parameters is in ode_def.cpp
 ============================================================================ */
#ifndef ODE_DEF_H
#define ODE_DEF_H

#include "filib_interval.h"
#include "tadiff.h" 
#include "fadiff.h"
#include "fadbad_aa.h"

using namespace std;

#define innerapprox 1

extern int systype;    // systype = 0 (ODE) or 1 (DDE) -- initialized in main.cpp / command line
extern int syschoice;  // to choose among the predefined systems of ODE or DDE -- initialized in main.cpp / command line

// dimension of the system of ODE/DDE to analyze:
extern int sysdim;
// dimension of uncertain input
extern int jacdim; // Jacobian will be dimension sysdim * jacdim
extern int sysdim_params;  // dimension of the vector of parameters params - Boost
extern vector<AAF> params;      // params of the ODE (nondeterministic disturbances)
extern vector<AAF> inputs;   // uncertain inputs and parameters : some will be used in initial condition, some as uncertain parameters
extern vector<AAF> center_inputs;
extern vector<interval> eps;

// for subdivisions of the initial domain to refine precision
extern int nb_subdiv_init; // number of subdivisiions
extern double recovering; // percentage of recovering between subdivisions
extern vector<vector<vector<interval>>> Xouter_print, Xouter_robust_print, Xouter_minimal_print, Xinner_print, Xinner_robust_print, Xinner_minimal_print, Xexact_print; // store results of subdivision
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

// for ODEs : initialize the state variable (and center for inner-approximation)
void set_initialconditions(vector<AAF> &x, vector<AAF> &xcenter, vector<vector<AAF>> &J);

// for DDEs : functions that initialize the DDE on first time period
vector<T<AAF>> Initfunc(const T<AAF>& t, vector<AAF> &beta);
vector <T<F<AAF>>> Initfunc(const  T<F<AAF>> &t, vector<T<F<AAF>>> &beta);

// defining analytical solution if any for comparison
vector <interval> AnalyticalSol(double t, vector<AAF> &beta, double d0);


// for ODEs and DDEs: define bounds for parameters and inputs, value of delay d0 if any, and parameters of integration (timestep, order of TM)
void init_system(double &t_begin, double &t_end, double &tau, double &d0, int &nb_subdiv, int &order);

// specific to subdivisions
void init_subdiv(int current_subdiv, vector<AAF> inputs_save, int param_to_subdivide);



// define here  your ODE system yp = \dot y = f(y)
class OdeFunc {
    public:
        template <class C>
        void operator()(vector<C> &yp, vector<C> y) {

            /**
             *
             * Colin This is hardcoded to be the Brusselator
             * Where is yp stored?
             */
            yp[0] = 1 - (params[1]+1)*y[0] + params[0]*y[0]*y[0]*y[1];
            yp[1] = params[1]*y[0] - params[0]*y[0]*y[0]*y[1];
    }
};


// define here  your DEE system yp = \dot y = f(y, y_prev, params)
class DdeFunc {
public:
    template <class C>
    void operator()(vector<C> &yp, vector<C> &y, vector<C> &y_prev, vector<AAF> &param_inputs) {
        yp[0] = y_prev[0] + y[1];
        yp[1] = y[0] - y_prev[0];
    }
};

// for DDE - variational equations for Jacobian with respect to the initial conditions have to be defined
class DdeJacFunc {
public:
    template <class C>
    //  void operator()(C yp[sysdim], C y[sysdim], int mode) {
    void operator()(vector<vector<C>> &Jp, vector<vector<C>> &J, vector<vector<C>> &J_prev, vector<C> &x,  vector<C> &x_prev) {
        // running example
        if (syschoice == 1)  // running example
        {
            Jp[0][0] = -J[0][0]*x_prev[0] - J_prev[0][0]*x[0];
        }
        else if (syschoice == 2) //  example 5.15
        {
            Jp[0][0] = J_prev[0][0] + J[1][0];
            Jp[0][1] = J_prev[0][1] + J[1][1];
            Jp[1][0] = J[0][0] - J_prev[0][0];
            Jp[1][1] = J[0][1] - J_prev[0][1];
        }
        else if (syschoice == 3) // Xue et al. 2017 ex 3
        {
            Jp[0][0] = 1.4*J[2][0] - 0.9*J_prev[0][0];
            Jp[0][1] = 1.4*J[2][1] - 0.9*J_prev[0][1];
            Jp[0][2] = 1.4*J[2][2] - 0.9*J_prev[0][2];
            Jp[0][3] = 1.4*J[2][3] - 0.9*J_prev[0][3];
            Jp[0][4] = 1.4*J[2][4] - 0.9*J_prev[0][4];
            Jp[0][5] = 1.4*J[2][5] - 0.9*J_prev[0][5];
            Jp[0][6] = 1.4*J[2][6] - 0.9*J_prev[0][6];
            
            Jp[1][0] = 2.5*J[4][0] - 1.5*J[1][0];
            Jp[1][1] = 2.5*J[4][1] - 1.5*J[1][1];
            Jp[1][2] = 2.5*J[4][2] - 1.5*J[1][2];
            Jp[1][3] = 2.5*J[4][3] - 1.5*J[1][3];
            Jp[1][4] = 2.5*J[4][4] - 1.5*J[1][4];
            Jp[1][5] = 2.5*J[4][5] - 1.5*J[1][5];
            Jp[1][6] = 2.5*J[4][6] - 1.5*J[1][6];
            
            Jp[2][0] = 0.6*J[6][0] - 0.8*J[2][0]*x[1] - 0.8*x[2]*J[1][0];
            Jp[2][1] = 0.6*J[6][1] - 0.8*J[2][1]*x[1] - 0.8*x[2]*J[1][1];
            Jp[2][2] = 0.6*J[6][2] - 0.8*J[2][2]*x[1] - 0.8*x[2]*J[1][2];
            Jp[2][3] = 0.6*J[6][3] - 0.8*J[2][3]*x[1] - 0.8*x[2]*J[1][3];
            Jp[2][4] = 0.6*J[6][4] - 0.8*J[2][4]*x[1] - 0.8*x[2]*J[1][4];
            Jp[2][5] = 0.6*J[6][5] - 0.8*J[2][5]*x[1] - 0.8*x[2]*J[1][5];
            Jp[2][6] = 0.6*J[6][6] - 0.8*J[2][6]*x[1] - 0.8*x[2]*J[1][6];
            
            Jp[3][0] = - 1.3*J[3][0]*x[2] - 1.3*x[3]*J[2][0];
            Jp[3][1] = - 1.3*J[3][1]*x[2] - 1.3*x[3]*J[2][1];
            Jp[3][2] = - 1.3*J[3][2]*x[2] - 1.3*x[3]*J[2][2];
            Jp[3][3] = - 1.3*J[3][3]*x[2] - 1.3*x[3]*J[2][3];
            Jp[3][4] = - 1.3*J[3][4]*x[2] - 1.3*x[3]*J[2][4];
            Jp[3][5] = - 1.3*J[3][5]*x[2] - 1.3*x[3]*J[2][5];
            Jp[3][6] = - 1.3*J[3][6]*x[2] - 1.3*x[3]*J[2][6];
            
            Jp[4][0] = 0.7*J[0][0] - J[3][0]*x[4] - x[3]*J[4][0];
            Jp[4][1] = 0.7*J[0][1] - J[3][1]*x[4] - x[3]*J[4][1];
            Jp[4][2] = 0.7*J[0][2] - J[3][2]*x[4] - x[3]*J[4][2];
            Jp[4][3] = 0.7*J[0][3] - J[3][3]*x[4] - x[3]*J[4][3];
            Jp[4][4] = 0.7*J[0][4] - J[3][4]*x[4] - x[3]*J[4][4];
            Jp[4][5] = 0.7*J[0][5] - J[3][5]*x[4] - x[3]*J[4][5];
            Jp[4][6] = 0.7*J[0][6] - J[3][6]*x[4] - x[3]*J[4][6];
            
            Jp[5][0] = 0.3*J[0][0] - 3.1*J[5][0];
            Jp[5][1] = 0.3*J[0][1] - 3.1*J[5][1];
            Jp[5][2] = 0.3*J[0][2] - 3.1*J[5][2];
            Jp[5][3] = 0.3*J[0][3] - 3.1*J[5][3];
            Jp[5][4] = 0.3*J[0][4] - 3.1*J[5][4];
            Jp[5][5] = 0.3*J[0][5] - 3.1*J[5][5];
            Jp[5][6] = 0.3*J[0][6] - 3.1*J[5][6];

            
            Jp[6][0] = 1.8*J[5][0] - 1.5*J[6][0]*x[1] - 1.5*x[6]*J[1][0];
            Jp[6][1] = 1.8*J[5][1] - 1.5*J[6][1]*x[1] - 1.5*x[6]*J[1][1];
            Jp[6][2] = 1.8*J[5][2] - 1.5*J[6][2]*x[1] - 1.5*x[6]*J[1][2];
            Jp[6][3] = 1.8*J[5][3] - 1.5*J[6][3]*x[1] - 1.5*x[6]*J[1][3];
            Jp[6][4] = 1.8*J[5][4] - 1.5*J[6][4]*x[1] - 1.5*x[6]*J[1][4];
            Jp[6][5] = 1.8*J[5][5] - 1.5*J[6][5]*x[1] - 1.5*x[6]*J[1][5];
            Jp[6][6] = 1.8*J[5][6] - 1.5*J[6][6]*x[1] - 1.5*x[6]*J[1][6];
        }
        else if (syschoice == 4) // Szczelina_1 2014
        {
            Jp[0][0] = -J[0][0] + J_prev[0][0]*x_prev[0];
        }
        else if (syschoice == 5) // Szczelina_2 2014
        {
            Jp[0][0] = -J[0][0] -3.2*J_prev[0][0] + 3.*J_prev[0][0]*x_prev[0]*x_prev[0];
        }
        else if (syschoice == 6)  // self-driving car
        {
            Jp[0][0] = J[1][0];
            Jp[0][1] = J[1][1];
            Jp[1][0] = - params[0]*J_prev[0][0] - params[1]*J_prev[1][0];
            Jp[1][1] = - params[0]*J_prev[0][1] - params[1]*J_prev[1][1];
        }
        else if (syschoice == 7)  // self-driving car with uncertain (but constant) coefficients
        {
            Jp[0][0] = J[1][0];
            Jp[0][1] = J[1][1];
            Jp[0][2] = J[1][2];
            Jp[0][3] = J[1][3];
            Jp[1][0] = -J[2][0]*x_prev[0] - x[2] *J_prev[0][0] + J[2][0] - J[3][0]*x_prev[1] - x[3]*J_prev[1][0];
            Jp[1][1] = -J[2][1]*x_prev[0] - x[2] *J_prev[0][1] + J[2][1]  - J[3][1]*x_prev[1] - x[3]*J_prev[1][1];
            Jp[1][2] = -J[2][2]*x_prev[0] - x[2] *J_prev[0][2] + J[2][2]  - J[3][2]*x_prev[1] - x[3]*J_prev[1][2];
            Jp[1][3] = -J[2][3]*x_prev[0] - x[2] *J_prev[0][3] + J[2][3]  - J[3][3]*x_prev[1] - x[3]*J_prev[1][3];
            for (int i=0 ; i<sysdim ; i++)
            {
                Jp[2][i] = 0;
                Jp[3][i] = 0;
            }
        }
        else if (syschoice == 8)  // self-driving car with uncertain (but constant) coefficients
        {
            Jp[0][0] = J[1][0];
            Jp[0][1] = J[1][1];
            Jp[0][2] = J[1][2];
            Jp[0][3] = J[1][3];
            Jp[1][0] = - inputs[2]*J_prev[0][0]   - inputs[3]*J_prev[1][0];
            Jp[1][1] = - inputs[2]*J_prev[0][1]   - inputs[3]*J_prev[1][1];
            Jp[1][2] = -x_prev[0] - inputs[2]*J_prev[0][2] + 1 - inputs[3]*J_prev[1][2];
            Jp[1][3] = - inputs[2]*J_prev[0][3] - x_prev[1] - inputs[3]*J_prev[1][3];
        }
        else if (syschoice == 9) // Ex 4 of Zou CAV 2015
        {
            Jp[0][0] = -3*x_prev[0]*x_prev[0]*J_prev[0][0];
        }
        else if (syschoice == 10) // platoon
        {
            double a = 2.5;
            //       yp[0] = 2 + (y[0]/5.-1)*(y[0]/5.-2)*(y[0]/5.-3)/6.; // x1' : '2+(x1/5-1)*(x1/5-2)*(x1/5-3)/6',
            Jp[0][0] = ((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][0]/30;
            Jp[0][1] = ((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][1]/30;
            Jp[0][2] = ((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][2]/30;
            Jp[0][3] = ((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][3]/30;
            Jp[0][4] = ((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][4]/30;
            Jp[0][5] = ((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][5]/30;
            Jp[0][6] = ((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][6]/30;
            Jp[0][7] = ((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][7]/30;
            Jp[0][8] = ((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][8]/30;
            //      yp[1] = y[2]; //  'x2' : 'dx2',
            Jp[1][0] = J[2][0];
            Jp[1][1] = J[2][1];
            Jp[1][2] = J[2][2];
            Jp[1][3] = J[2][3];
            Jp[1][4] = J[2][4];
            Jp[1][5] = J[2][5];
            Jp[1][6] = J[2][6];
            Jp[1][7] = J[2][7];
            Jp[1][8] = J[2][8];
            //      yp[2] = a*(2+(y[0]/5.-1)*(y[0]/5.-2)*(y[0]/5.-3)/6. - y_prev[2]); // 'dx2' : 'a*(2+(x1/5-1)*(x1/5-2)*(x1/5-3)/6-dx2(t-tau))',
            Jp[2][0] = a*(((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][0]/30 - J_prev[2][0]);
            Jp[2][1] = a*(((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][1]/30 - J_prev[2][1]);
            Jp[2][2] = a*(((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][2]/30 - J_prev[2][2]);
            Jp[2][3] = a*(((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][3]/30 - J_prev[2][3]);
            Jp[2][4] = a*(((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][4]/30 - J_prev[2][4]);
            Jp[2][5] = a*(((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][5]/30 - J_prev[2][5]);
            Jp[2][6] = a*(((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][6]/30 - J_prev[2][6]);
            Jp[2][7] = a*(((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][7]/30 - J_prev[2][7]);
            Jp[2][8] = a*(((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][8]/30 - J_prev[2][8]);
             //      yp[3] = y[4]; //  'x3' : 'dx3',
            Jp[3][0] = J[4][0];
            Jp[3][1] = J[4][1];
            Jp[3][2] = J[4][2];
            Jp[3][3] = J[4][3];
            Jp[3][4] = J[4][4];
            Jp[3][5] = J[4][5];
            Jp[3][6] = J[4][6];
            Jp[3][7] = J[4][7];
            Jp[3][8] = J[4][8];
            //     yp[4] = a*(y_prev[2]-y_prev[4]); //   'dx3' : 'a*(dx2(t-tau)-dx3(t-tau))',
            Jp[4][0] = a*(J_prev[2][0]-J_prev[4][0]);
            Jp[4][1] = a*(J_prev[2][1]-J_prev[4][1]);
            Jp[4][2] = a*(J_prev[2][2]-J_prev[4][2]);
            Jp[4][3] = a*(J_prev[2][3]-J_prev[4][3]);
            Jp[4][4] = a*(J_prev[2][4]-J_prev[4][4]);
            Jp[4][5] = a*(J_prev[2][5]-J_prev[4][5]);
            Jp[4][6] = a*(J_prev[2][6]-J_prev[4][6]);
            Jp[4][7] = a*(J_prev[2][7]-J_prev[4][7]);
            Jp[4][8] = a*(J_prev[2][8]-J_prev[4][8]);
             //     yp[5] = y[6]; //   'x4' : 'dx4',
            Jp[5][0] = J[6][0];
            Jp[5][1] = J[6][1];
            Jp[5][2] = J[6][2];
            Jp[5][3] = J[6][3];
            Jp[5][4] = J[6][4];
            Jp[5][5] = J[6][5];
            Jp[5][6] = J[6][6];
            Jp[5][7] = J[6][7];
            Jp[5][8] = J[6][8];
            //     yp[6] = a*(y_prev[4]-y_prev[6]);//   'dx4' : 'a*(dx3(t-tau)-dx4(t-tau))',
            Jp[6][0] = a*(J_prev[4][0]-J_prev[6][0]);
            Jp[6][1] = a*(J_prev[4][1]-J_prev[6][1]);
            Jp[6][2] = a*(J_prev[4][2]-J_prev[6][2]);
            Jp[6][3] = a*(J_prev[4][3]-J_prev[6][3]);
            Jp[6][4] = a*(J_prev[4][4]-J_prev[6][4]);
            Jp[6][5] = a*(J_prev[4][5]-J_prev[6][5]);
            Jp[6][6] = a*(J_prev[4][6]-J_prev[6][6]);
            Jp[6][7] = a*(J_prev[4][7]-J_prev[6][7]);
            Jp[6][8] = a*(J_prev[4][8]-J_prev[6][8]);
             //     yp[7] = y[8];//  'x5' : 'dx5',
            Jp[7][0] = J[8][0];
            Jp[7][1] = J[8][1];
            Jp[7][2] = J[8][2];
            Jp[7][3] = J[8][3];
            Jp[7][4] = J[8][4];
            Jp[7][5] = J[8][5];
            Jp[7][6] = J[8][6];
            Jp[7][7] = J[8][7];
            Jp[7][8] = J[8][8];
             //      yp[8] = a*(y_prev[6]-y[8]);//  'dx5' : 'a*(dx4(t-tau)-dx5)'
            Jp[8][0] = a*(J_prev[6][0]-J[8][0]);
            Jp[8][1] = a*(J_prev[6][1]-J[8][1]);
            Jp[8][2] = a*(J_prev[6][2]-J[8][2]);
            Jp[8][3] = a*(J_prev[6][3]-J[8][3]);
            Jp[8][4] = a*(J_prev[6][4]-J[8][4]);
            Jp[8][5] = a*(J_prev[6][5]-J[8][5]);
            Jp[8][6] = a*(J_prev[6][6]-J[8][6]);
            Jp[8][7] = a*(J_prev[6][7]-J[8][7]);
            Jp[8][8] = a*(J_prev[6][8]-J[8][8]);
        }
        else if (syschoice == 11) // platoon of 10 vehicles
        {
            double a = 2.5;
           
            for (int j=0; j<jacdim ; j++)
            {
                 //       yp[0] = 2 + (y[0]/5.-1)*(y[0]/5.-2)*(y[0]/5.-3)/6.; // x1' : '2+(x1/5-1)*(x1/5-2)*(x1/5-3)/6',
                Jp[0][j] = ((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][j]/30;
             //      yp[1] = y[2]; //  'x2' : 'dx2',
                Jp[1][j] = J[2][j];
            //      yp[2] = a*(2+(y[0]/5.-1)*(y[0]/5.-2)*(y[0]/5.-3)/6. - y_prev[2]); // 'dx2' : 'a*(2+(x1/5-1)*(x1/5-2)*(x1/5-3)/6-dx2(t-tau))',
                Jp[2][j] = a*(((x[0]/5-2)*(2*x[0]/5-4) +(x[0]/5-1)*(x[0]/5-3))*J[0][j]/30 - J_prev[2][j]);
            //      yp[3] = y[4]; //  'x3' : 'dx3',
                Jp[3][j] = J[4][j];
            //     yp[4] = a*(y_prev[2]-y_prev[4]); //   'dx3' : 'a*(dx2(t-tau)-dx3(t-tau))',
                Jp[4][j] = a*(J_prev[2][j]-J_prev[4][j]);
            //     yp[5] = y[6]; //   'x4' : 'dx4',
                Jp[5][j] = J[6][j];
            //     yp[6] = a*(y_prev[4]-y_prev[6]);//   'dx4' : 'a*(dx3(t-tau)-dx4(t-tau))',
                Jp[6][j] = a*(J_prev[4][j]-J_prev[6][j]);
            //     yp[7] = y[8];//  'x5' : 'dx5',
                Jp[7][j] = J[8][j];
            //      yp[8] = a*(y_prev[6]-y[8]);//  'dx5' : 'a*(dx4(t-tau)-dx5)'
                Jp[8][j] = a*(J_prev[6][j]-J_prev[8][j]);
                //
                Jp[9][j] = J[10][j];
                Jp[10][j] = a*(J_prev[8][j]-J_prev[10][j]);
                Jp[11][j] = J[12][j];
                Jp[12][j] = a*(J_prev[10][j]-J_prev[12][j]);
                Jp[13][j] = J[14][j];
                Jp[14][j] = a*(J_prev[12][j]-J_prev[14][j]);
                Jp[15][j] = J[16][j];
                Jp[16][j] = a*(J_prev[14][j]-J_prev[16][j]);
                Jp[17][j] = J[18][j];
                Jp[18][j] = a*(J_prev[16][j]-J[18][j]);
            }
        }
    }
};






#endif
