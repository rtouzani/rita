/*==============================================================================

                                 r  i  t  a

            An environment for Modelling and Numerical Simulation

  ==============================================================================

    Copyright (C) 2021 Rachid Touzani

    This file is part of rita.

    rita is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    rita is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

  ==============================================================================

                            Definition of class 'rita'

  ==============================================================================*/

#ifndef __RITA_H
#define __RITA_H

#include <stdlib.h>
#include <math.h>
#include <string>
#include <time.h>

#include <iostream>
using std::ostream;
using std::cin;
using std::cout;
using std::endl;
using std::cerr;
using std::string;

#include "ritaException.h"
#include "OFELI.h"
#include "equa.h"


namespace RITA {
/*!
 *  \addtogroup RITA
 *  @{
 */


class cmd;
class data;
class configure;
class transient;
class stationary;
class optim;
class integration;
class approximation;
class help;
class solve;
class mesh;

#define CATCH                                                   \
   catch(OFELIException &e) {                                   \
      std::cout << "OFELI error: " << e.what() << endl;         \
      return 1;						\
   }                                                            \
   catch(runtime_error &e) {                                    \
      std::cout << "OFELI Runtime error: " << e.what() << endl; \
      return 1;                                                 \
   }                                                            \
   catch( ... ) {                                               \
      std::cout << "OFELI Unexpected error: " << endl;          \
      return 1;                                                 \
   }

struct odae {
   bool isSet, log, isFct;
   vector<string> analytic, vars;
   vector<OFELI::Fct> theFct;
   OFELI::Vect<double> y;
   OFELI::Vect<string> J;
   double init_time, final_time, time_step;
   OFELI::TimeScheme scheme;
   int size, field, ind_fct;
   NonLinearIter nls;
   string fn;
   odae();
   void setVars(int opt);
};

enum type {
   NONE,
   STEADY_STATE,
   TRANSIENT,
   EIGEN,
   OPTIMIZATION,
   APPROXIMATION,
   INTEGRATION
};

enum objective_type {
   ANALYTIC_FUNCTION,
   PDE_BASED
};


class rita
{

 public:

    rita();
    ~rita();
    int run();
    void setVerbose(int verb) { _verb = verb; }
    void setInput(string file, int opt=1);
    void initConfig();
    int setIn();
    bool meshOK, solveOK, dataOK;
    data *_data;

    friend class mesh;
    friend class solve;
    friend class transient;
    friend class stationary;
    friend class optim;
    friend class eigen;
    friend class integration;
    friend class approximation;
    friend class data;
    friend class equa;

 private:

   bool _load, _obj_analytic;
   odae *_ae, *_ode;
   equa *_pde;
   string _script_file, _scheme;
   ofstream *_ofh, *_ofl, _ocf;
   ifstream _icf, *_in;
   cmd *_cmd;
   int _verb, _key, _ret, _opt;
   mesh *_mesh;
   help *_help;
   configure *_configure;
   solve *_solve;
   transient *_transient;
   stationary *_stationary;
   optim *_optim;
   approximation *_approx;
   integration *_integration;
   double _init_time, _time_step, _final_time;
   int _adapted_time_step, _nb_eigv, _nb_args;
   bool _eigen_vectors, _default_field;
   std::vector<equa *> PDE;
   std::vector<odae *> ALGEBRAIC, ODE;
   OFELI::Mesh* _theMesh;
   bool _analysis_ok, _exit_ok;
   int _nb_eq, _nb_fields, _ieq, _ifield, _nb_ae, _nb_ode, _nb_pde;
   int _dim, _analysis_type;

   typedef void (rita::* Mode_Ptr)();
   static Mode_Ptr MODE[22];

   void set(OFELI::Mesh* ms) { _theMesh = ms; }
   void set(cmd* com) { _cmd = com; }
   void set(std::ofstream* ofl, std::ofstream* ofh) { _ofl = ofl; _ofh = ofh; }
   OFELI::Mesh* getMesh() const { return _theMesh; }
   void setDim(int dim) { _dim = dim; }
   void set(data *d) { _data = d; }
   int getNbEq() const { return _nb_eq; }
   int getNbFields() const { return _nb_fields; }
   int runAE();
   int runODE();
   void runPDE();

   void getHelp();
   void getLicense();
   void setConfigure();
   void Load();
   void UnLoad();
   void setData();
   void setMesh();
   void setStationary();
   void setTransient();
   void setEigen();
   void setOptim();
   void setApproximation();
   void setIntegration();
   void setAlgebraic();
   void setODE();
   void setPDE();
   void setSummary();
   void setClear();
   void Finish();
   void setMesh(OFELI::Mesh* ms);
   void setSolve();
   int set_ls(string ls, string prec);
   int set_nls(string nls);
   int findField(const string& s);

   const vector<string> _kw = {"?","help","lic$ense","set","load","unload","data","mesh",
                               "stat$ionary","trans$ient","eigen","optim","approx$imation",
                               "integ$ration","algebraic","ode","pde","summary","solve","clear",
                               "exit","quit","EXIT"};
   map<string,OFELI::Iteration> Ls = {{"direct",OFELI::DIRECT_SOLVER},
                                      {"cg",OFELI::CG_SOLVER},
                                      {"cgs",OFELI::CGS_SOLVER},
                                      {"bicg",OFELI::BICG_SOLVER},
                                      {"bicg_stab",OFELI::BICG_STAB_SOLVER},
                                      {"gmres",OFELI::GMRES_SOLVER}};
   map<OFELI::Iteration,string> rLs = {{OFELI::DIRECT_SOLVER,"direct"},
                                       {OFELI::CG_SOLVER,"cg"},
                                       {OFELI::CGS_SOLVER,"cgs"},
                                       {OFELI::BICG_SOLVER,"bicg"},
                                       {OFELI::BICG_STAB_SOLVER,"bicg_stab"},
                                       {OFELI::GMRES_SOLVER,"gmres"}};
   map<std::string,NonLinearIter> NLs = {{"bisection",BISECTION},
                                         {"regula-falsi",REGULA_FALSI},
                                         {"picard",PICARD},
                                         {"secant",SECANT},
                                         {"newton",NEWTON}};
   map<std::string,OFELI::TimeScheme> _sch = {{"forward-euler",OFELI::FORWARD_EULER},
                                              {"backward-euler",OFELI::BACKWARD_EULER},
                                              {"crank-nicolson",OFELI::CRANK_NICOLSON},
                                              {"heun",OFELI::HEUN},
                                              {"newmark",OFELI::NEWMARK},
                                              {"leap-frog",OFELI::LEAP_FROG},
                                              {"AB2",OFELI::ADAMS_BASHFORTH},
                                              {"RK4",OFELI::RK4},
                                              {"RK3-TVD",OFELI::RK3_TVD},
                                              {"BDF2",OFELI::BDF2},
                                              {"builtin",OFELI::BUILTIN}};
   map<string,OFELI::Preconditioner> Prec = {{"ident",OFELI::IDENT_PREC},
                                             {"diag",OFELI::DIAG_PREC},
                                             {"dilu",OFELI::DILU_PREC},
                                             {"ilu",OFELI::ILU_PREC},
                                             {"ssor",OFELI::SSOR_PREC}};
   map<OFELI::Preconditioner,string> rPrec = {{OFELI::IDENT_PREC,"ident"},
                                              {OFELI::DIAG_PREC,"diag"},
                                              {OFELI::DILU_PREC,"dilu"},
                                              {OFELI::ILU_PREC,"ilu"},
                                              {OFELI::SSOR_PREC,"ssor"}};
   vector<int> _eq_type;
   int setSpaceDiscretization(string& sp);
};

} /* namespace RITA */

#endif
