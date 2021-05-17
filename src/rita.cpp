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

                            Main calling program and
                          Implementation of class 'rita'

  ==============================================================================*/

#include <fstream>
#include <ctime>

#include "rita.h"
#include "data.h"
#include "cmd.h"
#include "mesh.h"
#include "equa.h"
#include "solve.h"
#include "optim.h"
#include "integration.h"
#include "approximation.h"
#include "help.h"
#include "configure.h"

using std::cout;
using std::exception;

const string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf,sizeof(buf),"%Y-%m-%d  %X",&tstruct);
    return buf;
}


int main(int argc, char *argv[])
{
   std::time_t t = std::time(0);
   std::tm* now = std::localtime(&t);
   cout << "\n     R I T A     1.0\n";
   cout << "     Last update " << (now->tm_year+1900) << '-' << (now->tm_mon+1)
        << '-' <<  now->tm_mday << "\n";
   cout << "     Copyright (C) 2021, Rachid Touzani\n\n";
   cout << "Type \"help\" or \"license\" for more information." << endl;
   cout << "rita web site:   http://ofeli.org/rita\n" << endl;
   int ret = 0;

   try {
      RITA::rita r;
      if (argc>1) {
         if (string(argv[1])=="-h") {
            cout << RITA::H0;
            return 0;
         }
         else if (string(argv[1])=="--help") {
            cout << RITA::H1;
            return 0;
         }
         else if (string(argv[1])=="-v" || string(argv[1])=="--version") {
            cout << "rita, Release 1.0" << endl;
            return 0;
         }
         else if (argv[1][0]=='-') {
            cout << "Input error: Unknown argument." << endl;
            return 1;
         }
         else
            r.setInput(string(argv[1]),0);
      }
      ret = r.run();
   } CATCH_RITA_EXCEPTION
   return ret;
}


namespace RITA {

rita::Mode_Ptr rita::MODE [] = {
      nullptr,
      &rita::getHelp,
      &rita::getHelp,
      &rita::getLicense,
      &rita::setConfigure,
      &rita::Load,
      &rita::UnLoad,
      &rita::setData,
      &rita::setMesh,
      &rita::setStationary,
      &rita::setTransient,
      &rita::setEigen,
      &rita::setOptim,
      &rita::setApproximation,
      &rita::setIntegration,
      &rita::setAlgebraic,
      &rita::setODE,
      &rita::setPDE,
      &rita::setSummary,
      &rita::setSolve,
      &rita::setClear
};


rita::rita()
     : meshOK(false), solveOK(false), dataOK(false), _load(false), _ae(nullptr),
       _ode(nullptr), _pde(nullptr), _script_file(""), _in(nullptr), _verb(1), _ret(0),
       _default_field(true), _theMesh(nullptr), _exit_ok(false), _nb_eq(0),
       _nb_fields(0), _ieq(0), _ifield(0), _nb_ae(0), _nb_ode(0), _nb_pde(0),
       _analysis_type(NONE)
{
   _cmd = new cmd;
   _configure = new configure(this,_cmd);
   _mesh = new mesh(this,_cmd,_configure);
   _data = new data(this,_cmd,_configure);
   _help = new help;
   _solve = new solve(this,_cmd,_configure);
   ofl = _configure->getOStreamLog();
   ofh = _configure->getOStreamHistory();
   _solve->setSave(_configure->getSaveResults());
   _optim = new optim(this,_cmd,_configure);
   _integration = new integration(this,_cmd,_configure);
   _approx = new approximation(this,_cmd,_configure);
   _init_time = 0.;
   _final_time = 1.;
   _time_step = 0.1;
   _adapted_time_step = 0;
   _scheme = "backward-euler";
}


rita::~rita()
{
   if (_in!=nullptr)
      delete _in;
   delete _cmd;
   delete _mesh;
   delete _data;
   delete _solve;
   delete _optim;
   delete _integration;
   delete _help;
   delete _configure;
   if (_ae!=nullptr)
      delete _ae;
   if (_ode!=nullptr)
      delete _ode;
   if (_pde!=nullptr)
      delete _pde;
}


void rita::setInput(string file,
                    int    opt)
{
   _script_file = file;
   _opt = opt;
   Load();
}


int rita::run()
{
   while (1) {
      if (_ret>=100) {
         if (_verb>1)
            cout << "Terminating rita!" << endl;
         if (!_exit_ok) {
            *ofh << "exit" << endl;
            _exit_ok = true;
         }
         break;
      }
      if (_cmd->readline("rita> "))
         continue;
      _nb_args = _cmd->getNbArgs();
      _key = _cmd->getKW(_kw);
      if (_key<-1)
         continue;
      else if (_key==-1) {
         msg("","Unknown command: "+_cmd->token(),
             "Available commands: license, load, unload, data, mesh, stationary, transient\n"
             "                    eigen, optim, approximation, integration, algebraic, ode, pde\n"
             "                    summary, solve, clear\n"
             "Global commands:    help, ?, set, quit, exit");
         continue;
      }
      else if (_key>19) {
         if (!_exit_ok) {
            *ofh << "exit" << endl;
            _exit_ok = true;
            if (_key==23)
               _ret = 200;
         }
         break;
      }
      (this->*MODE[_key+1])();
      if (_ret>=100)
         break;
   }
   if (_script_file=="") {
      if (_ret>100)
         return 0;
      cout << "Are you sure (y/n) ? ";
      string ans="y";
      cin >> ans;
      if (ans!="y")
         return 0;
   }
   if (_verb>1)
      cout << "Terminating rita!" << endl;
   if (!_exit_ok) {
      *ofh << "exit" << endl;
      _exit_ok = true;
   }
   exit(0);
   return 1;
}


void rita::getLicense()
{
   _cmd->setNbArg(0);
   cout << "Copyright (C) 2021, Rachid Touzani\n";
   cout << "rita is free software: you can redistribute it and/or modify\n";
   cout << "it under the terms of the GNU General Public License as published by\n";
   cout << "the Free Software Foundation, either version 3 of the License, or\n";
   cout << "(at your option) any later version.\n\n";
   cout << "rita is distributed in the hope that it will be useful,\n";
   cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of\n";
   cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n";
   cout << "GNU General Public License for more details." << endl;
}


void rita::getHelp()
{
   _cmd->setNbArg(0);
   _help->setVerbose(_verb);
   _help->run();
}


void rita::setConfigure()
{
   _configure->run();
   _verb = _configure->getVerbose();
   _ret = _configure->ret();
}


void rita::Load()
{
   if (_script_file.size()==0) {
      if (_cmd->setNbArg(1,"Give input script file.")) {
         msg("load>","Missing script file to load.","",1);
         return;
      }
      _cmd->get(_script_file);
   }
   if (_in!=nullptr)
      delete _in, _in = nullptr;
   _in = new ifstream(_script_file);
   if (_in->is_open())
      _cmd->setIFStream(_in);
   else {
      msg("load>","Unable to open file: "+_script_file);
      if (_opt)
         _ret = 1;
      else
         _ret = 100;
      return;
   }
   _load = true;
   run();
}


void rita::UnLoad()
{
   if (_in!=nullptr)
      delete _in;
   _in = nullptr;
//   run();
}


void rita::setData()
{
   if (_verb>1)
      cout << "Entering mode 'data' ..." << endl;
   _cmd->setNbArg(0);
   _data->setVerbose(_verb);
   _ret = _data->run();
   _nb_fields = _data->getNbFields();
   if (_ret==0)
      dataOK = true;
   if (_ret==10) {
      dataOK = true;
      setMesh();
   }
   if (_verb>1)
      cout << "Leaving mode 'data' ..." << endl;
}


void rita::setMesh()
{
   if (_verb>1)
      cout << "Entering mode 'mesh' ..." << endl;
   _cmd->setNbArg(0);
   *ofh << _kw[_key] << endl;
   _mesh->set(_cmd);
   _mesh->setVerbose(_verb);
   _theMesh = new OFELI::Mesh;
   _ret = _mesh->run();
   _theMesh = _mesh->get();
   if (_theMesh)
      meshOK = true;
   if (_verb>1)
      cout << "Leaving mode 'mesh' ..." << endl;
}


void rita::setStationary()
{
   _analysis_type = STEADY_STATE;
}


void rita::setTransient()
{
   int ret = 0;
   _analysis_type = TRANSIENT;
   _scheme = "backward-euler";
   _adapted_time_step = 0;
   double it=_init_time, ft=_final_time, ts=_time_step;
   _ret = 0;
   static const string H = "transient [initial-time=it] [final-time=ft] [time-step=ts]  [scheme=s] [adapted]\n\n"
                           "it: Initial value of time. Default value is 0.\n"
                           "ft: Final (maximal) value of time. Default value is 1.\n"
                           "ts: Time step value. Default value is 0.1.\n"
                           "s: Time integration scheme. This is a string to choose among the values: forward-euler,\n"
                           "   backward-euler, crank-nicolson, heun, newmark, leap-frog, AB2 (Adams-Bashforth, 2nd Order),\n"
                           "   RK4 (Runge-Kutta, 4th Order), RK3-TVD (Runge-Kutta, 3rd order, TVD), BDF2 (Backward Difference\n"
                           "   Formula, 2nd Order), builtin (Any scheme built in the chosen PDE). The default value for this\n"
                           "   argument is backward-euler\n"
                           "adapted: Toggle meaning that adaptive time stepping is chosen.\n";
   static const vector<string> kw_scheme {"forward-euler","backward-euler","crank-nicolson","heun","newmark",
                                          "leap-frog","AB2","RK4","RK3-TVD","BDF2","builtin"};
   static const vector<string> kw {"help","?","set","initial$-time","final$-time","time$-step","adapted",
                                   "scheme","end","<","quit","exit","EXIT"};
   _cmd->set(kw);
   _nb_args = _cmd->getNbArgs();
   if (_nb_args==0) {
      msg("transient>","No argument for command.",H);
      _ret = 1;
      return;
   }
   if (_nb_args<1) {
      msg("transient>","No valid argument for command.");
      _ret = 1;
   }
   for (int k=0; k<_nb_args; ++k) {
      switch (_cmd->getArg()) {

         case 0:
         case 1:
            cout << H << endl;
            _ret = 0;

         case 3:
            it = _cmd->double_token();
            break;

         case 4:
            ft = _cmd->double_token();
            break;

         case 5:
            ts = _cmd->double_token();
            break;

         case 6:
            _adapted_time_step = _cmd->int_token();
            break;

         case 7:
            _scheme = _cmd->string_token();
            break;

         default:
            msg("transient>","Unknown argument: "+_cmd->Arg());
            _ret = 1;
            return;
      }
   }

   if (_nb_args>0) {
      _analysis_type = TRANSIENT;
      _time_step = ts;
      _init_time = it;
      _final_time = ft;
      if (ts<0)
         _time_step = -_time_step, _adapted_time_step = 1;
      if (find(kw_scheme.begin(),kw_scheme.end(),_scheme)==kw_scheme.end()) {
         msg("transient>","The scheme "+_scheme+" is unknown or unimplemented.");
         _ret = 1;
         return;
      }
      *ofh << "transient  initial-time=" << _init_time << "  final-time=" << _final_time
            << "  time-step=" << _time_step << "  adapted=" << _adapted_time_step
            << "  scheme=" << _scheme << endl;
   }
   else {
      if (_verb) {
         cout << "Default values:\n";
         cout << "   Initial time:            " << _init_time << endl;
         cout << "   Final time:              " << _final_time << endl;
         cout << "   Time step:               " << _time_step << endl;
         cout << "   Adapted time step ?      " << _adapted_time_step << endl;
         cout << "   Time integration scheme: " << _scheme << endl;
      }
      _cmd->setNbArg(0);
      *ofh << "transient" << endl;
      while (1) {
         int nb = _cmd->readline("rita>transient> ");
         if (nb<0)
            continue;
         switch (_key=_cmd->getKW(kw)) {

            case 0:
            case 1:
               cout << "\nAvailable Commands:\n";
               cout << "initial-time: Initial time value (Default = 0.)\n";
               cout << "final-time:   Final time value (Default = 1.)\n";
               cout << "time-step:    Time step (Default = 0.1, < 0: Adapted)\n";
               cout << "scheme:       Time integration scheme (Default = backward-euler)\n";
               cout << "              Available schemes: forward-euler, backward-euler, crank-nicolson, heun,\n";
               cout << "                                 newmark, leap-frog, AB2, RK4, RK3-TVD, BDF2, builtin\n";
               cout << "end or <:     back to higher level" << endl;
               break;

            case 2:
               setConfigure();
               return;

            case 3:
               if (_cmd->setNbArg(1,"Initial time to be given.")) {
                  msg("transient>initial-time>","Missing initial time value.","",1);
                  break;
               }
               ret = _cmd->get(it);
               if (!ret)
                  *ofh << "  initial-time " << it << endl;
               break;

            case 4:
               if (_cmd->setNbArg(1,"Final time to be given.")) {
                  msg("transient>final-time>","Missing final time value.","",1);
                  break;
               }
               ret = _cmd->get(ft);
               if (!ret)
                  *ofh << "  final-time " << ft << endl;
               break;

            case 5:
               if (_cmd->setNbArg(1,"Time step to be given.")) {
                  msg("transient>time-step>","Missing time step value.","",1);
                  break;
               }
               ret = _cmd->get(ts);
               if (!ret)
                  *ofh << "  time-step " << ts << endl;
               break;

            case 6:
               if (_cmd->setNbArg(1,"Time integration scheme.")) {
                  msg("transient>scheme>","Missing time integration scheme.","",1);
                  break;
               }
               ret = _cmd->get(kw_scheme,_scheme);
               if (ret<0) {
                  msg("transient>scheme>","Unknown time integration scheme.",
                      "Unknown time integration scheme.\n"
                      "Available values: forward-euler, backward-euler, crank-nicolson\n"
                      "                  heun, newmark, leap-frog, AB2, RK4, RK3-TVD, BDF2, builtin");
                  break;
               }
               *ofh << "  scheme " << _scheme << endl;
               break;

            case 7:
            case 8:
               *ofh << "  end" << endl;
               _analysis_type = TRANSIENT;
               _time_step = ts;
               _init_time = it;
               _final_time = ft;
               _adapted_time_step = 0;
               if (ts<0)
                  _time_step = -_time_step, _adapted_time_step = 1;
               _ret = 0;
               return;

            case -2:
            case -3:
            case -4:
               break;

            case 9:
            case 10:
               _ret = 100;
               return;

            case 11:
               _ret = 200;
               return;

            default:
               msg("transient>","Unknown command "+_cmd->token(),
                   "Available commands: initial-time, final-time, time-step, scheme, end, <\n"
                   "Global commands:    help, ?, set, quit, exit");
               break;
         }
      }
   }
}


void rita::setOptim()
{
   _optim->run();
}


void rita::setEigen()
{
   _analysis_type = EIGEN;
}


void rita::setAlgebraic()
{
   _nb_fields = _data->getNbFields();
   _ae = new odae;
   if (_analysis_type==NONE)
      _analysis_type = STEADY_STATE;
   _ret = runAE();
   if (_ret && _ae!=nullptr) {
      delete _ae;
      _ae = nullptr;
   }
   else {
      _ae->isSet = true;
      ALGEBRAIC.push_back(_ae);
      ODE.push_back(nullptr);
      PDE.push_back(nullptr);
      _ieq++, _nb_eq++;
      _eq_type.push_back(ALGEBRAIC_EQ);
   }
   if (_ret>=100)
      return;
}


void rita::setODE()
{
   _nb_fields = _data->getNbFields();
   _analysis_type = TRANSIENT;
   _ode = new odae;
   runODE();
   if (_ret && _ode!=nullptr) {
      delete _ode;
      _ode = nullptr;
   }
   else {
      ALGEBRAIC.push_back(nullptr);
      _ode->isSet = true;
      ODE.push_back(_ode);
      PDE.push_back(nullptr);
      _ieq++, _nb_eq++;
      _eq_type.push_back(ODE_EQ);
   }
   if (_ret>=100)
      return;
}


void rita::setPDE()
{
   string pde_name;
   int ret=0;
   _nb_fields = _data->getNbFields();
   if (_cmd->setNbArg(1,"Give PDE name.")) {
      msg("pde>","Missing pde name.","",1);
      _ret = 1;
      return;
   }
   if (_ieq==1) {
      msg("pde>","Only one equation is allowed.");
      return;
   }
   const vector<string> kw {"help","?","set","laplace","heat","wave","transport","linear-elasticity",
                            "truss","beam","incompressible-navier-stokes","clear","end","<","quit","exit","EXIT"};
   ret = _cmd->get(kw,pde_name);
   if (ret<0) {
      msg("pde>","Unknown pde "+pde_name,
          "Unknown pde name. Available pde's are:\n"
          "laplace, heat, wave, transport, linear-elasticity");
      _ret = 1;
      return;
   }
   _pde = new equa(this);
   _pde->set(pde_name,_theMesh);
   *ofh << "pde " << pde_name << endl;
   runPDE();
   if (_ret && _pde!=nullptr) {
      delete _pde; 
      _pde = nullptr;
   }
   else {
      _pde->set();
      PDE.push_back(_pde);
      ALGEBRAIC.push_back(nullptr);
      ODE.push_back(nullptr);
      _ieq++, _nb_eq++;
      _eq_type.push_back(PDE_EQ);
   }
   if (_ret>=100)
      return;
}


void rita::setApproximation()
{
   if (!_approx->run())
      _approx->go();
   _ret = 0;
   return;
}


void rita::setIntegration()
{
   if (!_integration->run())
      _integration->go();
   _ret = 0;
   return;
}


int rita::set_nls(string nls)
{
   auto it = NLs.find(nls);
   if (it==NLs.end()) {
      msg("equation>pde>nls>","Unknown nonlinear iterative solver: "+nls);
      return 1;
   }
   return 0;
}


int rita::set_ls(string ls,
                 string prec)
{
   auto it1 = Ls.find(ls);
   if (it1==Ls.end()) {
      msg("equation>pde>ls>","Unknown linear solver: "+ls);
      return 1;
   }
   auto it2 = Prec.find(prec);
   if (it2==Prec.end()) {
      msg("equation>pde>ls>","Unknown linear preconditioner: "+prec);
      return 1;
   }
   return 0;
}


int rita::setSpaceDiscretization(string& sp)
{
   _ret = 0;
   if (_cmd->setNbArg(1,"Missing space discretization method.")) {
      msg("equation>pde>space>","Missing space discretization method.","",1);
      return 1;
   }
   static const vector<string> kw {"fd","feP1","feP2","feQ1","fv"};
   _ret = _cmd->get(kw,sp);
   if (_ret<0) {
      msg("equation>pde>space>","Unknown space discretization method.",
          "Available Commands\n"
          "fd:   Finite Differences\n"
          "feP1: Finite elements, P1\n"
          "feP2: Finite elements, P2\n"
          "feQ1: Finite elements, Q1\n"
          "fv:   Finite volumes");
      return 1;
   }
   *ofh << "  space " << sp << endl;
   return 0;
}


void rita::setSummary()
{
   cout << "\nSUMMARY:" << endl;
   cout << "---------------------------------------------------------------" << endl;
   if (meshOK==false)
      msg("summary>","No mesh defined.");
   else {
      cout << "MESH DATA" << endl;
      cout << "Number of nodes:    " << _theMesh->getNbNodes() << endl;
      cout << "Number of elements: " << _theMesh->getNbElements() << endl;
      cout << "Number of sides:    " << _theMesh->getNbSides() << endl << endl;
   }
   cout << "MODEL ANALYSIS" << endl;
   if (_analysis_type == STEADY_STATE)
      cout << "Analysis is steady-state." << endl;
   else if (_analysis_type == TRANSIENT) {
      cout << "Analysis is transient." << endl;
      cout << "Initial time: " << _init_time << endl;
      cout << "Final time: " << _final_time << endl;
      cout << "Time step: " << _time_step << endl;
      cout << "Adapted time step ? " << _adapted_time_step << endl;
      cout << "Time integration scheme: " << _scheme << endl;
   }
   else if (_analysis_type == EIGEN) {
      cout << "Analysis is an eigenproblem analysis." << endl;
      cout << "Number of eigenvalues to extract: " << _nb_eigv << endl;
      cout << "Extract eigenvectors ? " << _eigen_vectors << endl;
   }
   else if (_analysis_type == OPTIMIZATION) {
      cout << "Analysis is optimization." << endl;
      //      cout << "Type of objective function: " << objective_type << endl;
      //      cout << "Optimization algorithm: " << opt_alg << endl;
   }
   for (int i=0; i<_nb_eq; ++i) {
      cout << "PDE: " << PDE[i]->eq << endl;
      //      cout << "Field: " << PDE[i]->field << endl;
      cout << "Space discretization method: " << PDE[i]->spD << endl;
      cout << "Linear system solver: " << rLs[PDE[i]->ls] << endl;
      cout << "Linear system preconditioner: " << rPrec[PDE[i]->prec] << endl;
   }
   cout << "---------------------------------------------------------------" << endl;
}


void rita::setClear()
{
}


void rita::setSolve()
{
   if (_verb>1)
      cout << "Entering mode 'solve' ..." << endl;
   _cmd->setNbArg(0);
   *ofh << _kw[_key] << endl;
   _solve->setVerbose(_verb);
   _solve->set(_theMesh);
   _ret = _solve->run();
   if (_ret==0)
      solveOK = true;
   if (_ret>=100)
      return;
   if (_verb>1)
      cout << "Leaving mode 'solve' ..." << endl;
}


void rita::msg(const string& loc, const string& m1, const string& m2, int c)
{
   if (c==0) {
   cout << "Error: " << m1 << endl;
   if (m2!="")
      cout << m2 << endl;
   }
   *ofl << "In rita>" + loc + ": " << m1 << endl;
}


odae::odae()
     : isSet(false), log(false), isFct(false), field(-1)
{
}


void odae::setVars(int opt)
{
   if (opt)
      vars.push_back("t");
   if (size==1)
      vars.push_back(fn);
   else {
      for (int i=0; i<size; ++i)
         vars.push_back(fn+to_string(i+1));
   }
}

} /* namespace RITA */
