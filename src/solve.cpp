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

                         Implementation of class 'solve'

  ==============================================================================*/

#include "rita.h"
#include "solve.h"
#include "configure.h"
#include "data.h"
#include "equa.h"
#include "stationary.h"
#include "transient.h"
#include "optim.h"
#include "util/macros.h"

namespace RITA {

solve::solve()
{
   _ret = 0;
   _key = 0;
}


solve::solve(rita *r, cmd *command, configure *config)
      : _rita(r), _set_analytic(false), _solved(false), _phase(false),
        _verb(1), _configure(config), _cmd(command)
{
   _ret = 0;
   _key = 0;
   _data = _rita->_data;
}

void solve::set(OFELI::Mesh *ms)
{
  _theMesh = ms;
}


int solve::run()
{
   _nb_eq = _rita->getNbEq();
   if ((_rita->_analysis_type==STEADY_STATE||_rita->_analysis_type==TRANSIENT) && _nb_eq==0) {
      _rita->msg("solve>","No equation defined.");
      _ret = 1;
      return _ret;
   }
   _nb_fields = _data->getNbFields();
   if (_nb_fields==0) {
      _rita->msg("solve>","No field defined.");
      _ret = 1;
      return _ret;
   }
   _fformat.resize(_nb_fields);
   _isave.resize(_nb_fields,0);
   _save_file.resize(_nb_fields);
   _phase_file.resize(_nb_fields);
   int eq=1, i=0;
   for (int e=0; e<_nb_eq; ++e) {
      if ((_rita->_eq_type)[e]==ALGEBRAIC_EQ)
         _rita->ALGEBRAIC[e]->analytic.resize(_rita->ALGEBRAIC[e]->size);
      else if ((_rita->_eq_type)[e]==ODE_EQ)
         _rita->ODE[e]->analytic.resize(_rita->ODE[e]->size);
      else if ((_rita->_eq_type)[e]==PDE_EQ) {
         if (_rita->PDE[e]->log.fail()) {
            if (_rita->PDE[e]->log.pde)
               _rita->msg("solve>","No PDE defined for equation "+to_string(e+1)+". Solution procedure aborted.");
            else if (_rita->PDE[e]->log.field)
               _rita->msg("solve>","Field definition incorrect for equation "+to_string(e+1)+". Solution procedure aborted.");
            else if (_rita->PDE[e]->log.spd)
               _rita->msg("solve>","No space discretization method available for PDE "+to_string(e+1)+". Solution procedure aborted.");
            else if (_rita->PDE[e]->log.ls)
               _rita->msg("solve>","No defined linear solver for PDE "+to_string(e+1)+". Solution procedure aborted.");
            else if (_rita->PDE[e]->log.nl)
               _rita->msg("solve>","No defined nonlinear solver for PDE "+to_string(e+1)+". Solution procedure aborted.");
            if (_verb>1)
               cout << "Getting back to higher level ..." << endl;
            _ret = 0;
            return _ret;
         }
         else
            _rita->PDE[e]->analytic.resize(1);
      }
   }

   while (1) {
      int nb = _cmd->readline("rita>solve> ");
      if (nb<0)
         continue;
      switch (_key=_cmd->getKW(_kw_solve)) {

         case 0:
         case 1:
            _cmd->setNbArg(0);
            cout << "\nAvailable Commands:\n";
            cout << "run:      Run the model\n";
            cout << "save:     Save output, can be executed before run\n";
            cout << "display:  Display solution and related data\n";
            cout << "plot:     Plot solution\n";
            cout << "analytic: Give analytic solution to test accuracy\n";
            cout << "error:    Compute error in various norms\n";
            cout << "post:     Post calculations\n";
            cout << "end or <: go back to higher level\n" << endl;
            cout << "Global commands:\n";
            cout << "help, ?, set, quit, exit" << endl;
            break;

         case 2:
            if (_verb)
               cout << "Setting configuration parameter(s) ..." << endl;
            _configure->setVerbose(_verb);
            _configure->run();
            _verb = _configure->getVerbose();
            break;

         case 3:
            if (_rita->_analysis_type==STEADY_STATE) {
               if (_verb)
                  cout << "Running stationary solver ..." << endl;
               *_rita->ofh << "  run" << endl;
               _ret = run_steady();
            }
            else if (_rita->_analysis_type==TRANSIENT) {
               if (_verb)
                  cout << "Running transient solver ..." << endl;
               *_rita->ofh << "  run" << endl;
               _ret = run_transient();
            }
            else if (_rita->_analysis_type==OPTIMIZATION) {
               if (_verb)
                  cout << "Running optimization problem solver ..." << endl;
               *_rita->ofh << "  run" << endl;
               _ret = run_optim();
            }
            else if (_rita->_analysis_type==EIGEN) {
               if (_verb)
                  cout << "Running eigen problem solver ..." << endl;
               *_rita->ofh << "  run" << endl;
               _ret = run_eigen();
            }
            else
               ;
            break;

         case 4:
            save();
            break;

         case 5:
            *_rita->ofh << "  display " << endl;
            display();
            break;

         case 6:
            _ret = plot();
            break;

         case 7:
            if (_verb)
               cout << "Setting analytical solution ..." << endl;
            *_rita->ofh << "  analytic" << endl;
            setAnalytic();
            _set_analytic = true;
            break;

         case 8:
            if (_verb)
               cout << "Computing error ..." << endl;
            eq = 1, i = 0;
            if (nb>1)
               _ret = _cmd->get(eq);
            if (nb>2)
               _ret = _cmd->get(i);
            *_rita->ofh << "  error  " << eq << "  " << i << endl;
            get_error(eq,i);
            break;

         case 9:
            if (_verb)
               cout << "Setting post calculations ..." << endl;
            *_rita->ofh << "  post" << endl;
            break;

         case 10:
         case 11:
            if (_verb>1)
               cout << "Getting back to higher level ..." << endl;
            *_rita->ofh << "  end" << endl;
            _ret = 0;
            return _ret;

         case 12:
         case 13:
            _ret = 100;
            return _ret;

         case 14:
            _ret = 200;
            return _ret;

         case -2:
            break;

         default:
            _rita->msg("solve>","Unknown command: "+_cmd->token(),
                       "Available commands for this mode:\n"
                       "help, ?, run, save, display, plot, analytic, error, post, end, <, exit");
            break;
      }
   }
}


int solve::run_steady()
{
   stationary st(_rita);
   st.setSave(_isave,_fformat,_save_file);
   int ret = 0;
   try {
      ret = st.run();
      _solved = true;
   } CATCH
   return ret;
}


int solve::run_transient()
{
   transient ts(_rita);
   ts.setSave(_isave,_fformat,_save_file,_phase,_phase_file);
   for (int e=0; e<_nb_eq; ++e) {
      if ((_rita->_eq_type)[e]==ALGEBRAIC_EQ) {
      }
      else if ((_rita->_eq_type)[e]==ODE_EQ) {
      }
      else if ((_rita->_eq_type)[e]==PDE_EQ)
         ts.setLinearSolver(_rita->PDE[e]->ls,_rita->PDE[0]->prec);
   }
   int ret = ts.run();
   for (int e=0; e<_nb_eq; ++e) {
      if (_isave[e] && theStep%_isave[e]==0) {
         if ((_rita->_eq_type)[e]==ALGEBRAIC_EQ) {
         }
         else if ((_rita->_eq_type)[e]==ODE_EQ) {
         }
         else if ((_rita->_eq_type)[e]==PDE_EQ)
           ;
      }
   }
   _solved = true;
   return ret;
}


int solve::run_eigen()
{
  //   return _rita->PDE[0].Eq->run();
   _rita->msg("solve>run>","Eigenproblem solver not implemented.");
   return 0;
}


int solve::run_optim()
{
   _optim = _rita->_optim;
   if (_optim->log) {
      _rita->msg("solve>run_optim>","Optimization problem undefined or improperly defined.");
      return 1;
   }
   int size=_optim->size;
   try {
      if (_optim->lp) {
         OFELI::LPSolver s;
         s.setSize(size,_optim->nb_lec,_optim->nb_gec,_optim->nb_eqc);
         s.set(*_data->u[_optim->ifield]);
         s.set(OFELI::LPSolver::OBJECTIVE,_optim->a,_optim->b);
         for (int i=0; i<_optim->nb_eqc; ++i)
            s.set(OFELI::LPSolver::EQ_CONSTRAINT,*_optim->a_eq[i],_optim->b_eq[i]);
         for (int i=0; i<_optim->nb_lec; ++i)
            s.set(OFELI::LPSolver::LE_CONSTRAINT,*_optim->a_le[i],_optim->b_le[i]);
         for (int i=0; i<_optim->nb_gec; ++i)
            s.set(OFELI::LPSolver::GE_CONSTRAINT,*_optim->a_ge[i],_optim->b_ge[i]);
         int ret = s.run();
         _optim->obj = s.getObjective();
	 if (ret==0)
            _optim->solved = true;
      }
      else {
         OFELI::OptSolver s(*_data->u[_optim->ifield]);
         s.setOptMethod(_optim->Alg);
         s.setObjective(*_optim->J_Fct);
         if (_optim->G_ok) {
            for (int i=0; i<size; ++i)
               s.setGradient(*_data->theFct[_optim->igrad+i],i+1);
         }
         if (_optim->H_ok) {
            for (int i=0; i<size*size; ++i)
               s.setHessian(*_data->theFct[_optim->ihess+i],i+1);
         }
         for (int i=0; i<_optim->nb_lec; ++i)
            s.setIneqConstraint(*_optim->inC_Fct[i],_optim->penal);
         for (int i=0; i<_optim->nb_eqc; ++i)
            s.setEqConstraint(*_optim->eqC_Fct[i],_optim->penal);
         s.setLowerBounds(_optim->lb);
         s.setUpperBounds(_optim->ub);
         if (!s.run()) {
            _optim->obj = s.getObjective();
            _optim->solved = true;
         }
      }
   } CATCH
   return 0;
}


void solve::save()
{
   int k=0, freq=1, eq=0, field_ok=0;
   string file="rita-1.pos", fformat="gmsh", ext="", fd="u", phase_f="";
   _phase = false;
   if (_nb_fields==1) {
      fd = _data->Field[0];
      eq = _data->FieldEquation[0];
      _fformat[0] = GNUPLOT;
      _save_file[0] = "u.dat";
      _isave[0] = 1;
      field_ok = 1;
      if ((_rita->_eq_type)[0]==PDE_EQ) {
         _fformat[0] = GMSH;
         _save_file[0] = "u.pos";
      }
   }
   _ret = 0;

   vector<string> kw = {"field","format","freq$uency","phase","file"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg("=");
      switch (n) {

         case  0:
            fd = _cmd->string_token();
            field_ok = 2;
            break;

         case 1:
            fformat = _cmd->string_token();
            break;

         case 2:
            freq = _cmd->int_token();
            break;

         case 3:
            phase_f = _cmd->string_token();
            _phase = true;
            break;

         case 4:
            file = _cmd->string_token();
            break;

         default:
            _rita->msg("solve>save>","Unknown argument: "+kw[n]);
            return;
      }
   }
   if (nb_args>0) {
      if (_data->getNbFields()==0) {
         _rita->msg("solve>","No field to save.");
         _ret = 1;
         return;
      }
      if (!field_ok) {
         _rita->msg("solve>:","No field name given.");
         _ret = 1;
         return;
      }
      k = _data->checkField(fd);
      if (k<0) {
         _rita->msg("solve>:","Unknown field: "+fd);
         _ret = 1;
         return;
      }
      eq = _data->FieldEquation[k];
      if ((_rita->_eq_type)[eq]==PDE_EQ) {
         _isave[k] = freq;
         _fformat[k] = GMSH;
         _save_file[k] = fd + ".pos";
      }
      else if ((_rita->_eq_type)[eq]==ODE_EQ ||
               (_rita->_eq_type)[eq]==ALGEBRAIC_EQ) {
         _isave[k] = freq;
         _fformat[k] = GNUPLOT;
         _save_file[k] = fd + ".dat";
      }
      eq = _data->FieldEquation[k];
      if ((_rita->_eq_type)[eq]!=PDE_EQ && fformat != "gnuplot") {
         _rita->msg("solve>save>","Only Gnuplot format is available for this type of equation.");
         _ret = 1;
         return;
      }
      if (!_ret)
         _fformat[k] = _ff[fformat];
      else {
         _rita->msg("solve>save>","Unknown file format: "+fformat);
         _ret = 1;
         return;
      }
      _save_file[k] = file;
      *_rita->ofh << "  save  field=" << fd << "  file=" << file << "  format=" << fformat
            << "  frequency=" << freq;
      if (_phase) {
         _phase_file[k] = phase_f;
         *_rita->ofh << "  phase=" << phase_f;
      }
      *_rita->ofh << endl;
   }
   _ret = 0;
   return;
}


void solve::setAnalytic()
{
   int eq=1, comp=1;
   bool exp_ok=false;
   string exp="";
   _ret = 0;

   vector<string> kw = {"eq$uation","comp$onent","def$inition"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg();
      switch (n) {

         case  0:
            eq = _cmd->int_token();
            break;

         case 1:
            comp = _cmd->int_token();
            break;

         case 2:
            exp = _cmd->string_token();
            exp_ok = true;
            break;

         default:
            _rita->msg("solve>analytic>","Unknown argument: "+kw[n]);
            _ret = 1;
            return;
      }
   }
   if (nb_args>0) {
      if (!exp_ok) {
         _rita->msg("solve>analytic>","No expression of analytical solution provided.");
         _ret = 1;
         return;
      }
      if ((_rita->_eq_type)[eq-1]==ALGEBRAIC_EQ)
         _rita->ALGEBRAIC[eq-1]->analytic[comp-1] = exp;
      else if ((_rita->_eq_type)[eq-1]==ODE_EQ)
         _rita->ODE[eq-1]->analytic[comp-1] = exp;
      else if ((_rita->_eq_type)[eq-1]==PDE_EQ)
         _rita->PDE[eq-1]->analytic[comp-1] = exp;
      *_rita->ofh << "analytic equation=" << eq << " component=" << comp << " expression=" << exp << endl;
   }
   _ret = 0;
   return;
}


void solve::display()
{
   for (int f=0; f<_nb_fields; ++f) {
      if (_data->FieldType[f]==ALGEBRAIC_EQ) {
         if (_data->u[f]->size()==1)
            cout << "Solution of algebraic equation: " << (*_data->u[f])[0] << endl;
         else {
            cout << "Solution of algebraic equation: " << endl;
            for (size_t i=0; i<_data->u[f]->size(); ++i)
               cout << (*_data->u[f])[i] << endl;
         }
      }
      else if (_data->FieldType[f]==ODE_EQ) {
         if (_data->u[f]->size()==1)
            cout << "Solution of ordinary differential equation: "
                 << (*_data->u[f])[0] << endl;
         else {
            cout << "Solution of ordinary differential equation: " << endl;
            for (size_t i=0; i<_data->u[f]->size(); ++i)
               cout << (*_data->u[f])[i] << endl;
         }
      }
      else if (_data->FieldType[f]==PDE_EQ) {
         cout << "Solution of partial differential equation: " << endl;
         cout << "Solution vector:\n" << *_data->u[f];
      }
      else if (_data->FieldType[f]==OPT) {
         if (!_optim->solved) {
            _rita->msg("solve>display>","Optimization problem ot properly solved.");
            _ret = 1;
            return;
         }
         if (_data->u[f]->size()==1)
            cout << "Solution of optimization problem: " << (*_data->u[f])[0] << endl;
         else {
            cout << "Solution of optimization problem: " << endl;
            for (size_t i=0; i<_data->u[f]->size(); ++i)
               cout << (*_data->u[f])[i] << endl;
         }
         cout << "Optimal objective: " << _rita->_optim->obj << endl;
      }
   }
}


int solve::plot()
{
   OFELI::saveField(*_data->u[0],"rita.pos",GMSH);
   string com = "gmsh rita.pos";
   int ret = system(com.c_str());
   remove("rita.pos");
   return ret;
}


void solve::get_error(int eq, int i)
{
   double err2=0., errI=0.;
   vector<double> y;
   if (!_set_analytic) {
      _rita->msg("solve>error>","No analytical solution given.");
      return;
   }
   if (!_solved) {
      _rita->msg("solve>error>","Problem has not been solved yet.");
      return;
   }

// Case of an algebraic equation
   if ((_rita->_eq_type)[eq-1]==ALGEBRAIC_EQ) {
      _var.resize(_rita->ALGEBRAIC[eq-1]->size);
      for (int i=0; i<_rita->ALGEBRAIC[eq-1]->size; ++i) {
         y.resize(_rita->ALGEBRAIC[eq-1]->size);
         for (int j=0; j<_rita->ALGEBRAIC[eq-1]->size; ++j) {
            _var[j] = _rita->ALGEBRAIC[eq-1]->fn+to_string(j+1);
            y[j] = _rita->ALGEBRAIC[eq-1]->y[j];
         }
         _theFct.set(_rita->ALGEBRAIC[eq-1]->analytic[i],_var);
         double u=_rita->ALGEBRAIC[eq-1]->y[i], v=_theFct(y);
         errI += (u-v)*(u-v);
      }
      cout << "Error: " << sqrt(errI) << endl;
   }

// Case of an ODE
   else if ((_rita->_eq_type)[eq-1]==ODE_EQ) {
       y.resize(_rita->ODE[eq-1]->size+1);
      _var.resize(_rita->ODE[eq-1]->size+1);
      _var[0] = "t";
      y[0] = _rita->_final_time;
      for (int i=0; i<_rita->ODE[eq-1]->size; ++i) {
         for (int j=0; j<_rita->ODE[eq-1]->size; ++j) {
            _var[j+1] = _rita->ODE[eq-1]->fn+to_string(j+1);
            y[j+1] = _rita->ODE[eq-1]->y[j];
         }
         _theFct.set(_rita->ODE[eq-1]->analytic[i],_var);
         double u=_rita->ODE[eq-1]->y[i], v=_theFct(y);
         errI += (u-v)*(u-v);
      }
      cout << "Error: " << sqrt(errI) << endl;
   }

// Case of a PDE
   else if ((_rita->_eq_type)[eq-1]==PDE_EQ) {
      int nb_dof = _theMesh->getNbDOF()/_theMesh->getNbNodes();
      for (int i=0; i<nb_dof; ++i) {
         _theFct.set(_rita->PDE[eq-1]->analytic[i]);
         for (_theMesh->topNode(); (OFELI::theNode=_theMesh->getNode());) {
            double u = (*_data->u[eq-1])(OFELI::theNodeLabel,i+1);
            double v = _theFct(TheNode.getCoord());
            err2 += (u-v)*(u-v);
            errI  = std::max(fabs(u-v),errI);
         }
      }
      err2 = sqrt(err2/_theMesh->getNbNodes());
      cout << "Error in L2-Norm:         " << err2 << endl;
      cout << "Error in L-Infinity-Norm: " << errI << endl;
   }

// Case of an Optimization problem
   else if ((_rita->_eq_type)[eq-1]==OPT)
      _rita->msg("solve>error>","Error computation is not available for optimization problems.");
}

} /* namespace RITA */
