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

                      Implementation of class 'transient'

  ==============================================================================*/

#include "OFELI_Config.h"
#include "util/util.h"
#include "io/IOField.h"
#include "io/saveField.h"
#include <iostream>

#include "transient.h"
#include "equa.h"

namespace RITA {

transient::transient(rita *r)
          : _ts_allocated(false), _ode_allocated(false), _nlas_allocated(false),
	    _phase(false), _rita(r), _rs(1), _nlas(nullptr), _ode(nullptr), _ts(nullptr)
{
   _data = _rita->_data;
   _nb_fields = _data->getNbFields();
   _nb_eq = _rita->getNbEq();
   _eq_type = &(_rita->_eq_type);
   theTimeStep = _time_step = _rita->_time_step;
   _init_time = _rita->_init_time;
   theFinalTime = _final_time = _rita->_final_time;
   _algebraic_eq = _rita->ALGEBRAIC;
   _ode_eq = _rita->ODE;
   for (int e=0; e<_nb_eq; ++e) {
      if ((*_eq_type)[e]==ALGEBRAIC_EQ) {
         _nlas = new OFELI::NLASSolver(_rita->ALGEBRAIC[e]->nls,_ode_eq[e]->size);
         _nlas_allocated = true;
      }
      else if ((*_eq_type)[e]==ODE_EQ) {
         _ode = new OFELI::ODESolver(_ode_eq[e]->scheme,_time_step,_final_time,_ode_eq[e]->size);
         _ode_allocated = true;
         for (int i=0; i<_ode_eq[e]->size; ++i)
            _ode->setInitial(_ode_eq[e]->y[i],i+1);
      }
      else if ((*_eq_type)[e]==PDE_EQ)
         setPDE(e);
   }
}


transient::~transient()
{
   if (_ts_allocated)
      delete _ts;
   if (_ode_allocated)
      delete _ode;
   if (_nlas_allocated)
      delete _nlas;
}


int transient::setPDE(int e)
{
   _pde_eq = _rita->PDE;
   try {
      _ts = new OFELI::TimeStepping(_rita->_sch[_rita->_scheme],_time_step,_final_time);
      _ts_allocated = true;
      _ts->setPDE(*(_pde_eq[e]->theEquation));
      _ts->setLinearSolver(_pde_eq[e]->ls,_pde_eq[e]->prec);
      _ts->setInitial(*_data->u[_pde_eq[e]->field[0]]);
      if (_pde_eq[e]->eq=="incompressible-navier-stokes" && _pde_eq[e]->spD=="feP1")
         _pde_eq[e]->theEquation->setInput(PRESSURE_FIELD,*_data->u[_pde_eq[e]->field[1]]);
   } CATCH
   return 0;
}


void transient::setLinearSolver(OFELI::Iteration      ls,
                                OFELI::Preconditioner prec)
{
  //   _ls = ls;
  //   _prec = prec;
}


void transient::setSave(vector<int>&    isave,
                        vector<int>&    fformat,
                        vector<string>& save_file,
			bool            phase,
                        vector<string>& phase_file)
{
   _isave = &isave;
   _fformat = &fformat;
   _save_file = &save_file;
   _phase = phase;
   _phase_file = &phase_file;
}


int transient::run()
{
   OFELI::Verbosity = 1;
   vector<ofstream> fs(_nb_fields), ffs(_nb_fields), pfs(_nb_fields);
   vector<OFELI::IOField> ff(_nb_fields);
   vector<string> fn(_nb_fields);
   for (int e=0; e<_nb_eq; ++e) {
      if ((*_eq_type)[e]==ODE_EQ) {
         int f = _ode_eq[e]->field;
         fn[f] = "rita-" + to_string(10*(e+1)+f+1) + ".sol";
      }
      else if ((*_eq_type)[e]==ALGEBRAIC_EQ) {
         int f = _algebraic_eq[e]->field;
         fn[f] = "rita-" + to_string(10*(e+1)+f+1) + ".sol";
      }
      else if ((*_eq_type)[e]==PDE_EQ) {
         for (int i=0; i<_pde_eq[e]->nb_fields; ++i) {
            int f = _pde_eq[e]->field[i];
            fn[f] = "rita-" + to_string(10*(e+1)+f+1) + ".sol";
         }
      }
   }
   if (_nb_eq==1) {
      if ((*_eq_type)[0]==ODE_EQ) {
         int f = _ode_eq[0]->field;
         _data->u[f]->resize(_ode_eq[0]->size);
         *_data->u[f] = _ode_eq[0]->y;
         if (_rs) {
            fs[f].open(fn[f].c_str(),std::fstream::out);
            fs[f] << "# Saved by rita: Solution of ODE, equation: 1" << endl;
            fs[f] << 0.;
            for (int i=0; i<_rita->_ode[0].size; ++i)
               fs[f] << "  " << _ode_eq[0]->y[i];
	    fs[f] << endl;
         }
         if ((*_isave)[0]) {
            ffs[f].open((*_save_file)[f].c_str());
            ffs[f] << "# Saved by rita: Solution of ODE, equation: 1" << endl;
            ffs[f] << 0.;
            for (int i=0; i<_rita->_ode[0].size; ++i)
               ffs[f] << "  " << _ode_eq[0]->y[i];
            ffs[f] << endl;
	    if (_phase) {
               pfs[f].open((*_phase_file)[f].c_str());
               pfs[f] << "# Saved by rita: Phase portrait of ODE, equation: 1" << endl;
            }
         }
         if (_ode_eq[0]->size==1)
            _ode->setInitial(_ode_eq[0]->y[0]);
         else
            _ode->setInitial(_ode_eq[0]->y);
      }
      else if ((*_eq_type)[0]==PDE_EQ && _rs) {
         for (int i=0; i<_pde_eq[0]->nb_fields; ++i) {
            int f = _pde_eq[0]->field[i];
            ff[f].open(fn[f],OFELI::IOField::OUT);
	 }
      }
   }
   else {
      for (int e=0; e<_nb_eq; ++e) {
         int f = _ode_eq[e]->field;
         _data->u[f]->resize(_ode_eq[e]->size);
         *_data->u[f] = _ode_eq[e]->y;
         if ((*_isave)[e]) {
            ffs[f].open((*_save_file)[e].c_str(),std::fstream::out);
            ffs[f] << "  " << _ode_eq[e]->y;
            if ((*_isave)[f]) {
               ffs[f].open((*_save_file)[f].c_str(),std::fstream::out);
               ffs[f] << "# Saved by rita: Solution of ODE, equation: " << e+1 << endl;
               ffs[f] << 0.;
               for (int i=0; i<_rita->_ode[e].size; ++i)
                  ffs[f] << "  " << _ode_eq[e]->y[i];
               ffs[f] << endl;
               if (_phase) {
                  pfs[f].open((*_phase_file)[f].c_str());
                  pfs[f] << "# Saved by rita: Phase portrait of ODE, equation: " << e+1 << endl;
               }
            }
         }
         if ((*_eq_type)[e]==ALGEBRAIC_EQ) {
            if (_algebraic_eq[e]->size==1)
               _nlas->setInitial(_ode_eq[e]->y[0]);
            else
               _nlas->setInitial(_ode_eq[e]->y);
         }
         if ((*_eq_type)[e]==ODE_EQ) {
            if (_ode_eq[e]->size==1)
               _ode->setInitial(_ode_eq[e]->y[0]);
            else
               _ode->setInitial(_ode_eq[e]->y);
         }
         else if ((*_eq_type)[e]==PDE_EQ && _rs) {
            for (int i=0; i<_pde_eq[e]->nb_fields; ++i) {
               int f = _pde_eq[e]->field[i];
               ff[_pde_eq[e]->field[f]].open(fn[f],OFELI::IOField::OUT);
            }
	 }
      }
   }
   theStep = 1;
   for (int e=0; e<_nb_eq; ++e) {
      if ((*_eq_type)[e]==ALGEBRAIC_EQ) {
         for (int i=0; i<_algebraic_eq[e]->size; ++i)
            _nlas->setf(_algebraic_eq[e]->theFct[i]);
      }
      else if ((*_eq_type)[e]==ODE_EQ) {
         for (int i=0; i<_ode_eq[e]->size; ++i)
            _ode->setF(_ode_eq[e]->theFct[i]);
      }
      else if ((*_eq_type)[e]==PDE_EQ)
         _ts->setLinearSolver(_pde_eq[e]->ls,_pde_eq[e]->prec);
   }

// Loop on time steps
   try {
      TimeLoop {

         if (_rita->_verb)
            cout << "Performing time step " << theStep <<", Time = " << theTime << endl;

//       Loop on equations to solve
         for (int e=0; e<_nb_eq; ++e) {

//          Case of an algebraic equation
            if ((*_eq_type)[e]==ALGEBRAIC_EQ) {
//               int f = _algebraic_eq[e]->field;
               cout << "No algebraic equation solver implemented." << endl;
            }

//          Case of an ODE
            else if ((*_eq_type)[e]==ODE_EQ) {
               int f = _ode_eq[e]->field;
               OFELI::Vect<double> z(_rita->_ode[e].size);
               _ode->runOneTimeStep();
               if (_rita->_ode[e].size==1)
                  _ode_eq[e]->y[0] = _ode->get();
               *_data->u[f] = _ode_eq[e]->y;
               if (_rs) {
                  fs[f] << theTime;
                  for (int i=0; i<_rita->_ode[e].size; ++i)
                     fs[f] << "  " << _ode_eq[e]->y[i];
                  fs[f] << endl;
               }
               if ((*_isave)[f]) {
                  if (theStep%(*_isave)[f]==0) {
                     ffs[f] << theTime;
                     for (int i=0; i<_rita->_ode[e].size; ++i)
                        ffs[f] << "  " << _ode_eq[e]->y[i];
                     ffs[f] << endl;
                     if (_phase) {
                        _ode->getTimeDerivative(z);
                        pfs[f] << _ode_eq[e]->y[0] << "  ";
                        for (int i=0; i<_rita->_ode[e].size; ++i)
                           pfs[f] << _ode->getTimeDerivative(i+1) << "  ";
                        pfs[f] << endl;
                     }
                  }
               }
            }

//          Case of a PDE
            else if ((*_eq_type)[e]==PDE_EQ) {
               equa *pde = _pde_eq[e];

//             Body force
               if (pde->set_bf) {
                  pde->bf.setTime(theTime);
                  if (pde->bf.withRegex(1))
                     pde->bf.set(pde->regex_bf);
                  _ts->setRHS(pde->bf);
                  pde->theEquation->setInput(BODY_FORCE,pde->bf);
               }

//             Boundary condition
               if (pde->set_bc) {
                  pde->bc.setTime(theTime);
                  if (pde->bc.withRegex(1)) {
                     for (auto const& v: pde->regex_bc)
                        pde->setNodeBC(v.first,v.second,theTime,pde->bc);
                  }
                  _ts->setBC(pde->bc);
               }

//             Boundary force
               if (pde->set_sf) {
                  pde->sf.setTime(theTime);
                  if (pde->sf.withRegex(1)) {
                     for (auto const& v: pde->regex_sf)
                        pde->setNodeBC(v.first,v.second,theTime,pde->sf);
                  }
            //            _ts->setSF(_data->sf[i]);
               }

//             Run
               _ts->runOneTimeStep();

//             Save in native OFELI format file
               for (int i=0; i<pde->nb_fields; ++i) {
                  int f = pde->field[i];
                  _data->u[f]->setTime(theTime);
                  _data->u[f]->setName(_data->Field[f]);
                  if (_rs)
                     ff[f].put(*_data->u[f]);
	       }
            }
         }
      }
   } CATCH

   theTime -= theTimeStep;
   for (int e=0; e<_nb_eq; ++e) {
      if ((*_eq_type)[e]==PDE_EQ) {
         for (int i=0; i<_pde_eq[e]->nb_fields; ++i) {
            int f = _pde_eq[e]->field[i];
            if (_rs) {
               ff[f].close();
               OFELI::saveFormat(*(_rita->_theMesh),fn[f],(*_save_file)[f],(*_fformat)[f],
                                 (*_isave)[f]);
            }
         }
      }
   }
   return 0;
}

} /* namespace RITA */
