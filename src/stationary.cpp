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

                      Implementation of class 'stationary'

  ==============================================================================*/

#include "stationary.h"
#include "io/IOField.h"
#include "io/saveField.h"
#include "equa.h"
#include <iostream>

using std::map;

namespace RITA {


stationary::stationary(rita *r)
           : _rita(r), _rs(1)
{
   _data = _rita->_data;
   _nb_fields = _rita->getNbFields();
   _nb_eq = _rita->getNbEq();
   _eq_type = &(_rita->_eq_type);
   _pde_eq.resize(_nb_eq);
   _pde_eq = _rita->PDE;
   _alg_eq.resize(_nb_eq);
   _alg_eq = _rita->ALGEBRAIC;
}


stationary::~stationary()
{
}


void stationary::setSave(vector<int>&    isave,
                         vector<int>&    fformat,
                         vector<string>& save_file)
{
   _isave = &isave;
   _fformat = &fformat;
   _save_file = &save_file;
}


int stationary::run()
{
   int ret = 0;
   vector<ofstream> fs(_nb_fields), ffs(_nb_fields), pfs(_nb_fields);
   try {
      vector<OFELI::IOField> ff(_nb_fields);
      vector<string> fn(_nb_fields);
      for (int e=0; e<_nb_eq; ++e) {
         if ((*_eq_type)[e]==ALGEBRAIC_EQ) {
            int f = _alg_eq[e]->field;
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
         if ((*_eq_type)[0]==ALGEBRAIC_EQ && _rs) {
            int f = _alg_eq[0]->field;
            _data->u[f]->resize(_alg_eq[0]->size);
            *_data->u[f] = _alg_eq[0]->y;
            fs[f].open(fn[f].c_str(),std::fstream::out);
            fs[f] << "# Saved by rita: Solution of Algebraic Equation, equation: 1" << endl;
            fs[f] << 0.;
            for (int i=0; i<_alg_eq[0]->size; ++i)
               fs[f] << "  " << _alg_eq[0]->y[i];
            fs[f] << endl;
            if ((*_isave)[0]) {
               ffs[f].open((*_save_file)[f].c_str());
               ffs[f] << "# Saved by rita: Solution of Algebraic Equation, equation: 1" << endl;
               ffs[f] << 0.;
               for (int i=0; i<_alg_eq[0]->size; ++i)
                  ffs[f] << "  " << _alg_eq[0]->y[i];
               ffs[f] << endl;
            }
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
            if ((*_eq_type)[e]==ALGEBRAIC_EQ) {
               int f = _alg_eq[e]->field;
               _data->u[f]->resize(_alg_eq[e]->size);
               *_data->u[f] = _alg_eq[e]->y;
               if ((*_isave)[f] && _rs) {
                  ffs[f].open((*_save_file)[e].c_str(),std::fstream::out);
                  ffs[f] << "  " << _alg_eq[e]->y;
                  if ((*_isave)[f]) {
                     ffs[f].open((*_save_file)[f].c_str(),std::fstream::out);
                     ffs[f] << "# Saved by rita: Solution of ODE, equation: " << e+1 << endl;
                     ffs[f] << 0.;
                     for (int i=0; i<_alg_eq[e]->size; ++i)
                        ffs[f] << "  " << _alg_eq[e]->y[i];
                  }
               }
	    }
            else if ((*_eq_type)[e]==PDE_EQ && _rs) {
               for (int i=0; i<_pde_eq[e]->nb_fields; ++i) {
                  int f = _pde_eq[e]->field[i];
                  ff[_pde_eq[e]->field[f]].open(fn[f],OFELI::IOField::OUT);
	       }
	    }
         }
      }

      for (int e=0; e<_nb_eq; ++e) {

//       Case of an algebraic equation
         if (_rita->_eq_type[e]==ALGEBRAIC_EQ) {
            NLASSolver nls(_alg_eq[e]->nls,_alg_eq[e]->size);
            if (_alg_eq[e]->size==1)
               nls.setInitial(_alg_eq[e]->y[0]);
            else
               nls.setInitial(_alg_eq[e]->y);
            for (int i=0; i<_alg_eq[e]->size; ++i)
               nls.setf(_alg_eq[e]->theFct[i]);
            nls.run();
            *_data->u[_alg_eq[e]->field] = _alg_eq[e]->y;
         }

//       Case of a PDE
         else if (_rita->_eq_type[e]==PDE_EQ) {
            equa *pde = _pde_eq[e];

//          Solution
            pde->theEquation->setInput(SOLUTION,*_data->u[pde->field[0]]);

//          Linear system solver
            pde->theEquation->setSolver(pde->ls,pde->prec);

//          Boundary condition
            if (pde->set_bc) {
               if (pde->bc.withRegex(1)) {
                  for (auto const& v: pde->regex_bc)
                     pde->bc.setNodeBC(v.first,v.second);
	       }
               pde->theEquation->setInput(BOUNDARY_CONDITION,pde->bc);
            }

//          Body force
            if (pde->set_bf) {
               if (pde->bf.withRegex(e+1))
                  pde->bf.set(pde->regex_bf);
               pde->theEquation->setInput(BODY_FORCE,pde->bf);
            }

//          Boundary force
            if (pde->set_sf) {
               if (pde->sf.withRegex(e+1)) {
                  for (auto const& v: pde->regex_sf)
                     pde->sf.setSideBC(v.first,v.second);
	       }
               pde->theEquation->setInput(BOUNDARY_FORCE,pde->sf);
            }

//          Run
            ret = pde->theEquation->run();

//          Save solution in file
            for (int i=0; i<pde->nb_fields; ++i) {
               int f = pde->field[i];
	    //            _data->u[f]->setName(_data->Field[f]);
               if (_rs) {
                  ff[f].open(fn[f],OFELI::IOField::OUT);
                  ff[f].put(*_data->u[f]);
                  ff[f].close();
                  if ((*_isave)[f])
                     OFELI::saveFormat(*(_rita->_theMesh),fn[f],(*_save_file)[f],(*_fformat)[f],
                                       (*_isave)[f]);
               }
            }
         }
      }
   } CATCH
   return ret;
}

} /* namespace RITA */
