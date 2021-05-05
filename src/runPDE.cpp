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

                       Implementation of function runPDE

  ==============================================================================*/

#include "rita.h"
#include "data.h"
#include "equa.h"
#include "cmd.h"
#include "configure.h"

using std::cout;
using std::exception;

namespace RITA {

void rita::runPDE()
{
   int nb_args = 0, nb_fields = 0, nb=0;
   bool field_ok = false;
   vector<string> field_name;
   string str = "", str1 = "";
   _pde->log.field = true;
   _ret = 0;
   if (_analysis_type==NONE)
      _analysis_type = STEADY_STATE;
   if (_theMesh==nullptr) {
      cout << "Error: No valid mesh created." << endl;
      *_ofl << "In rita>pde>: No mesh created";
      _pde->log.mesh = true;
      _ret = 1;
      return;
   }
   else if (_theMesh->getNbNodes()==0) {
      cout << "Error: Empty mesh." << endl;
      *_ofl << "In rita>pde>: Empty mesh";
      _pde->log.mesh = true;
      _ret = 1;
      return;
   }
   _pde->ls = OFELI::CG_SOLVER;
   _pde->prec = OFELI::DILU_PREC;
   _pde->spD = "feP1";
   const vector<string> kw = {"help","?","set","field","coef","in$it","bc","bf","source","sf","traction","space",
                              "ls","nls","clear","end","<","quit","exit","EXIT"};

   while (1) {
      if ((nb_args=_cmd->readline("rita>pde> "))<0)
         continue;
      switch (_key=_cmd->getKW(kw)) {

         case 0:
         case 1:
            cout << "\nAvailable Commands:\n";
            cout << "field:    Field name of an unknown of the equation\n";
            cout << "coef:     PDE coefficients\n";
            cout << "init:     Set initial condition or guess for pde\n";
            cout << "bc:       Set boundary conditions\n";
            cout << "source:   Set sources or body forces\n";
            cout << "sf:       Set side (boundary) forces\n";
            cout << "space:    Space discretization method\n";
            cout << "ls:       Set linear system solver\n";
            cout << "nls:      Set nonlinear system iteration procedure\n";
            cout << "clear:    Remove pde from model\n";
            cout << "end or <: go back to higher level" << endl;
            break;

         case 2:
            setConfigure();
            break;

         case 3:
            if (nb_fields==_pde->nb_fields) {
               cout << "Number of fields is larger than PDE requires it." << endl;
               *_ofl << "In rita>pde>field>: Too many fields for PDE." << endl;
               break;
            }
            if (_cmd->setNbArg(1,"Give name of an associated field.")) {
               *_ofl << "In rita>pde>field>: Missing name of an associated field." << endl;
               break;
            }
            if (!_cmd->get(str)) {
               field_ok = true, nb_fields++;
               field_name.push_back(str);
               *_ofh << "  field " << str << endl;
            }
            break;

         case 4:
            _pde->setCoef();
            break;

         case 5:
            _ret = _pde->setIn();
            break;

         case 6:
            _ret = _pde->setBC();
            break;

         case 7:
         case 8:
            _ret = _pde->setBF();
            break;

         case  9:
         case 10:
            _ret = _pde->setSF();
            break;

         case 11:
            if (_cmd->setNbArg(1,"Give space discretization method.")) {
               *_ofl << "In rita>pde>space>: Missing space discretization method." << endl;
               break;
            }
            _ret = setSpaceDiscretization(str);
            if (!_ret)
               _pde->spD = str;
	    else
               _pde->log.spd = true;
            break;

         case 12:
            if (_cmd->setNbArg(1,"Linear solver and optional preconditioner to be supplied.",1)) {
               *_ofl << "In rita>pde>ls>: Missing linear solver data." << endl;
               break;
            }
            nb = _cmd->getNbArgs();
            if (nb==0) {
               cout << "Missing linear system solver." << endl;
               *_ofl << "In rita>pde>ls>: Missing linear solver data." << endl;
            }
            _ret = _cmd->get(str);
            str1 = "ident";
            if (nb>1)
               _ret += _cmd->get(str1);
            if (!_ret) {
               *_ofh << "  ls " << str << " " << str1 << endl;
               if (!set_ls(str,str1)) {
                  _pde->ls = Ls[str];
                  _pde->prec = Prec[str1];
               }
            }
            else
               _pde->log.ls = true;
            break;

         case 13:
            if (_cmd->setNbArg(1,"Nonlinear solver to be supplied.",1)) {
               *_ofl << "In rita>pde>nls>: Missing nonlinear solver data." << endl;
               break;
            }
            _ret = _cmd->get(str);
            if (!_ret) {
               *_ofh << "  nls " << str << endl;
               _ret = set_nls(str);
               if (!_ret)
                  _pde->nls = str;
               else
                  _pde->log.nl = true;
            }
            break;

         case 14:
            cout << "PDE removed from model." << endl;
            *_ofh << "  clear" << endl;
            _ret = 10;
            return;

         case 15:
         case 16:
            _cmd->setNbArg(0);
	    if (!_ret) {
               if (!field_ok) {
                  cout << "Error: No field(s) defined for PDE." << endl;
                  *_ofl << "In rita>pde>end>: No field(s) defined for PDE." << endl;
                  break;
               }
               if (nb_fields!=_pde->nb_fields) {
                  cout << "Error: PDE requires exactly " << _pde->nb_fields << " field(s)." << endl;
                  *_ofl << "In rita>pde>end>: Incorrect number of fields for defined PDE." << endl;
                  break;
               }
	       for (int i=0; i<nb_fields; ++i) {
                  _ifield = _data->addField(field_name[i],NODES,0,_pde->nb_dof[i]);
                  _data->FieldEquation[_ifield] = _nb_eq;
                  _data->FieldType[_ifield] = PDE_EQ;
                  _pde->fn[i] = field_name[i];
                  _pde->field[i] = _ifield;
               }
               _nb_fields = _data->getNbFields();
               _pde->log.field = false;
               _pde->b.setSize(_theMesh->getNbEq());
               if (_verb) {
                  cout << "Summary of PDE attributes:" << endl;
                  cout << "   PDE name: " << _pde->eq << endl;
                  cout << "   PDE unknown field(s): ";
                  for (int i=0; i<_pde->nb_fields-1; ++i)
                     cout << _pde->fn[i] << ", ";
                  cout << _pde->fn[nb_fields-1] << endl;
                  cout << "   PDE space discretization: " << _pde->spD << endl;
                  cout << "   PDE linear solver: " << rLs[_pde->ls] << endl;
	       }
               _ret = 0;
	    }
	    else {
               cout << "Error: No PDE created." << endl;
               *_ofl << "In rita>pde>end>: No PDE data created." << endl;
	    }
            *_ofh << "  end" << endl;
            _nb_pde++;
            _ret = 0;
            return;

         case 17:
         case 18:
            if (_ret && _pde!=nullptr) {
               delete _pde; 
               _pde = nullptr;
            }
            _ret = 100;
            return;

         case 19:
            _ret = 200;
            return;

         case -2:
         case -3:
         case -4:
            break;

         default:
            cout << "Unknown Command: " << _cmd->token() << endl;
            cout << "Available commands: field, coef, in, bc, bf, sf, space, ls, nls, clear, end, <" << endl;
            cout << "Global commands:    help, ?, set, quit, exit" << endl;
            *_ofl << "In rita>pde>: Unknown Command " << _cmd->token() << endl;
            break;
      }
   }
   _ret = 0;
}

} /* namespace RITA */
