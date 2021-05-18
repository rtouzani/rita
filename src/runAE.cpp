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

                         Implementation of function runAE

  ==============================================================================*/

#include "rita.h"
#include "data.h"
#include "cmd.h"
#include "configure.h"

using std::cout;
using std::exception;

namespace RITA {

int rita::runAE()
{
   string str = "", var_name = "x", nls = "newton";
   bool field_ok = false;
   int ret=0, size=1, ind=-1, n=0;
   int count_field=0, count_fct=0, count_def=0, count_J=0, count_init=0;
   _ret = 0;
   vector<string> def, var, name;
   vector<double> init;
   _ae->isSet = false;
   OFELI::Vect<string> J(1,1);
   if (_analysis_type==NONE)
      _analysis_type = STEADY_STATE;

   const static vector<string> kw {"help","?","set","size","func$tion","def$inition","jacob$ian","init",
                                   "field","var$iable","nls","summary","clear","remove","end","<","quit",
                                   "exit","EXIT"};
   _cmd->set(kw);
   for (int k=0; k<_nb_args; ++k) {

      switch (_cmd->getArg("=")) {

         case 3:
            size = _cmd->int_token();
            break;

         case 4:
            name.push_back(_cmd->string_token());
            count_fct++, field_ok = true;
            break;

         case 5:
            def.push_back(_cmd->string_token());
            name.push_back(" ");
            count_def++;
            break;

         case 7:
            init.push_back(_cmd->double_token());
            count_init++;
            break;

         case 8:
         case 9:
            var_name = _cmd->string_token();
            field_ok = true;
            break;

         case 10:
            nls = _cmd->string_token();
            break;

         default:
            msg("algebraic>:","Unknown argument: "+_cmd->Arg());
            return 1;
      }
   }

   if (_nb_args>0) {
      if (size<=0) {
         msg("algebraic>","Illegal size value.");
         return 1;
      }
      if (count_fct && count_field) {
         msg("algebraic>","Function already defined in data module.");
         return 1;
      }
      if (count_field>1) {
         msg("algebraic>","Only one variable must be defined for an algebraic system.");
         return 1;
      }
      if (count_fct && count_def) {
         msg("algebraic>","Function already defined.");
         return 1;
      }
      if (count_fct>size || count_def>size) {
         msg("algebraic>","Number of function names is larger than system size.");
         return 1;
      }
      if (!field_ok) {
         msg("algebraic>","Missing a variable name.");
         return 1;
      }
      if (count_init>size) {
         msg("algebraic>","Number of initial guesses is larger than system size.");
         return 1;
      }
      if (count_fct>0 && count_def<size-1) {
         msg("algebraic>","Number of function definitions is larger than system size.");
         return 1;
      }
      if (count_init<size) {
         for (int i=count_init; i<size; ++i)
            init.push_back(0.);
      }
      *ofh << "algebraic";
      _ae->theFct.resize(size);
      if (count_fct) {
         _ae->isFct = true;
         if (count_fct<=size-1) {
            msg("algebraic>","Number of function names is lower than system size.");
            return 1;
         }
         for (int k=0; k<size; ++k) {
            ind = _data->checkFct(name[k]);
            if (ind==-1) {
               msg("algebraic>","Non defined function "+name[k]);
               return 1;
            }
            if (_ae->theFct[k].set(name[k],_data->theFct[ind]->expr,_data->theFct[ind]->var,1)) {
               msg("algebraic>","Error in function evaluation: "+_ae->theFct[k].getErrorMessage());
               return 1;
            }
            *ofh << " function=" << name[k];
         }
      }
      else {
         _ae->isFct = false;
         *ofh << " var=" << var_name;
         var.resize(size);
         var[0] = var_name;
         if (size>1) {
            for (int i=0; i<size; ++i)
               var[i] = var_name + to_string(i+1);
         }
         for (int i=0; i<size; ++i) {
            name[i] = "ff-" + to_string(++_data->iifct);
            _data->addFunction(name[i],def[i],var);
            _ae->theFct[i].set(name[i],def[i],var,1);
         }
         _ae->ind_fct = ind;
         for (int j=0; j<size; ++j)
            *ofh << " definition=" << def[j];
      }
      _ae->size = size;
      _ae->J.setSize(size,size);
      _ae->isSet = true;
      _ae->log = false;
      _data->addField(var_name,GIVEN_SIZE,size);
      _data->FieldEquation[_ifield] = _ieq;
      _nb_fields = _data->getNbFields();
      _ae->field = _ifield;
      _data->FieldType[_ifield] = ALGEBRAIC_EQ;
      for (const auto& v: init) {
         *ofh << " init=" << v;
         _ae->y.push_back(v);
      }
      _ae->nls = NLs[nls];
      _ae->isFct = false;
      *ofh << " nls=" << nls << endl;
      _nb_ae++;
   }

   else {
      *ofh << "algebraic " << endl;
      count_fct = count_init = count_field = count_def = count_J = 0;
      while (1) {
         if (_cmd->readline("rita>algebraic> ")<0)
            continue;
         switch (_key=_cmd->getKW(kw)) {

            case 0:
            case 1:
               _cmd->setNbArg(0);
               cout << "\nAvailable Commands:\n";
               cout << "size:       Size of algebraic system: (Number of equations)\n";
               cout << "function:   Name of already defined function\n";
               cout << "definition: Give function expression F to define the algebraic equation F(x)=0\n";
               cout << "jacobian:   Define jacobian of mapping (for the Newton's algorithm)\n";
               cout << "variable:   Variable (field) name as unknown of the equation\n";
               cout << "init:       Initial guess for iterations\n";
               cout << "nls:        Nonlinear equation iteration solver\n";
               cout << "summary:    Summary of Algebraic equation attributes\n";
               cout << "clear:      Remove equation from model\n";
               cout << "end or <:   go back to higher level" << endl;
               break;

            case 2:
               setConfigure();
               break;

            case 3:
               if (_cmd->setNbArg(1,"Size of algebraic system to be given.")) {
                  msg("algebraic>size>","Missing system size.","",1);
                  break;
               }
               if (!_cmd->get(size))
                  *ofh << "  size " << size << endl;
               _ret = 0;
               break;

            case 4:
               if (count_def) {
                  msg("algebraic>function>","Function already defined.");
                  break;
               }
               if (count_fct==size) {
                  msg("algebraic>function>","Number of functions is larger than system size.");
                  return 1;
               }
               if (_cmd->setNbArg(1,"Name of function to be given.")) {
                  msg("algebraic>function>","Missing function name.");
                  break;
               }
               str = _cmd->string_token();
               for (int i=0; i<_data->getNbFcts(); ++i) {
                  if (_data->theFct[i]->name==str)
                     ind = i;
	       }
               if (ind==-1) {
                  msg("algebraic>function>","Non defined function "+str);
                  ret = 1;
                  break;
               }
               if (!ret) {
                  *ofh << "  function " << str << endl;
                  name[count_fct++] = str;
                  field_ok = true;
                  if (_ae->theFct[count_fct].set(str,_data->theFct[ind]->expr,_data->theFct[ind]->var,1)) {
                     msg("algebraic>function>","Error in function evaluation: "+_ae->theFct[count_fct].getErrorMessage());
                     return 1;
                  }
               }
               _ret = 0;
               break;

            case 5:
               if (count_fct>0) {
                  msg("algebraic>definition>","Function already defined by its name.");
                  break;
               }
               if (count_def==size) {
                  msg("algebraic>definition>","Too many functions defining system.");
                  break;
               }
               if (_cmd->setNbArg(1,"Function F to define equation F(x)=0 to be given.")) {
                  msg("algebraic>definition>","Missing function expression.","",1);
                  break;
               }
               if (_cmd->get(str)) {
                  msg("algebraic>definition>","Error in equation definition.");
                  break;
               }
               else {
                  def.push_back(str);
                  count_def++;
                  *ofh << "  definition " << str << endl;
               }
               _ret = 0;
               break;

            case 6:
               if (count_J==size) {
                  msg("algebraic>jacobian>","Too many functions defining jacobian.");
                  break;
               }
               if (_cmd->setNbArg(size,"Partial derivatives of function defining equation to be given.")) {
                  msg("algebraic>jacobian>","Missing partial derivatives of function defining equation.","",1);
                  break;
               }
               for (int i=1; i<=size; ++i)
                  ret = _cmd->get(J(count_J+1,i));
               if (ret==0) {
                  *ofh << "  jacobian " << endl;
                  for (int i=1; i<=size; ++i)
                     *ofh << J(count_J+1,i) << "  ";
                  count_J++;
               }
               _ret = 0;
               break;

            case 7:
               if (_cmd->setNbArg(size,"Initial guesses to be given.")) {
                  msg("algebraic>initial>","Missing initial guesses.","",1);
                  break;
               }
               ret = 0;
               init.resize(size);
               count_init++;
               for (int i=0; i<size; ++i)
                  ret += _cmd->get(init[i]);
               if (!ret) {
                  *ofh << "  initial  ";
                  for (const auto& v: init)
                     *ofh << v << "  ";
                  *ofh << endl;
               }
               else {
                  msg("algebraic>initial>","Error in initial guess data.");
                  break;
               }
               _ret = 0;
               break;

            case 8:
            case 9:
               if (_ae->log) {
                  msg("algebraic>variable>","Equation must be defined first.");
                  break;
               }
               if (_cmd->setNbArg(1,"Give name of associated field.")) {
                  msg("algebraic>variable>","Missing name of associated variable.","",1);
                  break;
               }
               if (_cmd->get(var_name)) {
                  msg("algebraic>variable>","Unknown variable or field "+var_name);
                  break;
               }
               else {
                  if (size>1 && _verb) {
                     cout << "Unknowns are: ";
                     for (int i=1; i<size; ++i)
                        cout << var_name + to_string(i) << ", ";
                     cout << var_name + to_string(size) << endl;
                  }
                  count_field++, field_ok = true;
                  *ofh << "  variable " << var_name << endl;
               }
               break;

            case 10:
               if (_ae->log) {
                  msg("algebraic>nls>","equation must be set first.");
                  break;
               }
               if (_cmd->setNbArg(1,"Nonlinear solver to be supplied.",1)) {
                  msg("algebraic>nls>","Missing nonlinear solver data.","",1);
                  break;
               }
               ret = _cmd->get(nls);
               if (!ret)
                  *ofh << "  nls " << nls << endl;
               _ret = 0;
               break;

            case 11:
               cout << "Summary of algebraic system:\n";
               *ofh << "      summary" << endl;
               for (int e=0; e<_nb_eq; ++e) {
                  if (_eq_type[e]==ALGEBRAIC_EQ) {
                     odae *ae = ALGEBRAIC[e];
                     cout << "Equation Label: " << e << endl;
                     if (size>1)
                        cout << "Size:           " << size << endl;
                     if (_ae->isFct) {
                        if (size==1)
                           cout << "Equation defined by function: " << ae->theFct[0].name << endl;
                        else {
                           for (int i=0; i<size; ++i)
                              cout << "Equation: " << i+1 << ", defined by function: " << _ae->theFct[i].name << endl;
                        }
                     }
                     else {
                        if (size==1) {
                           cout << "Equation defined by: " << ae->theFct[0].expr << endl;
                           cout << "Variable is          " << ae->theFct[0].var[0] << endl;
                        }
                        else {
                           for (int i=0; i<size; ++i)
                              cout << "Equation: " << i+1 << ", defined by: " << ae->theFct[i].expr << endl;
                           for (int i=0; i<size; ++i)
                              cout << "Variable " << i+1 << " is " << ae->theFct[0].var[i] << endl;
                        }
                     }
                  }
               }
               _ret = 0;
               break;

            case 12:
               if (_ae->log) {
                  cout << "Algebraic equation has been removed from model." << endl;
                  *ofh << "  clear" << endl;
               }
               else {
                  msg("algebraic>clear>","No equation to remove.");
                  break;
               }
               _ret = 10;
               return _ret;

            case 13:
               if (_cmd->setNbArg(1,"Size of algebraic system to be given.")) {
                  msg("algebraic>remove>","Missing system size.","",1);
                  break;
               }
               if (!_cmd->get(n))
                  *ofh << "  remove " << n << endl;
               if (n<=0 || n>_nb_eq) {
                  msg("algebraic>remove>","Illegal equation label: "+to_string(n));
                  break;
               }
               for (int e=0; e<_nb_eq; ++e) {
                  if (e==n-1) {
                     if (_eq_type[e]==ALGEBRAIC_EQ) {
                        ALGEBRAIC[e] = nullptr;
                        _ifield--;
                     }
                  }
                  else {
                     msg("algebraic>remove>","Equation is not algebraic one.");
                     break;
                  }
               }
               _ret = 10;
               return _ret;

            case 14:
            case 15:
               _cmd->setNbArg(0);
               if (!field_ok) {
                  msg("algebraic>end>","Missing a variable name.");
                  return 1;
               }
               if ((count_fct>0 && count_fct<size) || (count_fct==0 && count_def<size)) {
                  msg("algebraic>end>","Insufficient number of functions defining system.");
                  *ofh << "  end" << endl;
                  break;
               }
               if (!field_ok) {
                  msg("algebraic>end>","No variable defined for algebraic system.");
                  *ofh << "  end" << endl;
                  break;
               }
               if (count_fct && count_field) {
                  msg("algebraic>end>","Function already defined in data module.");
                  return 1;
               }
               if (count_field>1) {
                  msg("algebraic>end>","Only one variable must be defined for an algebraic system.");
                  return 1;
               }
               *ofh << "  end" << endl;
               var.clear();
               if (size==1)
                  var.push_back(var_name);
               else {
                  for (int i=1; i<=size; ++i)
                     var.push_back(var_name+to_string(i));
               }
               if (!count_init) {
                  init.resize(size);
                  for (int i=0; i<size; ++i)
                     init[i] = 0.;
               }
               _ae->theFct.resize(size);
               if (!count_J) {
                  J.setSize(size,size);
                  for (int i=1; i<=size; ++i)
                     J(i,i) = 1.;
               }
               _ae->size = size;
               _ae->ind_fct = ind;
               _ae->isFct = count_fct;
               _ae->y.resize(size);
               for (int j=0; j<size; ++j) {
                  _ae->y[j] = init[j];
                  if (!count_fct) {
                     if (_ae->theFct[j].set(def[j],var,1)) {
                        msg("algebraic>end>","Error in function evaluation: "+_ae->theFct[j].getErrorMessage());
                        break;
                     }
                  }
               }
               if (count_J==size) {
                  for (int i=1; i<=size; ++i)
                     for (int j=1; j<=size; ++j)
                        _ae->J(i,j) = J(i,j);
               }
               _ifield = _data->addField(var_name,GIVEN_SIZE,size);
               _nb_fields = _data->getNbFields();
               _ae->field = _ifield;
               _data->FieldType[_ifield] = ALGEBRAIC_EQ;
               _ae->isSet = true;
               _ae->log = false;
               _ae->nls = NLs[nls];
               _data->FieldEquation[_ifield] = _ieq;
               _ae->isFct = false;
               if (count_fct)
                  _ae->isFct = true;
               _nb_ae++;
               _ret = 0;
               return _ret;
 
            case 16:
            case 17:
               return 100;

            case 18:
               return 200;

            case -2:
            case -3:
            case -4:
               break;

            default:
               msg("algebraic>","Unknown Command "+_cmd->token(),
                   "Available commands: size, function, definition, jacobian, init,\n"
                   "                    variable, nls, summary, clear, remove, end, <\n"
                   "Global commands:    help, ?, set, quit, exit");
               break;
         }
      }
   }
   _ret = 0;
   return _ret;
}

} /* namespace RITA */
