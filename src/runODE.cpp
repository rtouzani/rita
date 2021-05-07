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

                         Implementation of function runODE

  ==============================================================================*/

#include "rita.h"
#include "data.h"
#include "cmd.h"
#include "configure.h"

using std::cout;
using std::exception;

namespace RITA {

int rita::runODE()
{
   bool field_ok=false;
   int size=1, ret=0, count_fct=0, count_field=0, count_def=0, count_init=0, ind=-1;
   _init_time = 0., _time_step=0.1, _final_time=1.;

   vector<string> def, name, var;
   vector<double> init;
   _ode->isSet = false;
   _analysis_type = TRANSIENT;
   string str="", var_name="y", scheme="forward-euler";
   const vector<string> kw {"help","?","set","size","func$tion","def$inition","field","var$iable","init$ial",
                            "final$-time","time-step","scheme","summary","clear","end","<","quit","exit","EXIT"};
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

         case 6:
         case 7:
            var_name = _cmd->string_token();
            field_ok = true;
            count_field++;
            break;

         case 8:
            init.push_back(_cmd->double_token());
            count_init++;
            break;

         case 9:
            _final_time = _cmd->double_token();
            break;

         case 10:
            _time_step = _cmd->double_token();
            break;

         case 11:
            scheme = _cmd->string_token();
            break;

         default:
            msg("ode>","Unknown argument: "+_cmd->Arg());
	    return 1;
      }
   }

   if (_nb_args>0) {
      if (size<=0) {
         msg("ode>","Illegal size value.");
         return 1;
      }
      if (count_fct && count_field) {
         msg("ode>","Function already defined in data module.");
         return 1;
      }
      if (count_field>1) {
         msg("ode>","Only one variable must be defined for an ode system.");
         return 1;
      }
      if (count_fct && count_def) {
         msg("ode>","Function already defined.");
         return 1;
      }
      if (count_fct>size || count_def>size) {
         msg("ode>","Number of function names is larger than system size.");
         return 1;
      }
      if (_nb_ode>0 && count_field==0) {
         msg("ode>","No variable defined as unknown for new ode system.");
         return 1;
      }
      if (!field_ok) {
         msg("ode>","Missing a variable name.");
         return 1;
      }
      if (count_init>size) {
         msg("ode>","Number of initial conditions is larger than system size.");
         return 1;
      }
      if (count_fct>0 && count_def<size-1) {
         msg("ode>","Number of function definitions is larger than system size.");
         return 1;
      }
      if (count_fct && size>1) {
         msg("ode>","The option 'function' is not available with an ODE system.");
         return 1;
      }
      if (count_init<size) {
         for (int i=count_init; i<size; ++i)
            init.push_back(0.);
      }
      *ofh << "ode";
      _ode->theFct.resize(size);
      if (count_fct>0) {
         _ode->isFct = true;
         if (count_fct<=size-1) {
            msg("ode>","Number of function names is lower than system size.");
            return 1;
         }
         for (int k=0; k<size; ++k) {
            ind = _data->checkFct(name[k]);
            if (ind==-1) {
               msg("ode>","Non defined function "+name[k]);
               return 1;
            }
            if (_ode->theFct[k].set(name[k],_data->theFct[ind]->expr,_data->theFct[ind]->var,1)) {
               msg("ode>","Error in function evaluation: "+_ode->theFct[k].getErrorMessage());
               return 1;
            }
            *ofh << " function=" << name[k];
         }
      }
      else {
         _ode->isFct = false;
         *ofh << " var=" << var_name;
         var.resize(size+1);
         var[0] = "t";
         var[1] = var_name;
         if (size>1) {
            for (int i=1; i<=size; ++i)
               var[i] = var_name + to_string(i+1);
         }
         for (int i=0; i<size; ++i) {
            if (_ode->theFct[i].set(name[i],def[i],var,1)) {
               msg("ode>","Error in function evaluation: "+_ode->theFct[i].getErrorMessage());
               return 1;
            }
         }
         _ode->ind_fct = ind;
         for (int j=0; j<size; ++j)
            *ofh << " definition=" << def[j];
      }
      _ode->size = size;
      _ode->isSet = true;
      _ode->log = false;
      _ifield = _data->addField(var_name,GIVEN_SIZE,size);
      _data->FieldEquation[_ifield] = _ieq;
      _ode->field = _ifield;
      _ode->fn = var_name;
      _ode->log = false;
      _data->FieldType[_ifield] = ODE_EQ;
      _nb_fields = _data->getNbFields();
      for (const auto& v: init) {
         *ofh << " init=" << v;
         _ode->y.push_back(v);
      }
      _ode->scheme = _sch[scheme];
      *ofh << " scheme=" << scheme;
      *ofh << " time-step=" << _time_step << " final-time=" << _final_time << endl;
      _nb_ode++;
   }
   
   else {
      *ofh << "ode" << endl;
      while (1) {
         if (_cmd->readline("rita>ode> ")<0)
            continue;
         switch (_key=_cmd->getKW(kw)) {

            case 0:
            case 1:
               _cmd->setNbArg(0);
               cout << "\nAvailable Commands\n";
               cout << "size:       Size of differential system: (Number of equations)\n";
               cout << "function:   Give already defined function defining ode\n";
               cout << "definition: Give expression defining ode\n";
               cout << "variable:   Variable (field) name as unknown of the equation\n";
               cout << "initial:    Give an initial condition\n";
               cout << "final-time: Give final time\n";
               cout << "time-step:  Give time step\n";
               cout << "scheme:     Time integration scheme\n";
               cout << "summary:    Summary of ODE attributes\n";
               cout << "clear:      Remove ODE from model\n";
               cout << "end or <:   go back to higher level" << endl;
               break;

            case 2:
               setConfigure();
               break;

            case 3:
               if (_cmd->setNbArg(1,"Size of differential system to be given.")) {
                  msg("ode>size>","Missing system size.","",1);
                  break;
               }
               if (!_cmd->get(size))
                  *ofh << "  size " << size << endl;
               _ret = 0;
               break;

            case 4:
               if (count_def>0) {
                  msg("ode>function>","Function already defined by an expression.");
                  break;
               }
               if (count_fct==size)
                  msg("ode>function>","Too many functions defining ODE.");
               if (_cmd->setNbArg(1,"Function F to define equation y'(t) = F(t,y(t)) to be given.")) {
                  msg("ode>function>","Missing function expression.","",1);
                  break;
               }
               for (int i=0; i<_data->getNbFcts(); ++i) {
                  if (_data->theFct[i]->name==name[count_fct])
                     ind = i;
               }
               if (ind==-1) {
                  msg("ode>function>","Non defined function "+name[count_fct]);
                  ret = 1;
                  break;
               }
               ret = _cmd->get(name[count_fct]); 
               if (!ret) {
                  *ofh << "    function " << name[count_fct++] << endl;
                  field_ok = true;
                  _ode->theFct[count_fct].name = name[count_fct-1];
                  _ode->theFct[count_fct].set(_data->theFct[ind]->expr,_data->theFct[ind]->var);
                  count_fct++;
               }
               _ret = 0;
               break;

            case 5:
               if (count_fct>0) {
                  msg("ode>definition>","Function already defined by its name.");
                  break;
               }
               if (count_def==size)
                  msg("ode>definition>","Too many functions defining ODE.");
               if (_cmd->setNbArg(1,"Function F to define equation y'(t) = F(t,y(t)) to be given.")) {
                  msg("ode>definition>","Missing function expression.","",1);
                  break;
               }
               if (!_cmd->get(str)) {
                  def.push_back(str);
                  *ofh << "    definition " << str << endl;
                  count_def++;
               }
               _ret = 0;
               break;

            case 6:
            case 7:
               if (_ode->log) {
                  msg("ode>variable>","pde must be set first.");
                  break;
               }
               if (_cmd->setNbArg(1,"Give name of associated variable (field).")) {
                  msg("ode>variable>:","Missing name of associated variable (field).","",1);
                  break;
               }
               if (!_cmd->get(var_name)) {
                  field_ok = true, count_field++;
                  *ofh << "  variable " << var_name << endl;
               }
               else
                  msg("ode>variable>","Unknown variable (field) "+var_name);
               break;

            case 8:
               if (_cmd->setNbArg(size,"Initial conditions to be given.")) {
                  msg("ode>initial>","Missing initial conditions.","",1);
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
               else
                  msg("ode>initial>","Error in initial data.");
               _ret = 0;
               break;

            case 9:
               if (_cmd->setNbArg(1,"Final time to be given.")) {
                  msg("ode>final-time>","Missing final time value.","",1);
                  break;
               }
               if (!_cmd->get(_final_time))
                  *ofh << "    final-time " << _final_time << endl;
               _ret = 0;
               break;

            case 10:
               if (_cmd->setNbArg(1,"Missing time step.")) {
                  msg("ode>time-step>","Missing time step value.","",1);
                  break;
               }
               if (!_cmd->get(_time_step))
                  *ofh << "    time-step " << _time_step << endl;
               _ret = 0;
               break;

            case 11:
               if (_cmd->setNbArg(1,"Time integration scheme to be given.")) {
                  msg("ode>scheme>","Missing time integration scheme.","",1);
                  break;
               }
               if (!_cmd->get(scheme))
                  *ofh << "    scheme " << scheme << endl;
               _ret = 0;
               break;

            case 12:
               cout << "Summary of ODE attributes:\n";
               *ofh << "    summary" << endl;
               _ret = 0;
               break;

            case 13:
               _ode->log = false;
               cout << "ODE equation removed from model." << endl;
               *ofh << "    clear" << endl;
               _ret = 10;
               return _ret;

            case 14:
            case 15:
               if ((count_fct>0 && count_fct<size) || (count_fct==0 && count_def<size)) {
                  msg("ode>end>","Insufficient number of functions defining system.");
                  *ofh << "  end" << endl;
                  break;
               }
               if (!field_ok) {
                  msg("ode>end>","No variable defined for ode system.");
                  *ofh << "  end" << endl;
                  break;
               }
               if (count_fct && count_field) {
                  msg("ode>end>","Function already defined in data module.");
                  return 1;
               }
               if (count_field>1) {
                  msg("ode>end>","Only one variable must be defined for an ode system.");
                  return 1;
               }
               *ofh << "  end" << endl;
               var.clear();
               var.push_back("t");
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
               _ode->theFct.resize(size);
               _ode->size = size;
               _ode->ind_fct = ind;
               _ode->isFct = count_fct;
               _ode->y.resize(size);
               _ode->scheme = _sch[scheme];
               for (int j=0; j<size; ++j) {
                  _ode->y[j] = init[j];
                  if (!count_fct) {
                     if (_ode->theFct[j].set(def[j],var,1)) {
                        msg("ode>end>","Error in function evaluation: "+_ode->theFct[j].getErrorMessage());
                        break;
                     }
                  }
               }
               _data->addField(var_name,GIVEN_SIZE,size);
               _ode->field = _ifield;
               _data->FieldType[_ifield] = ODE_EQ;
               _nb_fields = _data->getNbFields();
               _ode->isSet = true;
               _ode->log = false;
               _data->FieldEquation[_ifield] = _ieq;
               _ode->isFct = false;
               if (count_fct)
                  _ode->isFct = true;
               _nb_ode++;
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
               msg("ode>","Unknown Command "+_cmd->token(),
                   "Available commands: size, function, definition, variable, initial, final-time, time-step,\n"
                   "                    scheme, summary, clear, end, <\n"
                   "Global commands:    help, ?, set, quit, exit");
               break;
         }
      }
   }
   _ret = 0;
   return _ret;
}

} /* namespace RITA */
