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

                      Implementation of class 'approx'

  ==============================================================================*/


#include "approximation.h"

namespace RITA {

approximation::approximation(rita* r, cmd* command, configure* config)
              : _rita(r), _configure(config), _cmd(command)
{
}


approximation::~approximation()
{
}


int approximation::run()
{
   _rita->_analysis_type = APPROXIMATION;
   string file="";
   int nb=0, nl=0, nb_tab=0;
   int file_count=0, lagrange_count=0, bspline_count=0, fitting_count=0, bezier_count=0;
   int nurbs_count=0, approx_count=0;
   const vector<string> _kw = {"help","?","set","file","lagrange","fitting","bspline","bezier","nurbs"
                               "end","<","quit","exit","EXIT"};
   _cmd->set(_kw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args<1) {
      cout << "Error: No argument for command." << endl;
      *_ofl << "In rita>approximation>: No argument for command." << endl;
      return 1;
   }
   for (int k=0; k<nb_args; ++k) {

      int n = _cmd->getArgs(nb);
      switch (n) {

         case 3:
            file = _cmd->string_token(0);
            nb_tab = 1;
            file_count++;
            break;

         case 4:
            nl = _cmd->int_token(0);
            lagrange_count++, approx_count++;
            break;

         case 5:
            approx_count++;
            break;

         case 6:
            approx_count++;
            break;

         case 7:
            fitting_count++, approx_count++;
            break;

         case 8:
            bspline_count++, approx_count++;
            break;

         case 9:
            bezier_count++, approx_count++;
            break;

         case 10:
            nurbs_count++, approx_count++;
            break;

         default:
            cout << "Error: Unknown argument: " << _cmd->Arg() << endl;
            *_ofl << "In rita>approximation>: Unknown argument: " << _cmd->Arg() << endl;
            return 1;
      }
   }

   if (nb_args>0) {
      if (file_count==0) {
         cout << "Error: No data file given." << endl;
         *_ofl << "In rita>approximation>: No data file given." << endl;
         return 1;
      }
      *_ofh << "approximation";
      tab.setFile(file);
      *_ofh << " file=" << file;
      if (lagrange_count++) {
         *_ofh << "lagrange=" << nl;
      }
      if (approx_count>1) {
         cout << "Error: More than one approximation method given." << endl;
         *_ofl << "In rita>approximation>: More than one approximation method given." << endl;
         return 1;
      }
   }
   /*   else {
      *_ofh << "approximation " << endl;
      while (1) {
      if (_cmd->readline("rita>approximation> ")<0)
         continue;
      switch (_key=_cmd->getKW(_kw_approx)) {

         case 0:
         case 1:
            _cmd->setNbArg(0);
            cout << "\nAvailable Commands:\n";
            cout << "file:         File containing tabulated data to approximate\n";
            cout << "lagrange:     \n";
            cout << "least-square: Least Square fitting\n";
            cout << "summary:      Summary of approximation problem attributes\n";
            cout << "clear:        Remove problem\n";
            cout << "end or <:     go back to higher level" << endl;
            break;

         case 2:
            setConfigure();
            break;

         case 3:
            if (!_cmd->get(size))
               *_ofh << "    size " << size << endl;
            _ret = 0;
            break;

         case 4:
            _ret = 0;
            break;

         case 5:
            _ret = 0;
            break;

         case 7:
            break;

         case 8:
            _ret = 0;
            break;

         case 9:
            break;

         case 11:
            cout << "Summary of optimization problem attributes:\n";
            *_ofh << "      summary" << endl;
            _ret = 0;
            break;
         case 13:
         case 14:
            _cmd->setNbArg(0);
            if (fct_ok<size) {
               cout << "Error: Insufficient number of functions defining system." << endl;
               *_ofl << "In rita>optimization>end>: Insufficient number of functions defining system." << endl;
               break;
            }
            if (field_ok==0) {
               cout << "Error: No field defined for optimization problem." << endl;
               *_ofl << "In rita>optimization>end>: No field defined for optimization problem." << endl;
               break;
            }
            *_ofh << "      end" << endl;
            _ret = 0;
            return;

         case 15:
         case 16:
            _ret = 100;
            return;

         case -2:
         case -3:
         case -4:
            break;

         default:
            cout << "Unknown Command: " << _cmd->token() << endl;
            cout << "Available commands: size, objective, gradient, hessian, constraint, init, field, algorithm" << endl;
	    cout << "                    summary, clear, end, <" << endl;
            cout << "Global commands:    help, ?, set, quit, exit" << endl;
            *_ofl << "In rita>optimization>: Unknown Command " << _cmd->token() << endl;
         break;
      }
   }
   _ret = 0;*/
   return 0;
}


int approximation::go()
{
   switch (method) {

      case LAGRANGE:
         break;

      case PIECEWISE_LAGRANGE:
         break;

      case HERMITE:
         break;

      case FITTING:
         break;

      case BSPLINE:
         break;

      case BEZIER:
         break;

      case NURBS:
         break;
   }
   return 0;
}

} /* namespace RITA */
