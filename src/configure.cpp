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
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

  ==============================================================================

                     Implementation of class 'configure'

  ==============================================================================*/

#include "configure.h"

namespace RITA {

configure::configure(cmd *command)
          : _verb(1), _save_results(1), _his_file(".rita.his"), _log_file(".rita.log"),
            _cmd(command)
{
   init();
}


configure::~configure()
{
   save();
}


void configure::save()
{
   _ocf.open((_HOME+"/.rita").c_str());
   _ocf << "# rita configuration file" << endl;
   _ocf << "# " << currentDateTime() << "\n#\n";
   _ocf << "verbosity " << _verb << endl;
   _ocf << "save-results " << _save_results << endl;
   _ocf << "history-file " << _his_file << endl;
   _ocf << "log-file " << _log_file << endl;
   _ocf.close();
}


void configure::init()
{
   _HOME = getenv("HOME");
   _icf.open((_HOME+"/.rita").c_str());
   if (_icf.fail()) {
      _icf.close();
      _ocf.open((_HOME+"/.rita").c_str());
   }
   else {
      read();
      _icf.close();
      _ocf.open((_HOME+"/.rita.backup").c_str());
   }
   _ocf << "# rita configuration file" << endl;
   _ocf << "# " << currentDateTime() << "\n#\n";
   _ocf << "verbosity " << _verb << endl;
   _ocf << "save-results " << _save_results << endl;
   _ocf.close();
   _ofl.open(_log_file);
   _ofl << "# rita Log file" << endl;
   _ofl << "# " << currentDateTime() << "\n#\n";
   _ofh.open(_his_file);
   _ofh << "# rita History file" << endl;
   _ofh << "# " << currentDateTime() << "\n#\n";
}


int configure::read()
{
   cmd com(_icf);
   com.readline();
   com.set(_kw);
   int nb_args = com.getNbArgs();
   if (nb_args < 0)
      return 1;

   for (int i=0; i<nb_args; ++i) {
      int n = com.getArg();
      switch (n) {

         case 0:
            _verb = com.int_token();
            break;

         case 1:
            _save_results = com.int_token();
            break;

         case 2:
            _his_file = com.string_token();
            break;

         case 3:
            _log_file = com.string_token();
            break;

         default:
            cout << "Unknown Setting: " << com.token() << endl;
            cout << "Available settings: verbosity, save-results, history, log" << endl;
            _ofl << "In rita>set>: Unknown setting: " << com.token() << endl;
            return 1;
      }
   }
   return 0;
}


int configure::run()
{
   bool verb_ok=false, hist_ok=false, log_ok=false, save_ok=false;
   string hfile, lfile, buffer;
   ifstream is;
   _cmd->set(_kw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args < 0)
      return 1;
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg();
      switch (n) {

         case 0:
            _verb = _cmd->int_token();
            verb_ok = true;
            break;

         case 1:
            _save_results = _cmd->int_token();
            save_ok = true;
            break;

         case 2:
            hist_ok = true;
            hfile = _his_file;
            _his_file = _cmd->string_token();
            break;

         case 3:
            log_ok = true;
            lfile = _log_file;
            _log_file = _cmd->string_token();
            break;

         default:
            cout << "Unknown Setting: " << _cmd->token() << endl;
            cout << "Available settings: verbosity, save-results, history, log" << endl;
            _ofl << "In rita>set>: Unknown setting: " << _cmd->token() << endl;
            return 1;
       }
   }
   if (nb_args>0) {
      _ofh << "set";
      if (verb_ok) {
         if (_verb<0 || _verb>10) {
            cout << "Error: Illegal value of verbosity: " << _verb << endl;
            _ofl << "In rita>set>: Illegal value of verbosity: " << _verb << endl;
            return 1;
         }
         _ofh << " verbosity=" << _verb;
      }
      if (save_ok) {
         if (_save_results<0) {
            cout << "Error: Illegal value of save: " << _save_results << endl;
            _ofl << "In rita>set>: Illegal value of save: " << _save_results << endl;
            return 1;
         }
	 _ofh << " save-results=" << _save_results;
      }
      if (hist_ok) {
         _ofh.close();
         is.open(hfile);
         _ofh.open(_his_file.c_str());
         while (!is.eof()) {
            getline(is,buffer);
            _ofh << buffer << endl;
         }
      }
      if (log_ok) {
         _ofl.close();
         is.open(lfile);
         _ofl.open(_log_file.c_str());
         while (!is.eof()) {
            getline(is,buffer);
            _ofl << buffer << endl;
         }
      }
   }
   return 0;
}

} /* namespace RITA */
