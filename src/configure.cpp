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

                     Implementation of class 'configure'

  ==============================================================================*/

#include "configure.h"

namespace RITA {

configure::configure(cmd *command)
          : _verb(1), _save_results(1), _his_file("rita.his"), _log_file("rita.log"),
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
   _ocf.open(".ritarc");
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
   _icf.open(".ritarc");
   if (_icf.fail())
      _icf.close();
   else {
      read();
      _icf.close();
      _ocf.open(".ritarc.backup");
      _ocf << "# rita configuration file" << endl;
      _ocf << "# " << currentDateTime() << "\n#\n";
      _ocf << "verbosity " << _verb << endl;
      _ocf << "save-results " << _save_results << endl;
      _ocf << "history-file " << _his_file << endl;
      _ocf << "log-file " << _log_file << endl;
      _ocf.close();
   }
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
   vector<string> kw;
   kw.push_back("verbosity");
   kw.push_back("save-results");
   kw.push_back("history-file");
   kw.push_back("log-file");
   int key=0;
   while (1) {
      int ret = com.readline();
      if (ret==-5)
         return 0;
      if (ret<0)
         continue;

      switch (key=com.getKW(kw)) {

         case 0:
            if (com.setNbArg(1,"Missing verbosity parameter."))
               return 0;
            com.get(_verb);
            break;

         case 1:
            if (com.setNbArg(1,"Missing result saving parameter."))
               return 0;
            com.get(_save_results);
            break;

         case 2:
            if (com.setNbArg(1,"Missing history file name."))
               return 0;
            com.get(_his_file);
            break;

         case 3:
            if (com.setNbArg(1,"Missing log file name."))
               return 0;
            com.get(_log_file);
            break;

         case -2:
            break;

         case -4:
            return 1;

         default:
            cout << "Unknown command in configuration file: " << com.token() << endl;
            return 0;
       }
   }
}


int configure::run()
{
   int ret = 0;
   string file, buffer;
   ifstream is;
   if (_cmd->setNbArg(2,"Give parameter to set.")) {
      cout << "Available Settings:\n";
      cout << "verbosity   : Give verbosity parameter\n";
      cout << "save-results: Give results saving parameter\n";
      cout << "history     : Give name of history file\n";
      cout << "log         : Give name of log file" << endl;
      return 0;
   }
   while (1) {
      switch (_key=_cmd->getKW(_kw)) {

         case 0:
         case 1:
            cout << "\nAvailable Settings:\n";
            cout << "verbosity   : Give verbosity parameter\n";
            cout << "save-results: Give results saving parameter\n";
            cout << "history     : Give name of history file\n";
            cout << "log         : Give name of log file" << endl;
            return 0;

         case 2:
            if (_cmd->setNbArg(2,"Give new verbosity parameter.")) {
               _ofl << "In rita>set>verbosity>: Missing verbosity value." << endl;
               return 0;
            }
            ret = _cmd->get(_verb);
            if (!ret)
               _ofh << "set verbosity " << _verb << endl;
            return 0;

         case 3:
            if (_cmd->setNbArg(2,"Give result saving parameter")) {
               _ofl << "In rita>set>save-results>: Missing result saving parameter." << endl;
               return 0;
            }
            ret = _cmd->get(_save_results);
            if (!ret)
               _ofh << "set save-results " << _save_results << endl;
            return 0;

         case 4:
            if (_cmd->setNbArg(2,"Give name of history file")) {
               _ofl << "In rita>set>history>: Missing name of history file." << endl;
               return 0;
            }
            file = _his_file;
            ret = _cmd->get(_his_file);
            if (!ret)
               _ofh << "set history " << _his_file << endl;
            _ofh.close();
            is.open(file);
            _ofh.open(_his_file.c_str());
            while (!is.eof()) {
               getline(is,buffer);
               _ofh << buffer << endl;
            }
            return 0;

         case 5:
            if (_cmd->setNbArg(2,"Give name of log file")) {
               _ofl << "In rita>set>log>: Missing name of log file." << endl;
               return 0;
            }
            file = _log_file;
            ret = _cmd->get(_log_file);
            if (!ret)
               _ofh << "set log " << _log_file << endl;
            _ofl.close();
            is.open(file);
            _ofl.open(_log_file.c_str());
            while (!is.eof()) {
               getline(is,buffer);
               _ofl << buffer << endl;
            }
            return 0;

         case 6:
         case 7:
            return 100;

         case 8:
            return 200;

         case -2:
            break;

         default:
            cout << "Unknown Setting: " << _cmd->token() << endl;
            cout << "Available settings: verbosity, save-results, history, log" << endl;
            _ofl << "In rita>set>: Unknown setting: " << _cmd->token() << endl;
            return 0;
       }
   }
}

} /* namespace RITA */
