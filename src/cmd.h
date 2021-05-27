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

                        Definition of class Cmd

  ==============================================================================*/

#pragma once

#define USE_CTRL_D

#include <stdlib.h>
#include <ctype.h>
#include <vector>
using std::vector;

#include <iostream>
using std::cout;
using std::endl;
using std::ostream;

#include <sstream>
#include <fstream>
using std::istringstream;
using std::ifstream;

#define SHORT_HAND "$"

#ifdef USE_CTRL_D
#include <termios.h> 
#include <csignal> 
#include <cstdlib>
#endif

#include <string>
using std::string;

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace RITA {

/*! \file cmd.h
 *  \brief Definition file for class cmd.
 */

/*! \class cmd
 *  \brief 
 *
 * \author Rachid Touzani
 * \copyright GNU Public License
 */

class cmd {

 private:
   ifstream *_is;
   string _prompt, _buffer, _arg, _tok;
   vector<string> _toks;
   bool _command, _comment;
   vector<string> _word;
   size_t _ind;
   int _nb_args, _script_line_nb;
   bool isValidNumber(const string& s);
   bool isNumeric(const string& s);
   int split();
#ifdef USE_CTRL_D
   struct termios _old_termios, _new_termios;
   static void handler(int sig);
#endif
   
 public:

   cmd();
   cmd(ifstream& is);
   ~cmd();
   void setIFStream(ifstream *is) { _is = is; _script_line_nb = 0; }
   ifstream *getIFStream() const { return _is; }
   int readline(string p="");
   void setPrompt(string p) { _prompt = p; }
   void setErrorMsg(string e);
   int setNbArg(int n, const string& s="", int opt=0);
   int setNbArg(int n, int m, const string& s="", int opt=0);
   void setNbArg() { _nb_args = -1; }
   int getNbArgs() const { return _nb_args; }
   int get(string& s);
   int get(int &i);
   int get(double& d);
   int get(const vector<string>& kw, string& s);
   int getKW(const vector<string> &kw);
   int getScriptLineNb() const { return _script_line_nb; }
   string token() const { return _word[_ind-1]; }
   string buffer() const { return _buffer; }
   string string_token() const { return _tok; }
   double double_token() const { return stod(_tok); }
   int int_token() const { return stoi(_tok); }
   string string_token(int i) const { return _toks[i]; }
   double double_token(int i) const { return stod(_toks[i]); }
   int int_token(int i) const { return stoi(_toks[i]); }
   int getArg(string delimiter="=");
   int getArgs(int &nb, string del1="=", string del2=",");
   string Arg() const { return _arg; }
   void set(const vector<string>& arg);
   const vector<string> *_kw;
   string& trim(string& s, const string& chars = "\t\n\v\f\r ");
   int find_kw(const string &arg);

   friend class rita;
};

} /* namespace RITA */

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
