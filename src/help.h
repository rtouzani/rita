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

                          Definition of class 'help'

  ==============================================================================*/

#pragma once

#include <string>
#include <iostream>

namespace RITA {

string H0 = "rita is an interactive application to solve modelling problems involving\n"
             "partial differential equations by the OFELI library.\n";
string H1 =  "rita is an interactive application to solve modelling problems involving\n"
             "partial differential equations by the OFELI library.\n"
             "rita is able to solve also basic numerical analysis problems.\n"
             "Solving a problem by rita consists in running various modes in order to achieve\n"
             "various steps of the solution procedure.\n\n"
             "rita can be either run interactively or using a script file. This latter option can be chosen\n"
             "simply either by executing 'rita script-file' or once in the rita execution by typing\n"
             "'load script-file'. The file 'script-file' must contain all commands to be executed,\n"
             "one on each line. It can contain comment lines (lines starting with a '#')\n\n"
             "Principal commands in rita are:\n"
             "data       Define data (fields, functions, ...).\n"
             "mesh       Define domain and generate its meshing\n"
             "stationary Set problem as stationary\n"
             "transient  Set problem as time dependent\n"
             "algebraic  Define an algebraic equation or system to solve\n"
             "ode        Define an ordinary differential equation (or system of equations) to solve\n"
             "pde        Define a partial differential equation (or system of equations) to solve\n"
             "optim      Set problem as an optimization one (Not yet implemented)\n"
             "eigen      Set problem as an eigenproblem (Not yet implemented)\n"
             "solve      Run the constructed model problem with all chosen options\n"
             "license    Print License Agreement of the software\n\n"
             "Global commands are commands that are available in all modes. These are:\n"
             "set        Modify configuration parameters for the present and future sessions\n"
             "help or ?  Display a help corresponding to current mode\n"
             "!          To execute a shell command, this one must be typed preceded by this mark\n"
             "end or <   Go back to higher level: Enables quiting a mode. Some commands in modes however\n"
             "           do this automatically.\n"
             "exit       exit the program.\n";

class help
{

 public:

    help() : _verb(1)
    {
       _H1 = "--------------------------------------------------------------------------------------------------\n"
             + H1 +
             "--------------------------------------------------------------------------------------------------\n";
    }

    ~help() { }

    int run()
    {
       if (_verb>0)
          std::cout << _H1 << std::endl;
       return 0;
    }

    void setVerbose(int verb) { _verb = verb; }

 private:

   int _verb;
   std::string _H1;

};

} /* namespace RITA */
