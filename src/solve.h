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

                            Definition of class 'solve'

  ==============================================================================*/

#ifndef __SOLVE_H
#define __SOLVE_H

#include <string>
#include <iostream>
using std::string;

#include "rita.h"
#include "linear_algebra/Matrix.h"
#include "solvers/LinearSolver.h"
#include "solvers/OptSolver.h"
#include "solvers/EigenProblemSolver.h"
#include "io/IOField.h"
#include "io/saveField.h"
#include "io/Fct.h"


namespace RITA {

class configure;
class data;
class optim;

class solve
{

 public:

    solve();
    solve(rita *r, cmd *command, configure *config);
    ~solve() { }
    void setVerbose(int verb) { _verb = verb; }
    int getVerbose() const { return _verb; }
    int getNbEq() const { return _nb_eq; }
    void set(cmd* com) { _cmd = com; }
    void set(OFELI::Mesh *ms);
    void setSave(int s) { _save_results = s; }
    int run();
    int ret() const { return _ret; }

 private:

    rita *_rita;
    bool _set_analytic, _solved, _phase;
    int _verb, _key, _ret, _save_results;
    int _nb_fields, _nb_eq;
    vector<string> _analytic_exp, _var;
    vector<int> _fformat, _isave;
    vector<string> _save_file, _phase_file;
    OFELI::Mesh *_theMesh;
    OFELI::Fct _theFct;
    configure *_configure;
    data *_data;
    optim *_optim;
    cmd *_cmd;

    void save();
    void display();
    int plot();
    int run_steady();
    int run_transient();
    int run_optim();
    int run_eigen();
    void get_error(int eq, int i);
    void setAnalytic();
    vector<string> _kw_solve = {"help","?","set","run","save","display","plot","analytic","error",
                                "post","end","<","quit","exit","EXIT"};
    vector<string> _kw_save = {"help","?","set","field","format","freq$uency","phase",
                               "file","end","<","quit","exit","EXIT"};
    vector<string> _kw_format = {"ofeli","gmsh","gnuplot","vtk","tecplot","matlab"};
    vector<string> _kw_display = {"help","?","field","sol$ution","phase","iter$ation","conv$ergence",
                                  "end","<","quit","exit","EXIT"};
    map<string,int> _ff = {{"gmsh",GMSH},{"gnuplot",GNUPLOT},{"vtk",VTK},{"tecplot",TECPLOT},{"matlab",MATLAB},
			   {"ofeli",OFELI_FF}};
    vector<string> _kw_analytic = {"help","?","set","eq$uation","comp$onent","exp$ression","end","<","quit","exit"};
};

} /* namespace RITA */

#endif
