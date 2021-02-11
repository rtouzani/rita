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

                        Definition of class 'transient'

  ==============================================================================*/

#ifndef __TRANSIENT_H
#define __TRANSIENT_H

#include "mesh/Mesh.h"
#include "solvers/LinearSolver.h"
#include "solvers/ODESolver.h"
#include "solvers/TimeStepping.h"
#include "solvers/NLASSolver.h"
#include "OFELI.h"
#include "rita.h"
#include "solve.h"
#include <map>

namespace RITA {

class equa;
struct odae;
 
class transient
{

 public:

    transient(rita *r);
    ~transient();
    void setLinearSolver(OFELI::Iteration ls, OFELI::Preconditioner prec);
    void setSave(vector<int>& isave, vector<int>& fformat, vector<string>& save_file,
                 bool phase, vector<string>& phase_file);
    int run();

 private:

    bool _ts_allocated, _ode_allocated, _nlas_allocated, _phase;
    rita *_rita;
    data *_data;
    double _init_time, _final_time, _time_step;
    int _nb_fields, _nb_eq, _rs;
    vector<int> *_fformat, *_isave;
    vector<string> *_save_file, *_phase_file;
    OFELI::NLASSolver *_nlas;
    OFELI::ODESolver *_ode;
    OFELI::TimeStepping *_ts;
    vector<int> *_eq_type;
    std::vector<equa *> _pde_eq;
    std::vector<odae *> _algebraic_eq, _ode_eq;
    int setPDE(int e);
};

} /* namespace RITA */

#endif
