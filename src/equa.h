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

                           Definition of class 'equa'

  ==============================================================================*/

#ifndef __EQUAT_H
#define __EQUAT_H

#include <string>
#include <fstream>
#include <vector>
#include <map>
using std::map;

#include "data.h"

#include "equations/Equa_impl.h"
#include "equations/Equation_impl.h"
#include "Laplace.h"
#include "Therm.h"
#include "Solid.h"
#include "Fluid.h"
#include "Electromagnetics.h"
#include "io/Fct.h"

using namespace OFELI;

namespace RITA {

class cmd;
class rita;

class equa
{

 public:

    struct Log {
      bool pde, field, spd, ls, nl, mesh;
      Log() { pde = field = spd = ls = nl = mesh = false; }
      bool fail() { return (pde || field || spd || ls || nl || mesh); }
    };

    equa(rita *r);
    ~equa();

    enum pde_eq {
       LAPLACE,
       HEAT,
       WAVE,
       TRANSPORT,
       LINEAR_ELASTICITY,
       TRUSS,
       BEAM,
       INCOMPRESSIBLE_NAVIER_STOKES,
       COMPRESSIBLE_EULER,
       COMPRESSIBLE_NAVIER_STOKES,
       INCOMPRESSIBLE_POROUS_1PHASE,
       EDDY_CURRENTS,
       MAXWELL,
       HELMHOLTZ,
       EIKONAL
    };

    enum sdm {
       FD,
       FE_P1,
       FE_P2,
       FE_Q1,
       FV,
       DG
    };

    int nb_fields, pde, Sdm, ieq;
    string eq, nls, spD;
    Iteration ls;
    Preconditioner prec;
    vector<string> analytic;
    vector<int> field;
    vector<string> fn;
    vector<int> nb_dof;
    Equa<double> *theEquation;
    void setFields();
    int set(string e, Mesh* ms);
    int setEq();
    void set();
    void set(data *d);
    void setCoef();
    int setIn();
    int setBC();
    int setBF();
    int setSF();
    void check();
    void set(cmd* cmd) { _cmd = cmd; }
    void setNodeBC(int code, string exp, double t, Vect<double>& v);
    void setSize(Vect<double>& v, dataSize s);
    Log log;
    bool set_u, set_bc, set_bf, set_sf, set_in;
    Vect<double> u, b, bc, bf, sf, *theSolution[5];
    std::map<int,string> regex_bc, regex_sf;
    string regex_bf, regex_u;

 private:
    rita *_rita;
    int _verb, _dim, _ret, _nb_dof, _nb_fields;
    cmd *_cmd;
    bool _rho_set, _Cp_set, _kappa_set, _mu_set, _sigma_set, _Mu_set, _epsilon_set, _omega_set;
    bool _beta_set, _v_set, _young_set, _poisson_set;
    map<string,int> pde_map = {{"laplace",LAPLACE},
                               {"heat",HEAT},
                               {"wave",WAVE},
                               {"transport",TRANSPORT},
                               {"linear-elasticity",LINEAR_ELASTICITY},
                               {"truss",TRUSS},
                               {"beam",BEAM},
                               {"incompressible-navier-stokes",INCOMPRESSIBLE_NAVIER_STOKES}, 
                               {"compressible-euler",COMPRESSIBLE_EULER},
                               {"incompressible-porous-1phase",INCOMPRESSIBLE_POROUS_1PHASE}};
    const vector<string> _kw = {"expression","value","file","save"};
    Mesh *_theMesh;
    string _rho_exp, _Cp_exp, _kappa_exp, _mu_exp,_sigma_exp, _Mu_exp, _epsilon_exp, _omega_exp;
    string _beta_exp, _v_exp, _young_exp, _poisson_exp;
    OFELI::Fct _theFct;
};

} /* namespace RITA */

#endif
