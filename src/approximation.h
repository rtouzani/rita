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

                       Definition of class 'approximation'

  ==============================================================================*/

#pragma once

#include "OFELI.h"
#include "io/Tabulation.h"
#include "rita.h"
#include "cmd.h"
#include <map>

namespace RITA {
 
class approximation
{

 public:

    enum ApproxType {
       LAGRANGE,
       PIECEWISE_LAGRANGE,
       HERMITE,
       FITTING,
       BSPLINE,
       BEZIER,
       NURBS
    };

    enum FitType {
       POLYNOMIAL,
       EXPONENTIAL,
       DEFINED
    };

    approximation(rita* r, cmd* command, configure* config);
    ~approximation();
    int set();
    int run();
    int go();

 private:

    rita *_rita;
    configure *_configure;
    cmd *_cmd;
    ApproxType _method;
    FitType ft;
    int _lagrange_degree, _hermite_degree;
    OFELI::Tabulation _tab;
    void lagrange();
    const vector<string> _kw = {"help","?","set","file","lagrange","fitting","bspline","bezier","nurbs"
                                "end","<","quit","exit","EXIT"};
    map<ApproxType,string> rApp = {{LAGRANGE,"lagrange"},
                                   {FITTING,"fitting"},
                                   {BSPLINE,"bspline"},
                                   {BEZIER,"bezier"},
                                   {NURBS,"nurbs"}};
};

} /* namespace RITA */
