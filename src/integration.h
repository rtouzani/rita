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

                          Definition of class 'integration'

  ==============================================================================*/

#pragma once

#include "rita.h"
#include "cmd.h"
#include "io/Fct.h"
#include <map>

namespace RITA {
 
class integration
{

 public:

    enum integration_formula {
       LRECTANGLE,
       RRECTANGLE,
       MIDPOINT,
       TRAPEZOIDAL,
       SIMPSON,
       GAUSS_LEGENDRE,
       GAUSS_LOBATTO
    };


    integration(rita *r, cmd* command, configure* config);
    ~integration();
    int set();
    int run();
    int go();
    integration_formula nim;
    vector<string> var;
    double xmin, xmax, ymin, ymax, zmin, zmax, res;
    int dim, nx, ny, nz, unif, ng;

 private:

    rita *_rita;
    std::ofstream *_ofh;
    configure *_configure;
    cmd *_cmd;
    vector<double> _x, _y, _z;
    OFELI::Fct *IFct;
    map<string,integration_formula> Nint = {{"left-rectangle",LRECTANGLE},
                                            {"right-rectangle",RRECTANGLE},
                                            {"mid-point",MIDPOINT},
                                            {"trapezoidal",TRAPEZOIDAL},
                                            {"simpson",SIMPSON},
                                            {"gauss-legendre",GAUSS_LEGENDRE},
                                            {"gauss-lobatto",GAUSS_LOBATTO}};
 };

} /* namespace RITA */
