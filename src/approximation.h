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

#ifndef __APPROXIMATION_H
#define __APPROXIMATION_H

#include "OFELI.h"
#include "OFELI.h"
#include "rita.h"
#include "cmd.h"
#include <map>

namespace RITA {
 
class approximation
{

 public:

    enum ApproxType {
       LAGRANGE,
       FITTING,
       BSPLINE,
       BEZIER,
       NURBS
    };
    approximation(rita* r, cmd* command, configure* config);
    ~approximation();
    int set();
    void set(std::ofstream* ofl, std::ofstream* ofh) { _ofl = ofl; _ofh = ofh; }
    int run();
    int go();

 private:

    rita *_rita;
    configure *_configure;
    cmd *_cmd;
    std::ofstream *_ofh, *_ofl;
    const vector<string> _kw = {"help","?","set","file","lagrange","fitting","bspline","bezier","nurbs"
                                "end","<","quit","exit","EXIT"};
    map<ApproxType,string> rApp = {{LAGRANGE,"lagrange"},
                                   {FITTING,"fitting"},
                                   {BSPLINE,"bspline"},
                                   {BEZIER,"bezier"},
                                   {NURBS,"nurbs"}};
};

} /* namespace RITA */

#endif
