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

                        Definition of class 'data'

  ==============================================================================*/

#ifndef __DATA_H
#define __DATA_H

#include <string>
#include <iostream>
using std::string;

#include "mesh/Mesh.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/Matrix.h"
#include "io/Fct.h"
#include "io/Tabulation.h"

namespace RITA {

class rita;
class configure;
class cmd;

enum eqType {
   ALGEBRAIC_EQ,
   ODE_EQ,
   PDE_EQ,
   OPT
};

enum dataSize {
   GIVEN_SIZE,
   GRID,
   NODES,
   ELEMENTS,
   SIDES,
   EDGES
};

class data
{

 public:

    data(rita *r, cmd *command, configure *config);
    ~data();
    int addField(string name, dataSize s, int n=0, int nb_dof=1);
    int checkField(string name);
    int checkFct(string name);
    int checkMesh(string name);
    int checkGrid(string name);
    int run();
    void setVerbose(int verb) { _verb = verb; }
    void setSave(int s) { _sr = s; }
    int getVerbose() const { return _verb; }
    void set(cmd* command) { _cmd = command; }
    int ret() const { return _ret; }
    int addFunction(const string &name, const string &def, const vector<string> &var);
    int addMesh(OFELI::Mesh* ms, string name);

    int getNbFields() const { return _nb_fields; }
    int getNbFcts() const { return _nb_fcts; }
    int getNbTabs() const { return _nb_tabs; }
    int getNbMeshes() const { return _nb_meshes; }
    vector<int> nb_dof;
    vector<OFELI::Vect<double> *> u;
    vector<string> Field;
    map<string,int> FieldName;
    vector<eqType> FieldType;
    vector<dataSize> FieldSizeType;
    void setNodeBC(int code, string exp, double t, OFELI::Vect<double>& v);
    std::vector<int> FieldEquation;
    vector<OFELI::Grid *> theGrid;
    vector<OFELI::Mesh *> theMesh;
    vector<OFELI::Tabulation *> theTab;
    vector<OFELI::Fct *> theFct;
    vector<OFELI::Matrix<double> *> theMatrix;
    vector<OFELI::Vect<double> *> theVector;
    vector<string> grid_name, tab_name, mesh_name, fct_name, vector_name, matrix_name;
    double obj, integral;
    bool ok;
    int iifct;

 private:

    rita *_rita;
    int _size, _ifield, _imesh, _igrid, _ifct, _itab;
    int _nb_fields, _default_field, _nb_fcts, _nb_meshes, _nb_tabs, _nb_grids, _nb_vectors, _nb_matrices; 
    int _nb_args, _verb, _key, _ret, _nb_dof, _nb, _sr;
    configure *_configure;
    cmd *_cmd;
    int _theMesh_alloc, _theTab_alloc, _theGrid_alloc, _theFct_alloc, _theVector_alloc, _theMatrix_alloc, _u_alloc;

    vector<string> _kw;
    OFELI::Mesh *_theMesh;
    OFELI::Tabulation *_theTab;
    OFELI::Grid *_theGrid;
    OFELI::Fct *_theFct;
    OFELI::Vect<double> *_theVector;
    OFELI::Matrix<double> *_theMatrix;
    OFELI::Vect<double> *_u;
    void getHelp();
    int setNbDOF();

    int setConfigure();
    int setVector();
    int setMatrix();
    int setGrid();
    int setField();
    int setTab();
    int setFunction();
    int setDerivative();
    //    void Clear();
    void Summary();
    vector<double> _xv;

    vector<string> _kw_data = {"help","?","set","grid","mesh","field","tab$ulation","func$tion",
                               "vect$or","matr$ix","clear","summary","end","<","quit","exit","EXIT"};
    vector<string> _kw_fct = {"help","?","set","name","var","exp","clear","end","<","quit","exit","EXIT"};
    vector<string> _kw_field = {"help","?","set","name","size","type","clear","end","<","quit","exit","EXIT"};
};

} /* namespace RITA */

#endif
