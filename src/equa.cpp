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

                       Implementation of class 'equa'

  ==============================================================================*/

#include "equa.h"
#include "cmd.h"
#include "rita.h"

namespace RITA {

equa::equa(rita *r)
     : eq("laplace"), nls(""), spD("feP1"),
       ls(CG_SOLVER), prec(DILU_PREC), _nb_fields(0), _theMesh(nullptr)
{
   _rita = r;
   _verb = 0;
   field.resize(6);
   fn.resize(6);
   ieq = pde_map["laplace"];
   _rho_set = _Cp_set = _kappa_set = _mu_set = _sigma_set = _Mu_set = false;
   _epsilon_set = _omega_set = _beta_set = _v_set = _young_set = _poisson_set = false;
   set_u = set_bc = set_bf = set_sf = set_in = false;
}


equa::~equa()
{
   if (theEquation!=nullptr)
      delete theEquation;
}


int equa::set(string e, Mesh* ms)
{
   ieq = pde_map[e];
   eq = e;
   _theMesh = ms;
   _dim = _theMesh->getDim();
   setFields();
   _nb_dof = _theMesh->getNbDOF()/_theMesh->getNbNodes();
   return setEq();
}


void equa::setFields()
{
   switch (ieq) {

      case LAPLACE:
         nb_fields = 1;
         fn = {"u"};
         nb_dof.push_back(1);
         break;

      case HEAT:
         nb_fields = 1;
         fn = {"u"};
         nb_dof.push_back(1);
         break;

      case WAVE:
         nb_fields = 1;
         fn = {"u"};
         nb_dof.push_back(1);
         break;

      case TRANSPORT:
         nb_fields = 1;
         fn = {"u"};
         nb_dof.push_back(1);
         break;

      case LINEAR_ELASTICITY:
	 nb_fields = 1;
	 fn = {"u"};
         nb_dof.push_back(_dim);
         break;

      case TRUSS:
         nb_fields = 1;
         fn = {"u"};
         nb_dof.push_back(2);
         break;

      case BEAM:
         nb_fields = 1;
         fn = {"u"};
         if (_dim==2)
            nb_dof.push_back(2);
         else if (_dim==3)
            nb_dof.push_back(6);
         break;

      case INCOMPRESSIBLE_NAVIER_STOKES:
         nb_fields = 2;
         fn = {"v","p"};
         nb_dof.push_back(_dim);
         nb_dof.push_back(1);
         break;
         
      case COMPRESSIBLE_EULER:
         nb_fields = 4;
         fn = {"u","rho","p","e"};
         nb_dof.push_back(_dim);
         nb_dof.push_back(1);
         nb_dof.push_back(1);
         nb_dof.push_back(1);
         break;
         
      case COMPRESSIBLE_NAVIER_STOKES:
         nb_fields = 4;
         fn = {"u","rho","p","e"};
         nb_dof.push_back(_dim);
         nb_dof.push_back(1);
         nb_dof.push_back(1);
         nb_dof.push_back(1);
         break;

      case INCOMPRESSIBLE_POROUS_1PHASE:
         nb_fields = 1;
         fn = {"p"};
         nb_dof.push_back(1);
         break;

      case EDDY_CURRENTS:
         nb_fields = 1;
         fn = {"A"};
         nb_dof.push_back(1);
         break;

      case MAXWELL:
         nb_fields = 1;
         fn = {"u"};
         nb_dof.push_back(_dim);
         break;

      case HELMHOLTZ:
         nb_fields = 1;
	 fn = {"p"};
         nb_dof.push_back(1);
         break;

      default:
         nb_fields = 1;
         fn = {"u"};
         nb_dof.push_back(1);
         break;
   }
}


void equa::setSize(Vect<double>& v, dataSize s)
{
   _theMesh = _rita->_theMesh;
   if (_theMesh==nullptr) {
      cout << "Error: No mesh data available." << endl;
      *_ofl << "In rita>pde>: No mesh data available." << endl;
      _ret = 1;
      return;
   }
   _nb_dof = 1;
   if (_theMesh->getDOFSupport()==NODE_DOF)
      _nb_dof = _theMesh->getNbDOF()/_theMesh->getNbNodes();
   else if (_theMesh->getDOFSupport()==SIDE_DOF)
      _nb_dof = _theMesh->getNbDOF()/_theMesh->getNbSides();
   else if (_theMesh->getDOFSupport()==ELEMENT_DOF)
      _nb_dof = _theMesh->getNbDOF()/_theMesh->getNbElements();
   else if (_theMesh->getDOFSupport()==EDGE_DOF)
      _nb_dof = _theMesh->getNbDOF()/_theMesh->getNbEdges();

   if (s==NODES) {
      if (_theMesh->getNbNodes()==0) {
         cout << "Error: Mesh has no nodes." << endl;
         *_ofl << "In rita>pde>: Mesh has no nodes" << endl;
         _ret = 1;
         return;
      }
      v.setMesh(*_theMesh,NODE_DOF,_nb_dof);
   }
   else if (s==ELEMENTS) {
      if (_theMesh->getNbElements()==0) {
         cout << "Error: Mesh has no elements." << endl;
         *_ofl << "In rita>pde>: Mesh has no elements" << endl;
         _ret = 1;
         return;
      }
      v.setMesh(*_theMesh,ELEMENT_DOF,_nb_dof);
   }
   else if (s==SIDES) {
      if (_theMesh->getNbSides()==0) {
         cout << "Error: Mesh has no sides." << endl;
         *_ofl << "In rita>pde>: Mesh has no sides" << endl;
         _ret = 1;
         return;
      }
      v.setMesh(*_theMesh,SIDE_DOF,_nb_dof);
   }
   else if (s==EDGES) {
      if (_theMesh->getNbEdges()==0) {
         cout << "Error: Mesh has no edges." << endl;
         *_ofl << "In rita>pde>: Mesh has no edges" << endl;
         _ret = 1;
         return;
      }
      v.setMesh(*_theMesh,EDGE_DOF,_nb_dof);
   }
}


int equa::setIn()
{
   u.setMesh(*_theMesh,NODE_DOF,_nb_dof);
   theSolution[0] = &u;
   bool val_ok=false, file_ok=false, save_ok=false;
   _ret = 0;
   string file="rita-init.dat", save="rita-init.out";
   vector<string> kw = {"val$ue","file","save"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args==0) {
      cout << "Error in command: No given arguments" << endl;
      *_ofl << "In rita>pde>initial>: Error in command: No arguments" << endl;
      return _ret;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg("=");
      switch (n) {

         case 0:
            regex_u = _cmd->string_token();
            val_ok = true;
            break;

         case 1:
            file = _cmd->string_token();
            file_ok = true;
            break;

         case 2:
            save = _cmd->string_token();
            save_ok = true;
            break;

         default:
            cout << "Error: Unknown argument: " << _kw[n] << endl;
            *_ofl << "In rita>pde>initial>: Unknown argument: " << _kw[n] << endl;
            return _ret;
      }
   }
   if (!val_ok) {
      cout << "Error: No value or expression given for initial condition." << endl;
      *_ofl << "In rita>pde>initial> No value or expression given for initial condition." << endl;
   }
   *_ofh << "  in  value=" << regex_u;
   u.setRegex(1);
   if (file_ok) {
      OFELI::IOField ffi(file,OFELI::IOField::IN);
      ffi.get(u);
      *_ofh << "  file=" << file;
   }
   if (save_ok) {
      OFELI::IOField ffo(save,OFELI::IOField::OUT);
      ffo.put(u);
      if (_verb)
         cout << "Initial condition saved in file: " << save << endl;
      *_ofh << "  save=" << save;
   }
   set_u = true;
   *_ofh << endl;
   return 0;
}


int equa::setBC()
{
   setSize(bc,NODES);
   _ret = 0;
   string file="rita-bc.dat", save="rita-bc.out", val="";
   bool code_ok=false, val_ok=false, file_ok=false, save_ok=false;
   int code=0;
   vector<string> kw = {"code","val$ue","file","save"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args==0) {
      cout << "Error in command: No given arguments" << endl;
      *_ofl << "In rita>pde>bc>: Error in command: No arguments" << endl;
      return _ret;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg("=");
      switch (n) {

         case 0:
            code = _cmd->int_token();
            code_ok = true;
            break;

         case 1:
            val = _cmd->string_token();
            val_ok = true;
            break;

         case 2:
            file = _cmd->string_token();
            file_ok = true;
            break;

         case 3:
            save = _cmd->string_token();
            save_ok = true;
            break;

         default:
            cout << "Error: Unknown argument: " << kw[n] << endl;
            *_ofl << "In rita>pde>bc>: Unknown argument: " << kw[n] << endl;
            return 1;
      }
   }
   if (!code_ok) {
      cout << "Error: No code given for boundary condition." << endl;
      *_ofl << "In rita>pde>bc> No code given for boundary condition." << endl;
      return 1;
   }
   if (!val_ok) {
      cout << "Error: No value or expression given for boundary condition." << endl;
      *_ofl << "In rita>pde>bc> No value or expression given for boundary condition." << endl;
      return 1;
   }
   if (code<=0) {
      cout << "Error: You cannot give a boundary condition for a nonpositive code." << endl;
      *_ofl << "In rita>pde>bc>: Illegal value of code: " << code << endl;
      return 1;
   }
   regex_bc.insert(pair<int,string>(code,val));
   bc.setRegex(1);
   set_bc = true;
   if (_verb)
      cout << "Nodes with code " << code << " have prescribed value by the expression: " << val << endl;
   *_ofh << "  bc  code=" << code << "  value=" << val;
   if (file_ok) {
      OFELI::IOField ffi(file,OFELI::IOField::IN);
      ffi.get(bc);
      *_ofh << "  file=" << file;
   }
   if (save_ok) {
      OFELI::IOField ffo(file,OFELI::IOField::OUT);
      ffo.put(bc);
      *_ofh << "  save=" << save;
      if (_verb)
         cout << "Boundary condition saved in file: " << save << endl;
   }
   *_ofh << endl;
   return 0;
}


int equa::setSF()
{
   setSize(sf,SIDES);
   int code=0;
   _ret = 0;
   string val="", file="rita-sf.dat", save="rita-sf.out";
   bool code_ok=false, val_ok=false, file_ok=false, save_ok=false;
   vector<string> kw = {"code","val$ue","file","save"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args==0) {
      cout << "Error in command: No given arguments" << endl;
      *_ofl << "In rita>pde>sf>: Error in command: No arguments" << endl;
      return _ret;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg("=");
      switch (n) {

         case 0:
            code = _cmd->int_token();
            code_ok = true;
            break;

         case 1:
            val = _cmd->string_token();
            val_ok = true;
            break;

         case 2:
            file = _cmd->string_token();
            file_ok = true;
            break;

         case 3:
            save = _cmd->string_token();
            save_ok = true;
            break;

         default:
            cout << "Error: Unknown argument: " << kw[n] << endl;
            *_ofl << "In rita>pde>sf>: Unknown argument: " << kw[n] << endl;
            return 1;
      }
   }
   if (!code_ok) {
      cout << "Error: No code given." << endl;
      *_ofl << "In rita>pde>sf> No code given." << endl;
      return 1;
   }
   if (!val_ok) {
      cout << "Error: No value or expression given for surface force." << endl;
      *_ofl << "In rita>pde>sf> No value or expression given for surface force." << endl;
      return 1;
   }
   if (code<=0) {
      cout << "Error: You cannot give a surface force for a nonpositive code." << endl;
      *_ofl << "In rita>pde>sf>: Illegal value of code: " << code << endl;
      return 1;
   }
   *_ofh << "  sf  code=" << code << "  value=" << val;
   regex_sf.insert(std::pair<int,string>(code,val));
   sf.setRegex(1);
   set_sf = true;
   if (file_ok) {
      OFELI::IOField ffi(file,OFELI::IOField::IN);
      ffi.get(sf);
      *_ofh << "  file=" << file;
   }
   if (save_ok) {
      OFELI::IOField ffo(save,OFELI::IOField::OUT);
      ffo.put(sf);
      *_ofh << "  save=" << save;
   }
   *_ofh << endl;
   return 0;
}


int equa::setBF()
{
   setSize(bf,NODES);
   _ret = 0;
   bool val_ok=false, file_ok=false, save_ok=false;
   string file="rita-source.dat", save;
   vector<string> kw = {"val$ue","file","save"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args==0) {
      cout << "Error in command: No given arguments" << endl;
      *_ofl << "In rita>pde>bf>: Error in command: No arguments" << endl;
      return _ret;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg("=");
      switch (n) {

         case 0:
            regex_bf = _cmd->string_token();
            val_ok = true;
            break;

         case 1:
            file = _cmd->string_token();
            {
               OFELI::IOField ffi(file,OFELI::IOField::IN);
               ffi.get(bf);
            }
            file_ok = true;
            break;

         case 2:
            save = _cmd->string_token();
            {
               OFELI::IOField ffo(file,OFELI::IOField::OUT);
               ffo.put(bf);
            }
            save_ok = true;
            break;

         default:
            cout << "Error: Unknown argument: " << kw[n] << endl;
            *_ofl << "In rita>pde>source>: Unknown argument: " << kw[n] << endl;
	    return 1;
      }
   }
   if (!val_ok) {
      cout << "Error: No value or expression given for source." << endl;
      *_ofl << "In rita>pde>source> No value or expression given for source." << endl;
   }
   *_ofh << "  source  value=" << regex_bf;
   bf.setRegex(1);
   set_bf = true;
   if (file_ok)
      *_ofh << "  file=" << file;
   if (save_ok)
      *_ofh << "  save=" << save;
   *_ofh << endl;
   return 0;
}


void equa::setNodeBC(int code, string exp, double t, Vect<double>& v)
{
   const vector<string> var = {"x","y","z","t"};
   _theFct.set(exp,var);
   for (size_t n=1; n<=_theMesh->getNbNodes(); ++n) {
      Node *nd = (*_theMesh)[n];
      for (size_t i=1; i<=nd->getNbDOF(); ++i) {
         if (nd->getCode(i)==code)
            v(nd->n(),i) = _theFct(nd->getCoord());
      }
   }
}


int equa::setEq()
{
   int ret = 0;
   if (_theMesh==nullptr) {
      cout << "Error: No mesh provided." << endl;
      *_ofl << "In rita>pde>: No mesh provided." << endl;
      log.mesh = true;
      return 1;
   }
   setFields();

   switch (ieq) {

//    Laplace equation
      case LAPLACE:
         switch (_dim) {
  
            case 1:
               if (spD!="feP1" && spD!="feP2") {
                  cout << "Error: This approximation of the equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Approximation of equation not implemented in rita." << endl;
                  log.pde = true;
               }
               break;

            case 2:
               if (spD!="feP1") {
                  cout << "Error: This approximation of the equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Approximation of equation not implemented in rita." << endl;
                  log.spd = true;
               }
               break;

            case 3:
               if (spD!="feP1") {
                  cout << "Error: This approximation of the equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Approximation of equation not implemented in rita." << endl;
                  log.spd = true;
               }
               break;
         }
         break;

//    Heat equation
      case HEAT:
         switch (_dim) {

            case 1:
               if (spD!="feP1") {
                  cout << "Error: This approximation of the equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Approximation of equation not implemented in rita." << endl;
                  log.pde = true;
               }
               break;

            case 2:
               if (spD!="feP1" && spD!="feP2") {
                  cout << "Error: This approximation of the equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Approximation of equation not implemented in rita." << endl;
                  log.spd = true;
               }
               break;

            case 3:
               if (spD!="feP1") {
                  cout << "Error: This approximation of the equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Approximation of equation not implemented in rita." << endl;
                  log.spd = true;
               }
                    break;
         }
         break;

//    Wave equation
      case WAVE:
         switch (_dim) {

            case 1:
               if (spD=="feP1") {
                  cout << "Error: This equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Equation not implemented in rita." << endl;
                  log.pde = true;
               }
               else {
                  cout << "Error: This approximation of the equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Approximation of equation not implemented in rita." << endl;
                  log.spd = true;
               }
               break;

            case 2:
               if (spD=="feP1") {
                  cout << "Error: This equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Equation not implemented in rita." << endl;
                  log.pde = true;
               }
               else {
                  cout << "Error: This approximation of the equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Approximation of equation not implemented in rita." << endl;
                  log.spd = true;
               }
               break;

            case 3:
               if (spD=="feP1") {
                  cout << "Error: This equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Equation not implemented in rita." << endl;
                  log.pde = true;
               }
               else {
                  cout << "Error: This approximation of the equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Approximation of equation not implemented in rita." << endl;
                  log.spd = true;
               }
               break;
         }
         break;

//    Linear transport equation
      case TRANSPORT:
         switch (_dim) {

            case 1:
               if (spD=="feP1") {
                  cout << "Error: This equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Equation not implemented in rita." << endl;
                  log.pde = true;
               }
               else {
                  cout << "Error: This approximation of the equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Approximation of equation not implemented in rita." << endl;
                  log.spd = true;
               }
               break;

            case 2:
               if (spD=="feP1") {
                  cout << "Error: This equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Equation not implemented in rita." << endl;
                  log.pde = true;
               }
               else {
                  cout << "Error: This approximation of the equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Approximation of equation not implemented in rita." << endl;
                  log.spd = true;
               }
               break;

            case 3:
               if (spD=="feP1") {
                  cout << "Error: This equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Equation not implemented in rita." << endl;
                  log.pde = true;
               }
               else {
                  cout << "Error: This approximation of the equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Approximation of equation not implemented in rita." << endl;
                  log.spd = true;
               }
               break;
         }
         break;

      case LINEAR_ELASTICITY:
         switch (_dim) {

            case 1:
               if (spD=="feP1") {
                  cout << "Error: This equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Equation not implemented in rita." << endl;
                  log.pde = true;
               }
               else {
                  cout << "Error: This approximation of the equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Approximation of equation not implemented in rita." << endl;
                  log.spd = true;
               }
               break;

            case 2:
               if (spD!="feP1" && spD!="feQ1") {
                  cout << "Error: This approximation of the equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Approximation of equation not implemented in rita." << endl;
                  log.pde = true;
               }
               break;

            case 3:
               if (spD!="feP1" && spD!="feQ1") {
                  cout << "Error: This approximation of the equation is not implemented in rita." << endl;
                  *_ofl << "In rita>pde>: Approximation of equation not implemented in rita." << endl;
                  log.spd = true;
               }
               break;
         }
         break;

      case TRUSS:
         switch (_dim) {

            case 1:
               cout << "Error: No implementation for 1-D available for bar equation." << endl;
               *_ofl << "In rita>pde>: No implementation for 1-D available for bar equation." << endl;
               log.pde = true;
               break;

            case 2:
               if (spD!="feP1") {
                  cout << "Error: Only 2-D P1 finite element is available for bar equation." << endl;
                  *_ofl << "In rita>pde>: Only 2-D P1 finite element is available for bar equation." << endl;
                  log.pde = true;
               }
               break;

            case 3:
               cout << "Error: No implementation for 3-D available for bar equation." << endl;
               *_ofl << "In rita>pde>: No implementation for 3-D available for bar equation." << endl;
               log.pde = true;
               break;
         }
         break;
         
      case BEAM:
         switch (_dim) {

            case 1:
               if (spD=="feP1") {
                  cout << "Error: No implementation for 1-D P1 finite element available for beam equation." << endl;
                  *_ofl << "In rita>pde>: No implementation for 1-D P1 finite element available for beam equation." << endl;
                  log.pde = true;
               }
               break;

            case 2:
               if (spD=="feP1") {
                  cout << "Error: No implementation for 2-D P1 finite element available for beam equation." << endl;
                  *_ofl << "In rita>pde>: No implementation for 2-D P1 finite element available for beam equation." << endl;
                  log.pde = true;
               }
               else if (spD=="feP2") {
                  cout << "Error: No implementation for 2-D P2 finite element available for beam equation." << endl;
                  *_ofl << "In rita>pde>: No implementation for 2-D P2 finite element available for beam equation." << endl;
                   log.pde = true;
               }
               break;

            case 3:
               if (spD=="feP1") {
                  cout << "Error: No implementation for 3-D P1 finite element available for beam equation." << endl;
                  *_ofl << "In rita>pde>: No implementation for 3-D P1 finite element available for beam equation." << endl;
                  log.pde = true;
               }
               break;
         }
         break;

      case INCOMPRESSIBLE_NAVIER_STOKES:
         switch (_dim) {

            case 1:
               if (spD=="feP1") {
                  cout << "Error: No implementation for 1-D available for incompressible Navier-Stokes equations." << endl;
                  *_ofl << "In rita>pde>: No implementation for 1-D available for incompressible Navier-Stokes equations." << endl;
                  log.pde = true;
               }
               break;

            case 2:
               if (spD!="feP1") {
                  cout << "Error: Only P1 finite element is implemented is available for incompressible Navier-Stokes equations." << endl;
                  *_ofl << "In rita>pde>: Only P1 finite element is implemented is available for incompressible Navier-Stokes equations." << endl;
                  log.pde = true;
               }
               break;

            case 3:
               if (spD=="feP1")
                  log.pde = true;
               break;
         }
         break;
         
      case COMPRESSIBLE_EULER:
         switch (_dim) {

            case 1:
               if (spD=="feP1") {
                  cout << "Error: This Space discretization method is not implemented for this PDE." << endl;
                  *_ofl << "In rita>pde>: Space discretization method not implemented for this PDE." << endl;
                   log.pde = true;
               }
               break;

            case 2:
               if (spD=="feP1") {
                  cout << "Error: This Space discretization method is not implemented for this PDE." << endl;
                  *_ofl << "In rita>pde>: Space discretization method not implemented for this PDE." << endl;
                  log.pde = true;
               }
               else if (spD=="feP2") {
                  cout << "Error: This Space discretization method is not implemented for this PDE." << endl;
                  *_ofl << "In rita>pde>: Space discretization method not implemented for this PDE." << endl;
                  log.pde = true;
               }
               break;

            case 3:
               if (spD=="feP1") {
                  cout << "Error: This Space discretization method is not implemented for this PDE." << endl;
                  *_ofl << "In rita>pde>: Space discretization method not implemented for this PDE." << endl;
                  log.pde = true;
               }
               break;
         }
         break;

      case INCOMPRESSIBLE_POROUS_1PHASE:
         switch (_dim) {

            case 1:
               if (spD=="feP1") {
                  cout << "Error: This Space discretization method is not implemented for this PDE." << endl;
                  *_ofl << "In rita>pde>: Space discretization method not implemented for this PDE." << endl;
                  log.pde = true;
               }
               break;

            case 2:
               if (spD=="feP1") {
                  cout << "Error: This Space discretization method is not implemented for this PDE." << endl;
                  *_ofl << "In rita>pde>: Space discretization method not implemented for this PDE." << endl;
                  log.pde = true;
               }
               else if (spD=="feP2") {
                  cout << "Error: This Space discretization method is not implemented for this PDE." << endl;
                  *_ofl << "In rita>pde>: Space discretization method not implemented for this PDE." << endl;
                  log.pde = true;
               }
               break;

            case 3:
               if (spD=="feP1") {
                  cout << "Error: This Space discretization method is not implemented for this PDE." << endl;
                  *_ofl << "In rita>pde>: Space discretization method not implemented for this PDE." << endl;
                  log.pde = true;
               }
               break;
         }
         break;

      default:
         cout << "Error: This equation is not implemented in rita." << endl;
         *_ofl << "In rita>pde>: Equation not implemented in rita." << endl;
         log.pde = true;
         break;

   }
   return ret;
}


void equa::setCoef()
{
   int ret = 0, key = 0;
   static const string H = "Command: coef [rho=x] [Cp=x] [kappa=x] [Mu=x] [sigma=x] [mu=x] [epsilon=x] [omega=x]\n"
                           "              [beta=x] [v=x] [young=x] [poisson=x]\n\n";
   const static vector<string> kw = {"help","?","set","rho","density","Cp","specific-heat","kappa","thermal-conductivity",
                                     "Mu","magnetic-permeability","sigma","electric-conductivity","mu","viscosity","epsilon",
                                     "electric-permittivity","omega","angular-frequency","beta","thermal-dilatation","v",
                                     "velocity","young","poisson","<","quit","exit","EXIT"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args==0) {
      cout << "Error: No argument for command." << H << endl;
      *_ofl << "In rita>pde>coef>: No argument for command." << endl;
      _ret = 1;
      return;
   }
   if (nb_args<1) {
      cout << "Error: No argument for command." << endl;
      *_ofl << "In rita>pde>coef>: No argument for command." << endl;
      _ret = 1;
      return;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg();
      switch (n) {

         case  3:
         case  4:
            _rho_exp = _cmd->string_token();
            _rho_set = true;
            break;

         case  5:
         case  6:
            _Cp_exp = _cmd->string_token();
            _Cp_set = true;
            break;

         case  7:
         case  8:
            _kappa_exp = _cmd->string_token();
            _kappa_set = true;
            break;

         case  9:
         case 10:
            _Mu_exp = _cmd->string_token();
            _Mu_set = true;
            break;

         case 11:
         case 12:
            _sigma_exp = _cmd->string_token();
            _sigma_set = true;
            break;

         case 13:
         case 14:
            _mu_exp = _cmd->string_token();
            _mu_set = true;
            break;

         case 15:
         case 16:
            _epsilon_exp = _cmd->string_token();
            _epsilon_set = true;
            break;

         case 17:
         case 18:
            _omega_exp = _cmd->string_token();
            _omega_set = true;
            break;

         case 19:
         case 20:
            _beta_exp = _cmd->string_token();
            _beta_set = true;
            break;

         case 21:
         case 22:
            _v_exp = _cmd->string_token();
            _v_set = true;
            break;

         case 23:
            _young_exp = _cmd->string_token();
            _young_set = true;
            break;

         case 24:
            _poisson_exp = _cmd->string_token();
            _poisson_set = true;
            break;

         default:
            cout << "Error: Unknown argument: " << _kw[n] << endl;
            *_ofl << "In rita>pde>coef>: Unknown argument: " << _kw[n] << endl;
            _ret = 1;
	    return;
      }
   }
   if (nb_args>0) {
      if (_rho_set) {
         if (ieq!=HEAT && ieq!=INCOMPRESSIBLE_NAVIER_STOKES && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
            cout << "Error: This PDE doesn't need density input" << endl;
            *_ofl << "In rita>pde>coef>: This PDE doesn't need density input" << endl;
            _ret = 1;
            return;
         }
      }
      if (_Cp_set) {
         if (ieq!=HEAT && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
            cout << "Error: This PDE doesn't need specific heat input" << endl;
            *_ofl << "In rita>pde>coef>: This PDE doesn't need specific heat input" << endl;
            _ret = 1;
            return;
         }
      }
      if (_kappa_set) {
         if (ieq!=HEAT && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
            cout << "Error: This PDE doesn't need thermal conductivity input" << endl;
            *_ofl << "In rita>pde>coef>: This PDE doesn't need thermal conductivity input" << endl;
            _ret = 1;
            return;
         }
      }
      if (_mu_set) {
         if (ieq!=INCOMPRESSIBLE_NAVIER_STOKES && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
            cout << "Error: This PDE doesn't need viscosity input" << endl;
            *_ofl << "In rita>pde>coef>: This PDE doesn't need viscosity input" << endl;
            _ret = 1;
            return;
         }
      }
      if (_sigma_set) {
         if (ieq!=EDDY_CURRENTS && ieq!=MAXWELL && ieq!=HELMHOLTZ) {
            cout << "Error: This PDE doesn't need electric conductivity input" << endl;
            *_ofl << "In rita>pde>coef>: This PDE doesn't need electric conductivity input" << endl;
            _ret = 1;
            return;
         }
      }
      if (_Mu_set) {
         if (ieq!=EDDY_CURRENTS && ieq!=MAXWELL && ieq!=HELMHOLTZ) {
            cout << "Error: This PDE doesn't need magnetic permeability input" << endl;
            *_ofl << "In rita>pde>coef>: This PDE doesn't need viscosity input" << endl;
            _ret = 1;
            return;
         }
      }
      if (_epsilon_set) {
         if (ieq!=MAXWELL && ieq!=HELMHOLTZ) {
            cout << "Error: This PDE doesn't need electric permittivity input" << endl;
            *_ofl << "In rita>pde>coef>: This PDE doesn't need electric permittivity input" << endl;
            _ret = 1;
            return;
         }
      }
      if (_omega_set) {
         if (ieq!=EDDY_CURRENTS && ieq!=MAXWELL && ieq!=HELMHOLTZ) {
            cout << "Error: This PDE doesn't need angular frequency input" << endl;
            *_ofl << "In rita>pde>coef>: This PDE doesn't need angular frequency input" << endl;
            _ret = 1;
            return;
         }
      }
      if (_beta_set) {
         if (ieq!=HEAT && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
            cout << "Error: This PDE doesn't need thermal expansion coefficient input" << endl;
            *_ofl << "In rita>pde>coef>: This PDE doesn't need thermal expansion coefficient input" << endl;
            _ret = 1;
            return;
         }
      }
      if (_v_set) {
         if (ieq!=WAVE && ieq!=TRANSPORT) {
            cout << "Error: This PDE doesn't need velocity input" << endl;
            *_ofl << "In rita>pde>coef>: This PDE doesn't need velocity input" << endl;
            _ret = 1;
            return;
         }
      }
      if (_young_set) {
         if (ieq!=LINEAR_ELASTICITY && ieq!=BEAM) {
            cout << "Error: This PDE doesn't need Young modulus input" << endl;
            *_ofl << "In rita>pde>coef>: This PDE doesn't need Young's modulus input" << endl;
            _ret = 1;
            return;
         }
      }
      if (_poisson_set) {
         if (ieq!=LINEAR_ELASTICITY && ieq!=BEAM) {
            cout << "Error: This PDE doesn't need Poisson ratio input" << endl;
            *_ofl << "In rita>pde>coef>: This PDE doesn't need Poisson ratio input" << endl;
            _ret = 1;
            return;
         }
      }
   }
   /*
   else {
      while (1) {
         if (_cmd->readline("rita>pde>coef> ")<0)
            continue;
         switch (key=_cmd->getKW(kw)) {

            case  0:
            case  1:
               _cmd->setNbArg(0);
               cout << "\nAvailable Commands:\n";
               cout << "rho:      Density\n";
               cout << "Cp:       Specific heat at constant pressure\n";
               cout << "kappa:    Thermal conductivity\n";
               cout << "mu:       Viscosity\n";
               cout << "sigma:    Electric conductivity\n";
               cout << "Mu:       Magnetic permeability\n";
               cout << "epsilon:  Electric permittivity\n";
               cout << "omega:    Angular frequency\n";
               cout << "beta:     Thermal dilatation coefficient\n";
               cout << "v:        Velocity\n";
               cout << "young:    Young modulus\n";
               cout << "poisson:  Poisson ratio\n";
               cout << "end or <: go back to higher level" << endl;
               break;

            case  2:
               _rita->setConfigure();
               break;

            case  3:
            case  4:
               if (ieq!=HEAT && ieq!=INCOMPRESSIBLE_NAVIER_STOKES && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
                  cout << "Error: This PDE doesn't need density input" << endl;
                  *_ofl << "In rita>pde>coef>rho>: This PDE doesn't need density input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for density as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>rho>: Missing regular expression for density" << endl;
                  break;
               }
               ret = _cmd->get(_rho_exp);
               if (!ret) {
                  *_ofh << "    rho " << _rho_exp << endl;
                  _rho_set = true;
	       }
               break;

            case  5:
            case  6:
               if (ieq!=HEAT && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
                  cout << "Error: This PDE doesn't need specific heat input" << endl;
                  *_ofl << "In rita>pde>coef>Cp>: This PDE doesn't need specific heat input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for specific heat as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>Cp>: Missing regular expression for specific heat" << endl;
                  break;
               }
               ret = _cmd->get(_Cp_exp);
               if (!ret)
                  *_ofh << "    Cp " << _Cp_exp << endl;
               break;

            case  7:
            case  8:
               if (ieq!=HEAT && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
                  cout << "Error: This PDE doesn't need thermal conductivity input" << endl;
                  *_ofl << "In rita>pde>coef>kappa>: This PDE doesn't need thermal conductivity input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for thermal conductivity as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>Cp>: Missing regular expression for thermal conductivity" << endl;
                  break;
               }
               ret = _cmd->get(_kappa_exp);
               if (!ret) {
                  *_ofh << "    kappa " << _kappa_exp << endl;
                  _kappa_set = true;
               }
               break;

            case  9:
            case 10:
               if (ieq!=INCOMPRESSIBLE_NAVIER_STOKES && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
                  cout << "Error: This PDE doesn't need viscosity input" << endl;
                  *_ofl << "In rita>pde>coef>mu>: This PDE doesn't need viscosity input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for viscosity as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>mu>: Missing regular expression for viscosity" << endl;
                  break;
               }
               ret = _cmd->get(_mu_exp);
               if (!ret) {
                  *_ofh << "    mu " << _mu_exp << endl;
                  _mu_set = true;
               }
               break;

            case 11:
            case 12:
               if (ieq!=EDDY_CURRENTS && ieq!=MAXWELL && ieq!=HELMHOLTZ) {
                  cout << "Error: This PDE doesn't need electric conductivity input" << endl;
                  *_ofl << "In rita>pde>coef>sigma>: This PDE doesn't need viscosity input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for electric conductivity as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>sigma>: Missing regular expression for electric conductivity" << endl;
                  break;
               }
               ret = _cmd->get(_sigma_exp);
               if (!ret) {
                  *_ofh << "    sigma " << _sigma_exp << endl;
                  _sigma_set = true;
               }
               break;

            case 13:
            case 14:
               if (ieq!=EDDY_CURRENTS && ieq!=MAXWELL && ieq!=HELMHOLTZ) {
                  cout << "Error: This PDE doesn't need magnetic permeability input" << endl;
                  *_ofl << "In rita>pde>coef>Mu>: This PDE doesn't need viscosity input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for magnetic permeability as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>Mu>: Missing regular expression for magnetic permeability" << endl;
                  break;
               }
               ret = _cmd->get(_Mu_exp);
               if (!ret) {
                  *_ofh << "    Mu " << _Mu_exp << endl;
                  _Mu_set = true;
               }
               break;

            case 15:
            case 16:
               if (ieq!=MAXWELL && ieq!=HELMHOLTZ) {
                  cout << "Error: This PDE doesn't need electric permittivity input" << endl;
                  *_ofl << "In rita>pde>coef>epsilon>: This PDE doesn't need electric permittivity input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for electric permittivity as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>epsilon>: Missing regular expression for electric permittivity" << endl;
                  break;
               }
               ret = _cmd->get(_epsilon_exp);
               if (!ret) {
                  *_ofh << "    epsilon " << _epsilon_exp << endl;
                  _epsilon_set = true;
               }
               break;

            case 17:
            case 18:
               if (ieq!=EDDY_CURRENTS && ieq!=MAXWELL && ieq!=HELMHOLTZ) {
                  cout << "Error: This PDE doesn't need angular frequency input" << endl;
                  *_ofl << "In rita>pde>coef>omega>: This PDE doesn't need angular frequency input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give value of angular frequency.")) {
                  *_ofl << "In rita>pde>coef>omega>: Missing value of angular frequency" << endl;
                  break;
               }
               ret = _cmd->get(_omega_exp);
               if (!ret) {
                  *_ofh << "    omega " << _omega_exp << endl;
                  _omega_set = true;
               }
               break;

            case 19:
            case 20:
               if (ieq!=HEAT && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
                  cout << "Error: This PDE doesn't need thermal expansion coefficient input" << endl;
                  *_ofl << "In rita>pde>coef>beta>: This PDE doesn't need thermal expansion coefficient input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for thermal expansion coefficient as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>beta>: Missing regular expression for thermal expansion coefficient" << endl;
                  break;
               }
               ret = _cmd->get(_beta_exp);
               if (!ret) {
                  *_ofh << "    beta " << _beta_exp << endl;
                  _beta_set = true;
               }
               break;

            case 21:
            case 22:
               if (ieq!=WAVE && ieq!=TRANSPORT) {
                  cout << "Error: This PDE doesn't need velocity input" << endl;
                  *_ofl << "In rita>pde>coef>v>: This PDE doesn't need velocity input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for velocity as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>v>: Missing regular expression for velocity" << endl;
                  break;
               }
               ret = _cmd->get(_v_exp);
               if (!ret) {
                  *_ofh << "    v " << _v_exp << endl;
                  _v_set = true;
               }
               break;

            case 23:
               if (ieq!=LINEAR_ELASTICITY && ieq!=BEAM) {
                  cout << "Error: This PDE doesn't need Young modulus input" << endl;
                  *_ofl << "In rita>pde>coef>young>: This PDE doesn't need Young's modulus input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for Young's modulus as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>young>: Missing regular expression for Young's modulus" << endl;
                  break;
               }
               ret = _cmd->get(_young_exp);
               if (!ret) {
                  *_ofh << "    young " << _young_exp << endl;
                  _young_set = true;
               }
               break;

            case 24:
               if (ieq!=LINEAR_ELASTICITY && ieq!=BEAM) {
                  cout << "Error: This PDE doesn't need Poisson ratio input" << endl;
                  *_ofl << "In rita>pde>coef>poisson>: This PDE doesn't need Poisson ratio input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for Poisson ratio as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>poisson>: Missing regular expression for Poisson ratio" << endl;
                  break;
               }
               ret = _cmd->get(_poisson_exp);
               if (!ret) {
                  *_ofh << "    poisson " << _poisson_exp << endl;
                  _poisson_set = true;
               }
               break;

            case 25:
            case 26:
               _ret = 0;
               *_ofh << "    end" << endl;
               return;

            case 27:
            case 28:
               _ret = 100;
               return;

            case 29:
               _ret = 200;
               return;

            case -2:
            case -3:
            case -4:
               break;

            default:
               cout << "Unknown Command: " << _cmd->token() << endl;
               cout << "Available commands: rho, Cp, kappa, mu, sigma, Mu, epsilon" << endl;
	       cout << "                    omega, beta, v, young, poisson, end, <" << endl;
               cout << "Global commands:    help, ?, set, quit, exit" << endl;
               *_ofl << "In rita>pde>coef>: Unknown PDE Coefficient " << _cmd->token() << endl;
               break;
         }
      }
      }*/
}


void equa::set()
{
   int ret = 0;
   switch (ieq) {

      case LAPLACE:
         switch (_dim) {
  
            case 1:
               if (spD=="feP1")
                  theEquation = new Laplace1DL2(*_theMesh);
               else if (spD=="feP2")
                  theEquation = new Laplace1DL3(*_theMesh);
               break;

            case 2:
               if (spD=="feP1")
                  theEquation = new Laplace2DT3(*_theMesh);
               break;

            case 3:
               if (spD=="feP1")
                  theEquation = new Laplace3DT4(*_theMesh);
               break;
         }
         break;

      case HEAT:

         switch (_dim) {

            case 1:
               if (spD=="feP1") {
                  theEquation = new DC1DL2(*_theMesh);
                  theEquation->setTerms(LUMPED_CAPACITY|DIFFUSION);
               }
               break;

            case 2:
               if (spD=="feP1")
                  theEquation = new DC2DT3(*_theMesh);
               else if (spD=="feP2")
                  theEquation = new DC2DT6(*_theMesh);
               theEquation->setTerms(LUMPED_CAPACITY|DIFFUSION);
               break;

            case 3:
               if (spD=="feP1") {
                  theEquation = new DC3DT4(*_theMesh);
                  theEquation->setTerms(LUMPED_CAPACITY|DIFFUSION);
	       }
               break;
         }
         if (_rho_set)
            theEquation->set_rho(_rho_exp);
         if (_Cp_set)
            theEquation->set_Cp(_Cp_exp);
         if (_kappa_set)
            theEquation->set_kappa(_kappa_exp);
         break;


      case LINEAR_ELASTICITY:
         switch (_dim) {

            case 2:
               if (spD=="feP1")
                  theEquation = new Elas2DT3(*_theMesh);
               else if (spD=="feQ1")
                  theEquation = new Elas2DQ4(*_theMesh);
               break;

            case 3:
               if (spD=="feP1")
                  theEquation = new Elas3DT4(*_theMesh);
               else if (spD=="feQ1")
                  theEquation = new Elas3DH8(*_theMesh);
               break;
         }
         if (_rho_set)
            theEquation->set_rho(_rho_exp);
         if (_young_set)
            theEquation->set_young(_young_exp);
         if (_poisson_set)
            theEquation->set_poisson(_poisson_exp);
         break;

      case TRUSS:
         switch (_dim) {

            case 2:
               if (spD=="feP1")
                  theEquation = new Bar2DL2(*_theMesh);
               break;
         }
         if (_rho_set)
            theEquation->set_rho(_rho_exp);
         if (_young_set)
            theEquation->set_young(_young_exp);
         break;

   case INCOMPRESSIBLE_NAVIER_STOKES:
         switch (_dim) {

            case 2:
               if (spD=="feP1")
                  theEquation = new TINS2DT3S(*_theMesh);
               break;

            case 3:
               if (spD=="feP1")
                  theEquation = new TINS3DT4S(*_theMesh);
               break;
         }
         if (_rho_set)
            theEquation->set_rho(_rho_exp);
         if (_mu_set)
            theEquation->set_mu(_mu_exp);
         if (_beta_set)
            theEquation->set_beta(_beta_exp);
         break;
   }
   if (ret==0 && theEquation->SolverIsSet()==false)
      ls = CG_SOLVER, prec = DILU_PREC;
}

} /* namespace RITA */
