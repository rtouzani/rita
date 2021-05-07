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

                         Implementation of class 'data'

  ==============================================================================*/

#include "data.h"
#include "configure.h"
#include "rita.h"
#include "linear_algebra/Matrix.h"
#include "io/IOField.h"

using std::cout;
using std::endl;
using std::map;
using std::pair;

namespace RITA {

data::data(rita*      r,
	   cmd*       command,
	   configure* config)
     : _rita(r), _size(0), _ifield(0), _imesh(0), _igrid(0), _ifct(0), _itab(0),
       _nb_fields(0), _default_field(0), _nb_fcts(0), _nb_meshes(0), _nb_tabs(0), _nb_grids(0),
       _verb(1), _configure(config), _cmd(command),
       _theMesh_alloc(0), _theTab_alloc(0), _theGrid_alloc(0), _theFct_alloc(0),
       _theVector_alloc(0), _theMatrix_alloc(0), _u_alloc(0),
       _theMesh(nullptr), _theTab(nullptr), _theGrid(nullptr), _theFct(nullptr), 
       _theVector(nullptr), _theMatrix(nullptr)
{
   _u = nullptr;
   _ret = 0; _key = 0;
   _nb_grids = _nb_tabs = _nb_meshes = _nb_fcts = iifct = _nb_vectors = _nb_matrices = 0;
   ok = false;
}


data::~data()
{
   if (_theGrid_alloc)
      delete _theGrid;
   if (_theMesh_alloc)
      delete _theMesh;
   if (_theFct_alloc)
      delete _theFct;
   if (_theTab_alloc)
      delete _theTab;
   if (_theVector_alloc)
      delete _theVector;
   if (_theMatrix_alloc)
      delete _theMatrix;
   if (_u_alloc)
      delete _u;
}


int data::checkField(string name)
{
   auto it = find(Field.begin(),Field.end(),name);
   if (it != Field.end())
      return distance(Field.begin(),it);
   else
      return -1;
}


int data::checkMesh(string name)
{
   auto it = find(mesh_name.begin(),mesh_name.end(),name);
   if (it != mesh_name.end())
      return distance(mesh_name.begin(),it);
   else
      return -1;
}


int data::checkGrid(string name)
{
   auto it = find(grid_name.begin(),grid_name.end(),name);
   if (it != grid_name.end())
      return distance(grid_name.begin(),it);
   else
      return -1;
}


int data::checkFct(string name)
{
   for (auto it=theFct.begin(); it!=theFct.end(); ++it) {
      if ((*it)->name==name)
         return distance(theFct.begin(),it);
   }
   return -1;
}


int data::addFunction(const string&         name,
		      const string&         def,
		      const vector<string>& var)
{
   string nm = name;
   if (nm!="") {
      int k = checkFct(name);
      if (k>=0) {
         _rita->msg("rita>data>","Function "+nm+" exists already.");
         return 1;
      }
   }
   if (nm=="")
      nm = "F-" + to_string(_nb_fcts+1);
   _theFct = new OFELI::Fct(nm,def,var);
   if (_theFct->check()) {
      _rita->msg("data>",_theFct->getErrorMessage());
      delete _theFct;
      return 1;
   }
   _theFct_alloc = 1;
   theFct.push_back(_theFct);
   _nb_fcts++;
   return 0;
}


int data::addMesh(OFELI::Mesh* ms,
                  string       name)
{
   theMesh.push_back(ms);
   mesh_name.push_back(name);
   return ++_nb_meshes;
}


int data::addField(string   name,
                   dataSize s,
                   int      n,
                   int      nb_dof)
{
   bool new_field = true;
   _ifield = _nb_fields;
   int k = checkField(name);
   if (k>=0) {
      _ifield = k;
      new_field = false;
   }

   if (s==GIVEN_SIZE)
      _size = n;
   _nb_dof = 1;
   if (new_field) {
      Field.push_back(name);
      FieldName[name] = _ifield;
      FieldSizeType.push_back(s);
      FieldEquation.push_back(0);
      FieldType.push_back(ALGEBRAIC_EQ);
   }
   bool wm = (s==NODES || s==ELEMENTS || s==SIDES || s==EDGES);
   if (s==GIVEN_SIZE) {
      _u = new OFELI::Vect<double>;
      _u_alloc = 1;
      _u->setSize(n);
      _u->setName(name);
      u.push_back(_u);
      _nb_dof = 1;
      if (new_field)
         _nb_fields++;
      return _ifield;
   }
   else if (s==GRID) {
      _u = new OFELI::Vect<double>;
      _u_alloc = 1;
      _u->setSize(n);
      _u->setName(name);
      u.push_back(_u);
      return _ifield;
   }
   else if (wm) {
      if (_rita->_theMesh==nullptr) {
         _rita->msg("data>","No mesh data available.");
         _ret = -1;
         return _ret;
      }
   }
   if (s==NODES) {
      _theMesh = _rita->_theMesh;
      if (_theMesh->getNbNodes()==0) {
         _rita->msg("data>","Mesh has no nodes");
         _ret = 1;
         return _ret;
      }
      _nb_dof = nb_dof;
   }
   else if (s==ELEMENTS) {
      _theMesh = _rita->_theMesh;
      if (_theMesh->getNbElements()==0) {
         _rita->msg("data>","Mesh has no elements.");
         _ret = 1;
         return _ret;
      }
      _nb_dof = nb_dof;
   }
   else if (s==SIDES) {
      _theMesh = _rita->_theMesh;
      if (_theMesh->getNbSides()==0) {
         _rita->msg("data>","Mesh has no sides");
         _ret = 1;
         return _ret;
      }
      _nb_dof = nb_dof;
   }
   else if (s==EDGES) {
      _theMesh = _rita->_theMesh;
      if (_theMesh->getNbEdges()==0) {
         _rita->msg("data>","Mesh has no edges");
         _ret = 1;
         return _ret;
      }
      _nb_dof = nb_dof;
   }
   _u = new OFELI::Vect<double>;
   _u_alloc = 1;
   _u->setName(name);
   u.push_back(_u);
   if (s==NODES)
      _u->setMesh(*_theMesh,NODE_DOF,nb_dof);
   else if (s==ELEMENTS)
      _u->setMesh(*_theMesh,ELEMENT_DOF,nb_dof);
   else if (s==SIDES)
      _u->setMesh(*_theMesh,SIDE_DOF,nb_dof);
   else if (s==EDGES)
      _u->setMesh(*_theMesh,EDGE_DOF,nb_dof);
   if (new_field)
      _nb_fields++;
   return _ifield;
}


int data::run()
{
   *_rita->ofh << "data" << endl;
   while (1) {
      _cmd->readline("rita>data> ");
      if (_nb < 0)
         continue;
      switch (_key=_cmd->getKW(_kw_data)) {

         case 0:
         case 1:
            _cmd->setNbArg(0);
            cout << "\nAvailable Commands:\n";
            cout << "grid:       Define a grid\n";
            cout << "mesh:       Define a finite element mesh\n";
            cout << "field:      Define a field\n";
            cout << "tabulation: Define a tabulated function\n";
            cout << "function:   Define a function\n";
            cout << "vector:     Define a vector\n";
            cout << "matrix:     Define a matrix\n";
            cout << "summary:    Summary of prescribed data\n";
            cout << "end or <:   go back to higher level" << endl;
            break;

         case 2:
            _rita->setConfigure();
            break;

         case 3:
            _ret = setGrid();
            break;

         case 4:
            _ret = 10;
            return _ret;

         case 5:
            _ret = setField();
            break;

         case 6:
            _ret = setTab();
            break;

         case 7:
            _ret = setFunction();
            break;

         case 8:
            _ret = setVector();
            break;

         case 9:
            _ret = setMatrix();
            break;

         case 10:
	   //            Clear();
            break;

         case 11:
            Summary();
            break;

         case 12:
         case 13:
            _ret = 0;
            ok = true;
            *_rita->ofh << "  end" << endl;
            return _ret;

         case 14:
         case 15:
            _ret = 100;
            return _ret;

         case 16:
            _ret = 200;
            return _ret;

         case -4:
            return 1;

         default:
            _rita->msg("data>","Unknown Command "+_cmd->token(),
                       "Available commands: grid, mesh, field, tabulation, function, vector, matrix, summary\n"
                       "Global commands:    help, ?, set, <, end, quit, exit");
            break;
       }
   }
   return 0;
}


void data::getHelp()
{
   cout << "In command data, the following data can be defined:\n";
   cout << "field: Define a field (unknown)\n";
   cout << "function: Define a function\n";
   cout << "mesh: Clear all defined fields\n";
   cout << "summary: Display this summary\n\n";
   cout << "Global commands: \n";
   cout << "help or ?: Display this help\n";
   cout << "set: Set configuration data\n";
   cout << "end or <: Back to higher level\n";
   cout << "exit: End the program\n" << endl;
}


int data::setConfigure()
{
   _configure->setVerbose(_verb);
   int ret = _configure->run();
   _verb = _configure->getVerbose();
   return ret;
}


int data::setGrid()
{
   double xmin=0., xmax=1., ymin=0., ymax=1., zmin=0., zmax=1.;
   int nb=0, d1=-1, d2=-1, d3=-1, dim=0, nx=10, ny=10, nz=10;
   string name="G-"+to_string(_nb_grids+1);
   vector<string> kw {"name","min","max","ne"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args<=0) {
      _rita->msg("data>grid>","Error in command.");
      return 1;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case 0:
            name = _cmd->string_token(0);
            break;

         case 1:
            if (nb==0 || nb>3) {
               _rita->msg("data>grid>","Illegal number of arguments");
               return 1;
            }
            dim = d1 = nb;
            xmin = _cmd->double_token(0);
            if (d1>1)
               ymin = _cmd->double_token(1);
            if (d1>2)
               zmin = _cmd->double_token(2);
            break;

         case 2:
            if (nb==0 || nb>3) {
               _rita->msg("data>grid>","Illegal number of arguments");
               return 1;
            }
            dim = d2 = nb;
            xmax = _cmd->double_token(0);
            if (d2>1)
               ymax = _cmd->double_token(1);
            if (d2>2)
               zmax = _cmd->double_token(2);
            break;

         case 3:
            if (nb==0 || nb>3) {
               _rita->msg("data>grid>","Illegal number of arguments");
               return 1;
            }
            dim = d3 = nb;
            nx = _cmd->int_token(0);
            if (d3>1)
               ny = _cmd->int_token(1);
            if (d3>2)
               nz = _cmd->int_token(2);
            break;

         default:
            _rita->msg("data>grid>","Unknown argument: "+kw[n]);
            return 1;
      }
   }
   if ((d1!=d2 || d1!=d3 || d2!=d3) && d1>=0 && d2>=0 && d3>=0) {
      _rita->msg("data>grid>","Dimensions do not match.");
      return 1;
   }
   if (xmin>=xmax || ymin>=ymax || zmin>=zmax) {
      _rita->msg("data>grid>","Domain definition is incorrect.");
      return 1;
   }
   if (nx<1 || ny<1 || nz<1) {
      _rita->msg("data>grid>","Number of subdivisions is incorrect.");
      return 1;
   }
   *_rita->ofh << "  grid  name=" << name;
   if (dim==1) {
      _theGrid = new OFELI::Grid(xmin,xmax,nx);
      *_rita->ofh << "  min=" << xmin << "  max=" << xmax << "  ne=" << nx;
   }
   else if (dim==2) {
      _theGrid = new OFELI::Grid(xmin,xmax,ymin,ymax,nx,ny);
      *_rita->ofh << "  min=" << xmin << "," << ymin << "  max=" << xmax << "," << ymax
                 << "  ne=" << nx << "," << ny;
   }
   else if (dim==3) {
      _theGrid = new OFELI::Grid(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz);
      *_rita->ofh << "  min=" << xmin << "," << ymin << "," << zmin << "  max=" << xmax
                 << "," << ymax << "," << zmax << "  ne=" << nx << "," << ny << "," << nz;
   }
   *_rita->ofh << endl;
   _theGrid_alloc = 1;
   theGrid.push_back(_theGrid);
   grid_name.push_back(name);
   _nb_grids++;
   return 0;
}


int data::setVector()
{
   int size=0, nb=0;
   string name="vect-"+to_string(_nb_vectors+1);
   vector<string> kw = {"name","size","def$ine","set"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args<=0) {
      _rita->msg("data>vector>","Error in command.","Available arguments: name, size, define, set.");
      return 1;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case 0:
            name = _cmd->string_token(0);
            break;

         case 1:
            size = _cmd->int_token(0);
            break;

         case 2:
            break;

         case 3:
            break;

         default:
            _rita->msg("data>vector>","Unknown argument: "+kw[n]);
            return 1;
      }
   }
   if (nb_args>0) {
      *_rita->ofh << endl;
   }
   theVector.push_back(_theVector);
   vector_name.push_back(name);
   _nb_vectors++;
   return 0;
}


int data::setMatrix()
{
   string name="mat-"+to_string(_nb_matrices+1);
   vector<string> kw = {"name","storage","size","def$ine","set"};
   theMatrix.push_back(_theMatrix);
   matrix_name.push_back(name);
   _nb_matrices++;
   return 0;
}


int data::setField()
{
   int nb=0;
   double umin=0., umax=0.;
   string name="u", type="size", arg=" ", mn="";
   vector<string> kw = {"name","size","mesh","grid","nbdof","type","uniform"};
   vector<string> types = {"size","nodes","elements","sides","edges"};
   map<string,dataSize> tt = {{"size",GIVEN_SIZE}, {"nodes",NODES}, {"elements",ELEMENTS},
                              {"sides",SIDES}, {"edges",EDGES}};
   if (_rita->_theMesh!=nullptr)
      _theMesh = _rita->_theMesh;
   if (_default_field==1) {
      _ifield = 0;
      FieldName.clear();
   }
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args<=0) {
      _rita->msg("data>field>","Error in command.",
                 "Available arguments: name, size, mesh, grid, nbdof, type, uniform.");
      return 1;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg();
      switch (n) {

         case 0:
            name = _cmd->string_token();
            break;

         case 1:
            _size = _cmd->int_token();
            type = "size";
            break;

         case 2:
            mn = _cmd->string_token();
	    //            mesh_ok = 1;
            break;

         case 3:
            _nb_dof = _cmd->int_token();
            break;

         case 4:
            _nb_dof = _cmd->int_token();
            break;

         case 5:
            type = _cmd->string_token();
            if (type!="size" && type!="nodes" && type!="elements" &&
                type!="sides" && type!="edges") {
               _rita->msg("data>field>","Unknown type: "+type);
               return 0;
            }
            break;

         case 6:
            umin = _cmd->double_token();
            if (nb>1)
               umax = _cmd->double_token();
	    //            uniform_ok = nb;
            break;

         default:
            _rita->msg("data>field>","Unknown argument: "+kw[n]);
            return 1;
      }
   }
   if (nb_args>0) {
      int k = addField(name,tt[type],_size);
      if (k==0) {
         if (_verb>1) {
            if (type=="size") 
               cout << "New field: " << name << ", size: " << _size << endl;
            else
               cout << "New field: " << name << ", Nb. of DOF: " << _nb_dof << endl;
            cout << "Total number of fields: " << _nb_fields << endl;
         }
         *_rita->ofh << "  field  name=" << name;
         if (_size)
            *_rita->ofh << " size=" << _size;
         *_rita->ofh << " type=" << type << "  nbdof=" << _nb_dof << " min=" << umin << " max=" << umax;
         nb_dof[_ifield] = _nb_dof;
      /*      if (uniform==1) {
         _rita->msg("data>field>","Minimal and maximal values must be given for field.");
	 return 1;
	 }*/
         *_rita->ofh << endl;
      }
   }
   return 0;
}


int data::setFunction()
{
   _ret = 0;
   bool var_ok=false, def_ok=false, name_ok=false;
   string vv="", def="", name="f";
   int nb=1, ret=0;
   vector<string> var;
   vector<string> kw {"name","var","field","nb","def","definition"};

   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args<0) {
      _rita->msg("data>function>:","Error in command.");
      return 1;
   }
   for (int i=0; i<nb_args; ++i) {

      switch(_cmd->getArg("=")) {

         case 0:
            name = _cmd->string_token();
            name_ok = true;
            break;

         case 1:
         case 2:
            vv = _cmd->string_token();
            var_ok = true;
            break;

         case 3:
            nb = _cmd->int_token();
            break;

         case 4:
         case 5:
            def = _cmd->string_token();
            def_ok = true;
            break;

         default:
            _rita->msg("data>function>","Unknown argument: "+_cmd->Arg());
            return 1;
      }
      if (!name_ok) {
         _rita->msg("data>function>","Missing function name.");
         return 1;
      }
   }
   if (nb_args==0) {
      _rita->msg("data>function>","No command argument given.");
      return 1;
   }
   else {
      if (!var_ok) {
         _rita->msg("data>function>","No variable defined.");
         return 1;
      }
      if (!def_ok) {
         _rita->msg("data>function>","No function definition given.");
         return 1;
      }
      for (const auto& v: theFct) {
         if (v->name==name) {
            _rita->msg("data>function>","Function "+name+" exists already.");
            ret = 1;
         }
      }
      if (ret)
         return 1;
      else {
         if (nb==1)
            var.push_back(vv);
         else {
            for (int i=0; i<nb; ++i)
               var.push_back(vv+to_string(i+1));
         }
         addFunction(name,def,var);
         if (name_ok)
            *_rita->ofh << "  function  name=" << name;
         for (const auto& v: var)
            *_rita->ofh << " var=" << v;
         *_rita->ofh << "  def=" << def << endl;
      }
   }
   return 0;
}


int data::setNbDOF()
{
   if (_cmd->setNbArg(1,"Give number of degrees of freedom.")) {
      _rita->msg("data>>nbdof>","Missing value of nbdof.");
      return 1;
   }
   _ret = _cmd->get(_nb_dof);
   if (!_ret)
      *_rita->ofh << "  nbdof " << _nb_dof;
   return _ret;
}

/*
void data::Clear()
{
   *_rita->ofh << "  clear" << endl;
   _nb_fields = 0;
   Field[0] = "u";
   FieldName["u"] = 0;
   FieldEquation[0] = 0;
   if (_verb)
      cout << "Field data cleared." << endl;
   if (u[_ifield].size()>0) {
      u[_ifield] = 0.;
   }
   if (bc[_ifield].size()>0) {
      bc[_ifield] = 0.;
      if (_verb)
         cout << "Boundary condition vector cleared." << endl;
   }
   if (sf[_ifield].size()>0) {
      sf[_ifield] = 0.;
      if (_verb)
         cout << "Surface force condition vector cleared." << endl;
   }
   if (bf[_ifield].size()>0) {
      bf[_ifield] = 0.;
      if (_verb)
         cout << "Body force vector cleared." << endl;
   }
   _ifield = 0;
   _ret = 90;
}*/


int data::setTab()
{
   int dim1=0, dim2=0, dim3=0;
   string file="";
   string name = "T" + to_string(theTab.size()+1);
   double xmin=0., xmax=1., ymin=0., ymax=1., zmin=0., zmax=1.;
   int nb=0, nx=10, ny=10, nz=10, grid_ok=0, file_ok=0, field_ok=0;
   vector<string> kw = {"name","file","min","max","ne","field"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args<0) {
      _rita->msg("data>tabulation>","Error in command.");
      return 1;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case 0:
            name = _cmd->string_token(0);
            break;

         case 1:
            file = _cmd->string_token(0);
            file_ok = 1;
            break;

         case 2:
            dim1 = nb;
            xmin = _cmd->double_token(0);
            if (nb>1)
               ymin = _cmd->double_token(1);
            if (nb>2)
               zmin = _cmd->double_token(2);
            grid_ok += 1;
            break;

         case 3:
            dim2 = nb;
            xmax = _cmd->double_token(0);
            if (nb>1)
               ymax = _cmd->double_token(1);
            if (nb>2)
               zmax = _cmd->double_token(2);
            grid_ok += 10;
            break;

         case 4:
            dim3 = nb;
            nx = _cmd->int_token(0);
            if (nb>1)
               ny = _cmd->int_token(1);
            if (nb>2)
               nz = _cmd->int_token(2);
            grid_ok += 100;
            break;

         case 5:
	   //            fd = _cmd->string_token();
            field_ok = true;
            break;

         default:
            _rita->msg("data>tabulation>","Unknown argument: "+kw[n]);
            return 1;
      }
   }
   if (!file_ok && grid_ok<111) {
      _rita->msg("data>tabulation>","No grid data given.");
      return 1;
   }
   if (!field_ok && !file_ok) {
      _rita->msg("data>tabulation>","No associated field given.");
      return 1;
   }
   if (!file_ok && (dim1!=dim2 || dim1!=dim3 || dim2!=dim3)) {
      _rita->msg("data>tabulation>","Incompatible space dimensions as given by grid data.");
      return 1;
   }
   *_rita->ofh << "  tabulation  name=" << name;
   if (file_ok) {
      ifstream ip(file);
      if (ip.is_open())
         ip.close();
      else {
         _rita->msg("data>tabulation>","Unable to open file: "+file);
         return 1;
      }
      _theTab = new OFELI::Tabulation(file);
      _theTab_alloc = 1;
      *_rita->ofh << "  file=" << file;
   }
   else {
      if (dim1==1)
         _theGrid = new OFELI::Grid(xmin,xmax,nx);
      else if (dim1==2)
         _theGrid = new OFELI::Grid(xmin,xmax,ymin,ymax,nx,ny);
      if (dim1==3)
         _theGrid = new OFELI::Grid(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz);
   }
   tab_name.push_back(name);
   theTab.push_back(_theTab);
   _nb_tabs++;
   *_rita->ofh << endl;
   return 0;
}


void data::Summary()
{
   int nb = 0;
   cout << "SUMMARY OF DATA:" << endl;
   cout << "===============================================================" << endl;
   cout << "Number of fields: " << _nb_fields << endl;
   for (const auto& v: u)
      cout << "Field: " << v->getName() << ", size: " << v->size() << endl;
   cout << "===============================================================" << endl;
   cout << "Number of functions: " << _nb_fcts << endl;
   for (const auto& v: theFct) {
      cout << "Function: " << v->name << ", Variable(s): ";
      for (size_t j=0; j<v->nb_var-1; ++j)
         cout << v->var[j] << ",";
      cout << v->var[v->nb_var-1] << ", Definition: " << v->expr << endl;
   }
   cout << "---------------------------------------------------------------" << endl;
   nb = 0;
   cout << "\nNumber of fields: " << _nb_fields << endl;
   for (int i=0; i<_nb_fields; ++i) {
      if (FieldSizeType[i]==GIVEN_SIZE)
         cout << "Field: " << Field[i] << ", Size: " << u[i]->size() << endl;
      else
         cout << "Field: " << Field[i] << ", Number of degrees of freedom: " << nb_dof[i] << endl;
   }
   cout << "---------------------------------------------------------------" << endl;
   cout << "\nNumber of tabulated functions: " << _nb_tabs << endl;
   cout << "---------------------------------------------------------------" << endl;
   cout << "\nNumber of grids: " << _nb_grids << endl;
   nb = 0;
   for (const auto& v: theGrid) {
      cout << "Grid No.            " << nb++ << endl;
      cout << "Grid name:          " << grid_name[nb-1] << endl;
      cout << "Space dimension:    " << v->getDim() << endl;
      if (v->getDim()==1)
         cout << "Domain:                    (" << v->getX(1) << ","
              << v->getX(v->getNx()+1) << ")" << endl;
      if (v->getDim()==2)
         cout << "Domain:                    (" << v->getX(1) << ","
              << v->getX(v->getNx()+1) << ")x(" << v->getY(1) << ","
              << v->getY(v->getNy()+1) << ")" << endl;
      if (v->getDim()==3)
         cout << "Domain:                    (" << v->getX(1) << ","
              << v->getX(v->getNx()+1) << ")x(" << v->getY(1) << ","
              << v->getY(v->getNy()+1) << ")x(" << v->getZ(1) << ","
              << v->getZ(v->getNz()+1) << ")" << endl;
      cout << "Number of x-intervals:    " << v->getNx() << endl;
      cout << "Number of y-intervals:    " << v->getNy() << endl;
      cout << "Number of z-intervals:    " << v->getNz() << endl;
   }
   cout << "---------------------------------------------------------------" << endl;
   cout << "\nNumber of meshes: " << _nb_meshes << endl;
   nb = 0;
   for (const auto& v: theMesh) {
      cout << "Mesh No.            " << nb++ << endl;
      cout << "Mesh name:          " << mesh_name[nb-1] << endl;
      cout << "Space dimension:    " << v->getDim() << endl;
      cout << "Number of nodes:    " << v->getNbNodes() << endl;
      cout << "Number of elements: " << v->getNbElements() << endl;
      cout << "Number of sides:    " << v->getNbSides() << endl;
   }
   cout << "===============================================================" << endl;
}

} /* namespace RITA */
