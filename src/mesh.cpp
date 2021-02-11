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
 
                      Implementation of class 'mesh'

  ==============================================================================*/

#include <iostream>
#include <stdlib.h>
#include "mesh.h"
#include "mesh/saveMesh.h"
#include "rita.h"
#include "data.h"

#ifdef USE_GMSH
#include <gmsh.h>
#endif

namespace RITA {

mesh::MeshData_Ptr mesh::MESH_DATA [] = { nullptr,
                                          &mesh::getHelp,
                                          &mesh::getHelp,
                                          &mesh::setConfigure,
                                          &mesh::set1D,
                                          &mesh::setRectangle,
                                          &mesh::setCube,
                                          &mesh::setPoint,
                                          &mesh::setCurve,
                                          &mesh::setSurface,
                                          &mesh::setVolume,
                                          &mesh::setContour,
                                          &mesh::setCode,
                                          &mesh::Generate,
                                          &mesh::setNbDOF,
                                          &mesh::List,
                                          &mesh::Plot,
                                          &mesh::Clear,
                                          &mesh::Save,
                                          &mesh::Read
                                        };

mesh::mesh(rita*      r,
           cmd*       command,
           configure* config)
     : _rita(r), _saved(false), _generated(false), _geo(false), _ret(0), _generator(0),
       _theMesh(nullptr), _theDomain(nullptr), _nb_Ccontour(0), _nb_Scontour(0), _nb_Vcontour(0),
       _nb_sub_domain(0), _nb_point(0), _nb_curve(0), _nb_surface(0), _nb_volume(0), _configure(config),
       _cmd(command)
{
}


mesh::~mesh()
{
   if (_theMesh!=nullptr)
      delete _theMesh, _theMesh = nullptr;
   if (_theDomain!=nullptr)
      delete _theDomain, _theDomain = nullptr;
}


int mesh::run()
{
   _nb_dof = 1;
   _data = _rita->_data;
   const vector<string> kw = {"help","?","set","1d","rect$angle","cube","point","curve",
                              "surface","volume","contour","code","gen$erate","nbdof",
                              "list","plot","clear","save","read","end","<","quit","exit","EXIT"};
#ifndef USE_GMSH
   _theDomain = new OFELI::Domain;
#endif
   while (1) {
      int nb_args = _cmd->readline("rita>mesh> ");
      if (nb_args < 0)
         continue;
      _key = _cmd->getKW(kw);
      if (_key<-1)
         continue;
      if (_key==-1) {
         cout << "Unknown Command: " << _cmd->token() << endl;
         cout << "Available commands:\n";
         cout << "1d, rectangle, cube, point, curve, surface, volume, contour, code\n";
         cout << "generate, list, plot, clear, save, read, end, <" << endl;
         cout << "Global commands:\n";
         cout << "help, ?, set, quit, exit" << endl;
         *_ofl << "In rita>mesh>: Unknown Command " << _cmd->token() << endl;
         continue;
      }
      else if (_key==19 || _key==20) {
          *_ofh << "  end" << endl;
          _ret = 0;
         return _ret;
      }
      else if (_key>20) {
         _ret = 100;
         return _ret;
      }
      (this->*MESH_DATA[_key+1])();
      if (_ret>=90)
         return _ret;
   }
   return 1;
}


void mesh::getHelp()
{
   _cmd->setNbArg(0);
   cout << "\nAvailable commands:" << endl;
   cout << "1d        : Data to generate a 1-D mesh" << endl;
   cout << "rectangle : Data to mesh a rectangle" << endl;
   cout << "cube      : Data to mesh a cube (parallelepiped)" << endl;
   cout << "point     : Define a point" << endl;
   cout << "curve     : Define a curve" << endl;
   cout << "surface   : Define a surface" << endl;
   cout << "volume    : Define a volume" << endl;
   cout << "contour   : Define a contour as a sequence of curves or surfaces" << endl;
//   cout << "subdomain : Define a subdomain" << endl;
   cout << "code      : Set code for points, lines, surfaces, volumes" << endl;
   cout << "generate  : Generate mesh of a polygon" << endl;
   cout << "list      : List mesh data" << endl;
   cout << "plot      : Plot mesh" << endl;
   cout << "clear     : Clear mesh" << endl;
   cout << "read      : Read mesh from file" << endl;
   cout << "save      : Save mesh in file" << endl;
   cout << "end or <  : Return to higher level" << endl;
}


void mesh::List()
{
   if (_generated==0) {
      cout << "Geometry data\n" << endl;
      cout << "Space dimension:    " << _theMesh->getDim() << endl;
      cout << "Number of nodes:    " << _theMesh->getNbNodes() << endl;
      cout << "Number of elements: " << _theMesh->getNbElements() << endl;
      cout << "Number of sides:    " << _theMesh->getNbSides() << endl;
   }
   else {
      cout << "Mesh data\n" << endl;
      cout << "Space dimension:    " << _theMesh->getDim() << endl;
      cout << "Number of nodes:    " << _theMesh->getNbNodes() << endl;
      cout << "Number of elements: " << _theMesh->getNbElements() << endl;
      cout << "Number of sides:    " << _theMesh->getNbSides() << endl;
   }
}


void mesh::setConfigure()
{
   _configure->setVerbose(_verb);
   _configure->run();
   _verb = _configure->getVerbose();
}


void mesh::set1D()
{
   _dim = 1;
   _nb_dof = 1;
   double xmin=0., xmax=1.;
   int nb=0, ret=0, ne=10, cmin=0, cmax=0;
   _mesh_file = "rita-1d.m";
   _saved = false;
   _ret = 0;
   const vector<string> kw = {"help","?","set","domain","ne","codes","nbdof","save","end",
                              "<","quit","exit","EXIT"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case 3:
            if (nb==1)
               xmin = xmax = _cmd->double_token(0);
            else if (nb==2) {
               xmin = _cmd->double_token(0);
               xmax = _cmd->double_token(1);
            }
            else {
               cout << "Error: This argument requires 1 or 2 parameters." << endl;
               *_ofl << "In rita>mesh>1d>: This argument requires 1 or 2 parameters." << endl;
               _ret = 1;
               return;
            }
            break;

         case 4:
            ne = _cmd->int_token(0);
            break;

         case 5:
            if (nb==1)
               cmin = cmax = _cmd->int_token(0);
            else if (nb==2) {
               cmin = _cmd->int_token(0);
               cmax = _cmd->int_token(1);
            }
            else {
               cout << "Error: This argument requires 1 or 2 parameters." << endl;
               *_ofl << "In rita>mesh>1d>: This argument requires 1 or 2 parameters." << endl;
               _ret = 1;
               return;
            }
            break;

         case 6:
            _nb_dof = _cmd->int_token(0);
            break;

         case 7:
            _mesh_file = _cmd->string_token(0);
            break;

         default:
	    cout << "Error: Unknown argument: " << kw[n] << endl;
            *_ofl << "In rita>mesh>1d>: Unknown argument: " << kw[n] << endl;
	    return;
      }
   }
   if (nb_args>0) {
      if (xmax<=xmin) {
         cout << "Error: xmax must be > xmin: " << xmin << endl;
         *_ofl << "In rita>mesh>1d>: Error in values of xmin and xmax: " << xmin << " " << xmax << endl;
         _ret = 1;
         return;
      }
      if (ne<2) {
         cout << "Error: Number of elements must be > 1" << endl;
         *_ofl << "In rita>mesh>1d>: Illegal number of elements." << endl;
         _ret = 1;
         return;
      }
      if (_theMesh!=nullptr)
         delete _theMesh, _theMesh = nullptr;
      _theMesh = new OFELI::Mesh(xmin,xmax,ne,cmin,cmax,1,size_t(_nb_dof));
      _data->addMesh(_theMesh,"M-"+to_string(_data->getNbMeshes()));
      _theMesh->removeImposedDOF();
      _saved = true;
      _generator = 1;
      _generated = true;
      _theMesh->put(_mesh_file);
      _data->mesh_name.push_back("M"+to_string(_data->theMesh.size()));
      _data->theMesh.push_back(_theMesh);
      if (_verb)
         cout << "1-D mesh complete and saved in file " << _mesh_file << endl;
      *_ofh << "  1d  domain=" << xmin << "," << xmax << "  codes=" << cmin
            << "," << cmax << "  ne=" << ne << "  nbdof=" << _nb_dof
            << "  save=" << _mesh_file << endl;
   }

   else {
      if (_verb) {
         cout << "Default interval: (0,1)\n";
         cout << "Default number of elements: 10\n";
         cout << "Default codes for end nodes: 0 0\n";
         cout << "Default nb of dof: 1" << endl;
      }
      *_ofh << "  1d" << endl;
      while (1) {
         if (_cmd->readline("rita>mesh>1d> ")<0)
            continue;
         switch (_key=_cmd->getKW(kw)) {

            case 0:
            case 1:
               _cmd->setNbArg(0);
               cout << "\nAvailable Commands:\n";
               cout << "domain   : Enter xmin and xmax\n";
               cout << "ne       : Enter number of elements to generate\n";
               cout << "codes    : Codes to associate to end nodes\n";
               cout << "nbdof    : Number of dof per node\n";
               cout << "save     : Save generated 1-D mesh to file and return to higher level\n";
               cout << "end or < : go back to higher level" << endl;
               break;

            case 2:
               if (_verb)
                  cout << "Setting configuration parameters ..." << endl;
               setConfigure();
               return;

            case 3:
               if (_verb>1)
                  cout << "Setting interval bounds ..." << endl;
               if (_cmd->setNbArg(2,"Give xmin and xmax.")) {
                  *_ofl << "In rita>mesh>1d>interval>: Missing values of xmin and xmax." << endl;
                  break;
               }
               ret  = _cmd->get(xmin);
               ret += _cmd->get(xmax);
               if (xmax<=xmin) {
                  cout << "Error: xmax must be > xmin: " << xmin << endl;
                  *_ofl << "In rita>mesh>1d>interval>: Error in values of xmin and xmax: " << xmin
                        << " " << xmax << endl;
                  break;
               }
               if (!ret)
                  *_ofh << "    interval " << xmin << " " << xmax << endl;
               break;

            case 4:
               if (_verb>1)
                  cout << "Setting interval subdivision ..." << endl;
               if (_cmd->setNbArg(1,"Give number of elements.")) {
                  *_ofl << "In rita>mesh>1d>ne>: Missing number of elements." << endl;
                  break;
               }
               _cmd->get(ne);
               *_ofh << "    ne " << ne << endl;
               break;

            case 5:
               if (_verb>1)
                  cout << "Setting boundary codes ..." << endl;
               if (_cmd->setNbArg(2,"Give cmin and cmax.")) {
                  *_ofl << "In rita>mesh>1d>codes>: Missing codes for end nodes." << endl;
                  break;
               }
               ret  = _cmd->get(cmin);
               ret += _cmd->get(cmax);
               if (!ret)
                  *_ofh << "    codes " << cmin << " " << cmax << endl;
               break;

            case 6:
               if (_verb>1)
                  cout << "Setting number of degrees of freedom ..." << endl;
               if (_cmd->setNbArg(1,"Give number of dof per node.")) {
                  *_ofl << "In rita>mesh>1d>nbdof>: Missing number of dof." << endl;
                  break;
               }
               ret = _cmd->get(_nb_dof);
               if (!ret)
                  *_ofh << "    nbdof " << _nb_dof << endl;
               break;

            case 7:
               if (_verb>1)
                  cout << "Saving mesh in file ..." << endl;
               if (_saved) {
                  cout << "Warning: You are trying to delete an existing mesh.\n";
                  cout << "retype command 'save' to confirm." << endl;
                  *_ofl << "In rita>mesh>1d>save>: Trying to delete a created mesh." << endl;
                  _saved = false;
                  break;
               }
               if (_theMesh!=nullptr)
                  delete _theMesh, _theMesh = nullptr;
               _theMesh = new OFELI::Mesh(xmin,xmax,ne,cmin,cmax,1,size_t(_nb_dof));
               _theMesh->removeImposedDOF();
               _saved = true;
               _generator = 1;
               _generated = true;
               _mesh_file = "rita-1d.m";
               if (_cmd->getNbArgs()>0)
                  ret = _cmd->get(_mesh_file);
               if (!ret) {
                  _theMesh->put(_mesh_file);
                  _data->mesh_name.push_back("M"+to_string(_data->theMesh.size()));
                  _data->theMesh.push_back(_theMesh);
                  *_ofh << "  save " << _mesh_file << endl;
                  cout << "1-D mesh complete and saved in file " << _mesh_file << endl;
               }
               break;

            case 8:
            case 9:
               if (_verb>1)
                  cout << "Getting back to higher level ..." << endl;
               if (!_saved) {
                  if (_theMesh!=nullptr)
                     delete _theMesh, _theMesh = nullptr;
                  _theMesh = new OFELI::Mesh(xmin,xmax,ne,cmin,cmax,1,size_t(_nb_dof));
                  _theMesh->removeImposedDOF();
                  _saved = true;
               }
               *_ofh << "  end" << endl;
               if (_verb && !_saved)
                  cout << "Mesh '1d' complete." << endl;
               _ret = 90;
               return;

            case 10:
            case 11:
               _ret = 100;
               return;

            case 12:
               _ret = 200;
               return;

            case -2:
               break;

            default:
               cout << "Unknown Command: " << _cmd->token() << endl;
               cout << "Available commands:\ninterval, ne, codes, nbdof, save, end, <" << endl;
               cout << "Global commands:\nhelp, ?, set, quit, exit" << endl;
               *_ofl << "In rita>mesh>1d>: Unknown command: " << _cmd->token() << endl;
               break;
         }
      }
   }
}


void mesh::setRectangle()
{
   double xmin=0., xmax=1., ymin=0., ymax=1.;
   int c[4] = {0,0,0,0}, cv[4] = {0,0,0,0};
   int nb=0, nx=10, ny=10, ret=0, ret1=0, ret2=0, ret3=0, ret4=0;
   _nb_dof = 1;
   _dim = 2;
   _ret = 0;
   _saved = false;
   _mesh_file = "rita-rectangle.m";
   const vector<string> kw = {"help","?","set","min","max","ne","codes","nbdof","save","end","<",
                              "quit","exit","EXIT"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case  3:
            if (nb==1)
               xmin = ymin = _cmd->double_token(0);
            if (nb==2) {
               xmin = _cmd->double_token(0);
               ymin = _cmd->double_token(1);
            }
	    else {
               cout << "Error: This argument requires 1 or 2 parameters." << endl;
               *_ofl << "In rita>mesh>rectangle>: This argument requires 1 or 2 parameters." << endl;
               _ret = 1;
               return;
	    }
            break;

         case  4:
            if (nb==1)
               xmax = ymax = _cmd->double_token(0);
            if (nb==2) {
               xmax = _cmd->double_token(0);
               ymax = _cmd->double_token(1);
            }
	    else {
               cout << "Error: This argument requires 1 or 2 parameters." << endl;
               *_ofl << "In rita>mesh>rectangle>: This argument requires 1 or 2 parameters." << endl;
               _ret = 1;
               return;
	    }
            break;

         case  5:
            if (nb==1)
               nx = ny = _cmd->int_token(0);
            else {
               nx = _cmd->int_token(0);
               ny = _cmd->int_token(1);
	    }
            break;

         case  6:
            if (nb==1)
               c[0] = c[1] = c[2] = c[3] = _cmd->int_token(0);
            else if (nb==4) {
               c[0] = _cmd->int_token(0);
               c[1] = _cmd->int_token(1);
               c[2] = _cmd->int_token(2);
               c[3] = _cmd->int_token(3);
            }
            else {
               cout << "Error: This argument requires 1 or 4 parameters." << endl;
               *_ofl << "In rita>mesh>rectangle>: This argument requires 1 or 4 parameters." << endl;
               _ret = 1;
               return;
            }
            break;

         case  7:
            _nb_dof = _cmd->int_token(0);
            break;

         case 8:
            _mesh_file = _cmd->string_token(0);
            break;

         default:
	    cout << "Error: Unknown argument: " << kw[n] << endl;
            *_ofl << "In rita>mesh>rectangle>: Unknown argument: " << kw[n] << endl;
	    return;
      }
   }
   if (nb_args>0) {
      if (xmax<=xmin) {
         cout << "Error: xmax must be > xmin: " << xmin << endl;
         *_ofl << "In rita>mesh>rectangle>: Error in values of xmin and xmax: " << xmin << " " << xmax << endl;
         return;
      }
      if (ymax<=ymin) {
         cout << "Error: ymax must be > ymin: " << xmin << endl;
         *_ofl << "In rita>mesh>rectangle>: Error in values of ymin and ymax: " << ymin << " " << ymax << endl;
         return;
      }
      if (!_saved) {
         if (_theMesh!=nullptr)
            delete _theMesh, _theMesh = nullptr;
         _theMesh = new OFELI::Mesh(xmin,xmax,ymin,ymax,nx,ny,c[3],c[1],c[0],c[2],TRIANGLE,size_t(_nb_dof));
         _data->mesh_name.push_back("M"+to_string(_data->theMesh.size()));
         _data->theMesh.push_back(_theMesh);
         _saved = _generated = true;
      }
      _generator = 1;
      _generated = true;
      *_ofh << "  rectangle min=" << xmin << "," << ymin << " max=" << xmax << "," << ymax;
      *_ofh << "  ne=" << nx << "," << ny << "  codes=" << c[0] << "," << c[1] << "," << c[2];
      *_ofh << "," << c[3] << "  nbdof=" << _nb_dof << "  save=" << _mesh_file << endl;
      _ret = 0;
      return;
   }

   else {
      if (_verb) {
         cout << "Default rectangle: (0,1)*(0,1)\n";
         cout << "Default node codes: 0\n";
         cout << "Default subdivision: 10*10\n";
         cout << "Default nb of dof: 1" << endl;
      }
      *_ofh << "  rectangle" << endl;

      while (1) {
         if (_cmd->readline("rita>mesh>rectangle> ")<0)
            continue;
         switch (_key=_cmd->getKW(kw)) {

            case 0:
            case 1:
               _cmd->setNbArg(0);
               cout << "\nAvailable commands:\n";
               cout << "min      : Values of xmin and ymin\n";
               cout << "max      : Values of xmax and ymax\n";
               cout << "ne       : Number of elements along x- and y-axes\n";
               cout << "codes    : Code to assign to nodes (or sides if < 0) on boundary lines\n";
               cout << "nbdof    : Number of dof per node\n";
               cout << "save     : Save generated mesh to file and return to higher level\n";
               cout << "end or < : go back to higher level" << endl;
               break;

            case 2:
               if (_verb)
                  cout << "Setting configuration parameter(s) ..." << endl;
               setConfigure();
               return;

            case 3:
               if (_verb>1)
                  cout << "Setting xmin and ymin ..." << endl;
               if (_cmd->setNbArg(2,"Give values of xmin and ymin.")) {
                  *_ofl << "In rita>mesh>rectangle>min>: Missing xmin and ymin values" << endl;
                  break;
               }
               ret  = _cmd->get(xmin);
               ret += _cmd->get(ymin);
               if (!ret)
                  *_ofh << "    min " << xmin << " " << ymin << endl;
               break;

            case 4:
               if (_verb>1)
                  cout << "Setting xmax and ymax ..." << endl;
               if (_cmd->setNbArg(2,"Give values of ymax and ymax.")) {
                  *_ofl << "In rita>mesh>rectangle>max>: Missing xmax and ymax values" << endl;
                  break;
               }
               ret  = _cmd->get(xmax);
               ret += _cmd->get(ymax);
               if (!ret)
                  *_ofh << "    max " << xmax << " " << ymax << endl;
               break;

            case 5:
               if (_verb>1)
                  cout << "Setting mesh subdivisions ..." << endl;
               if (_cmd->setNbArg(2,"Give subdivision in the x and y-directions.")) {
                  *_ofl << "In rita>mesh>rectangle>ne>: Missing subvdivisions in x and y directions" << endl;
                  break;
               }
               ret  = _cmd->get(nx);
               ret += _cmd->get(ny);
               if (!ret)
                  *_ofh << "    ne " << nx << " " << ny << endl;
               break;

            case 6:
               if (_verb>1)
                  cout << "Setting boundary codes ..." << endl;
               if (_cmd->setNbArg(4,"Code to assign on the line y=ymin.")) {
                  *_ofl << "In rita>mesh>rectangle>codes>: Missing code to assign on the line y=ymin." << endl;
                  break;
               }
               ret1 = _cmd->get(c[0]);
               ret2 = _cmd->get(c[1]);
               ret3 = _cmd->get(c[2]);
               ret4 = _cmd->get(c[3]);
               if (!ret1 && !ret2 && !ret3 && !ret4) {
                  *_ofh << "    codes " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << endl;
                  if (cv[0]==0)
                     cv[0] = c[0];
                  if (cv[1]==0)
                     cv[1] = c[0];
                  if (cv[1]==0)
                     cv[1] = c[1];
                  if (cv[2]==0)
                     cv[2] = c[1];
                  if (cv[2]==0)
                     cv[2] = c[2];
                  if (cv[3]==0)
                     cv[3] = c[2];
                  if (cv[3]==0)
                     cv[3] = c[3];
                  if (cv[0]==0)
                     cv[0] = c[3];
               }
               break;

            case 7:
               if (_verb>1)
                  cout << "Setting number of degrees of freedom ..." << endl;
               if (_cmd->setNbArg(1,"Give number of dof per node.")) {
                  *_ofl << "In rita>mesh>rectangle>nbdof>: Missing number of dof." << endl;
                  break;
               }
               ret = _cmd->get(_nb_dof);
               if (!ret)
                  *_ofh << "    nbdof " << _nb_dof << endl;
               break;

            case 8:
               if (_verb>1)
                  cout << "Saving mesh in file ..." << endl;
               if (_saved) {
                  cout << "Warning: You are trying to delete an existing mesh.\n";
                  cout << "retype command 'save' to confirm." << endl;
                  *_ofl << "In rita>mesh>rectangle>save>: Trying to delete a created mesh." << endl;
                  _saved = false;
                  break;
               }
               if (_theMesh!=nullptr)
                  delete _theMesh, _theMesh = nullptr;
               _theMesh = new OFELI::Mesh(xmin,xmax,ymin,ymax,nx,ny,c[3],c[1],c[0],c[2],TRIANGLE,size_t(_nb_dof));
               _saved = true;
               _generator = 2;
               if (_cmd->getNbArgs()>0)
                  ret = _cmd->get(_mesh_file);
               if (!ret) {
                  _theMesh->put(_mesh_file);
                  _data->mesh_name.push_back("M"+to_string(_data->theMesh.size()));
                  _data->theMesh.push_back(_theMesh);
                  *_ofh << "    save " << _mesh_file << endl;
                  _saved = true;
                  cout << "2-D mesh of rectangle complete and saved in file " << _mesh_file << endl;
               }
               break;

            case 9:
            case 10:
               if (_verb>1)
                  cout << "Getting back to higher level ..." << endl;
               *_ofh << "    end" << endl;
               if (!_saved) {
                  if (_theMesh!=nullptr)
                     delete _theMesh, _theMesh = nullptr;
                  _theMesh = new OFELI::Mesh(xmin,xmax,ymin,ymax,nx,ny,c[3],c[1],c[0],c[2],TRIANGLE,size_t(_nb_dof));
                  _data->mesh_name.push_back("M"+to_string(_data->theMesh.size()));
                  _data->theMesh.push_back(_theMesh);
                  _saved = _generated = true;
               }
               if (_verb && !_saved)
                  cout << "Mesh 'rectangle' complete." << endl;
               _ret = 90;
               return;

            case 11:
            case 12:
               _ret = 100;
               return;

            case 13:
               _ret = 200;
               return;

            case -2:
               break;

            default:
               cout << "Unknown Command!" << endl;
               cout << "Available commands: min, max, ne, codes, nbdof, save, end, <" << endl;
               cout << "Global commands: help, ?, set, quit, exit" << endl;
               *_ofl << "In rita>mesh>rectangle>: Unknown command: " << _cmd->token() << endl;
               break;
         }
      }
   }
   _ret = 0;
   return;
}


void mesh::setCube()
{
   double xmin=0., xmax=1., ymin=0., ymax=1., zmin=0., zmax=1.;
   int nb=0, ret=0, nx=10, ny=10, nz=10, cxmin=0, cxmax=0, cymin=0, cymax=0, czmin=0, czmax=0;
   _nb_dof = 1;
   _dim = 3;
   *_ofh << "  cube" << endl;
   _saved = false;
   _ret = 0;
   if (_verb) {
      cout << "Default cube: (0,1)*(0,1)*(0,1)\n";
      cout << "Default node codes: 0\n";
      cout << "Default subdivision: 10*10*10\n";
      cout << "Default nb of dof: 1" << endl;
   }
   const vector<string> kw = {"help","?","set","min","max","ne","codes","nbdof","save","end","<",
                              "quit","exit","EXIT"};

   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case  3:
            if (nb==1)
               xmin = ymin = zmin = _cmd->double_token(0);
            if (nb==3) {
               xmin = _cmd->double_token(0);
               ymin = _cmd->double_token(1);
               zmin = _cmd->double_token(2);
            }
	    else {
               cout << "Error: This argument requires 1 or 3 parameters." << endl;
               *_ofl << "In rita>mesh>cube>: This argument requires 1 or 3 parameters." << endl;
               _ret = 1;
               return;
	    }
            break;

         case  4:
            if (nb==1)
               xmax = ymax = zmax = _cmd->double_token(0);
            if (nb==3) {
               xmax = _cmd->double_token(0);
               ymax = _cmd->double_token(1);
               zmax = _cmd->double_token(2);
            }
	    else {
               cout << "Error: This argument requires 1 or 3 parameters." << endl;
               *_ofl << "In rita>mesh>cube>: This argument requires 1 or 3 parameters." << endl;
               _ret = 1;
               return;
	    }
            break;

         case  5:
            if (nb==1)
               nx = ny = nz = _cmd->int_token(0);
            else {
               nx = _cmd->int_token(0);
               ny = _cmd->int_token(1);
               nz = _cmd->int_token(2);
	    }
            break;

         case  8:
            if (nb==1)
               cxmin = cxmax = cymin = cymax = czmin = czmax = _cmd->int_token(0);
            else if (nb==6) {
               cxmin = _cmd->int_token(0);
               cxmax = _cmd->int_token(1);
               cymin = _cmd->int_token(2);
               cymax = _cmd->int_token(3);
               czmin = _cmd->int_token(4);
               czmax = _cmd->int_token(5);
            }
            else {
               cout << "Error: This argument requires 1 or 6 parameters." << endl;
               *_ofl << "In rita>mesh>cube>: This argument requires 1 or 6 parameters." << endl;
               _ret = 1;
               return;
            }
            break;

         case  9:
            _nb_dof = _cmd->int_token(0);
            break;

         case 10:
            _mesh_file = _cmd->string_token(0);
            break;

         default:
	    cout << "Error: Unknown argument: " << kw[n] << endl;
            *_ofl << "In rita>mesh>cube>: Unknown argument: " << kw[n] << endl;
	    return;
      }
   }
   if (nb_args>0) {
      if (xmax<=xmin) {
	cout << "Error: xmax must be > xmin: " << xmin << ", " << xmax << endl;
         *_ofl << "In rita>mesh>cube>: Error in values of xmin and xmax." << endl;
         return;
      }
      if (ymax<=ymin) {
	cout << "Error: ymax must be > ymin: " << ymin << ", " << ymax << endl;
         *_ofl << "In rita>mesh>cube>: Error in values of ymin and ymax." << endl;
         return;
      }
      if (zmax<=zmin) {
	cout << "Error: zmax must be > zmin: " << zmin << ", " << zmax << endl;
         *_ofl << "In rita>mesh>cube>: Error in values of zmin and zmax." << endl;
         return;
      }
      if (!_saved) {
         if (_theMesh!=nullptr)
            delete _theMesh, _theMesh = nullptr;
         _theMesh = new OFELI::Mesh(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,cxmin,cxmax,cymin,
                                    cymax,czmin,czmax,HEXAHEDRON,size_t(_nb_dof));
         _data->mesh_name.push_back("M"+to_string(_data->theMesh.size()));
         _data->theMesh.push_back(_theMesh);
         _saved = _generated = true;
      }
      _generator = 1;
      _generated = true;
      *_ofh << "  cube min=" << xmin << "," << ymin << "," << zmin << " max=" << xmax << "," << ymax << "," << zmax;
      *_ofh << " ne=" << nx << "," << ny << "," << nz << " codes=" << cxmin << "," << cxmax << "," << cymin << ",";
      *_ofh << cymax << "," << czmin << "," << czmax << " nbdof=" << _nb_dof << " save=" << _mesh_file << endl;
      _ret = 0;
      return;
   }

   else {
      while (1) {
         if (_cmd->readline("rita>mesh>cube> ")<0)
            continue;
         switch (_key=_cmd->getKW(kw)) {

            case 0:
            case 1:
               _cmd->setNbArg(0);
               cout << "\nAvailable commands:\n";
               cout << "min      : Values of xmin, ymin and zmin\n";
               cout << "max      : Values of xmax, ymax and zmax\n";
               cout << "ne       : Number of elements along x- y- and z-axes\n";
               cout << "codes    : Codes to assign to nodes (or sides if < 0) on boundary faces\n";
               cout << "nbdof    : Number of dof per node\n";
               cout << "save     : Save generated mesh\n";
               cout << "end or < : go back to higher level" << endl;
               break;

            case 2:
               if (_verb)
                  cout << "Setting configuration parameter(s) ..." << endl;
               setConfigure();
               break;

            case 3:
               if (_verb>1)
                  cout << "Setting xmin, ymin and zmin ..." << endl;
               if (_cmd->setNbArg(3,"Give values of xmin, ymin and zmin.")) {
                  *_ofl << "In rita>mesh>cube>min>: Missing values of xmin, ymin and zmin." << endl;
                  break;
               }
               ret  = _cmd->get(xmin);
               ret += _cmd->get(ymin);
               ret += _cmd->get(zmin);
               if (!ret)
                  *_ofh << "    min " << xmin << " " << ymin << " " << zmin << endl;
               break;

            case 4:
               if (_verb>1)
               cout << "Setting xmax, ymax and zmax ..." << endl;
               if (_cmd->setNbArg(3,"Give values of xmax, ymax and zmax.")) {
                  *_ofl << "In rita>mesh>cube>max>: Missing values of xmax, ymax and zmax." << endl;
                  break;
               }
               ret  = _cmd->get(xmax);
               ret += _cmd->get(ymax);
               ret += _cmd->get(zmax);
               if (!ret)
                  *_ofh << "    max " << xmax << " " << ymax << " " << zmax << endl;
               break;

            case 6:
               if (_verb>1)
                  cout << "Setting mesh subdivisions ..." << endl;
               if (_cmd->setNbArg(3,"Give subdivision in the x-, y- and z-directions.")) {
                  *_ofl << "In rita>mesh>cube>ne>: Missing subvdivisions in x, y and z directions" << endl;
                  break;
               }
               ret  = _cmd->get(nx);
               ret += _cmd->get(ny);
               ret += _cmd->get(nz);
               if (!ret)
                  *_ofh << "    ne " << nx << " " << ny << " " << nz << endl;
               break;
   
            case 7:
               if (_verb>1)
                  cout << "Setting boundary codes ..." << endl;
               if (_cmd->setNbArg(6,"Codes to assign to faces.")) {
                  *_ofl << "In rita>mesh>cube>codes>: Missing codes to assign to faces." << endl;
                  break;
               }
               ret  = _cmd->get(cxmin);
               ret += _cmd->get(cxmax);
               ret += _cmd->get(cymin);
               ret += _cmd->get(cymax);
               ret += _cmd->get(czmin);
               ret += _cmd->get(czmax);
               if (!ret)
                  *_ofh << "    codes " << cxmin << " " << cxmax << " " << cymin << " " << cymax
                        << " " << czmin << " " << czmax << endl;
               break;

            case 8:
               if (_verb>1)
                  cout << "Setting number of degrees of freedom ..." << endl;
               if (_cmd->setNbArg(1,"Give number of dof per node.")) {
                  *_ofl << "In rita>mesh>cube>nbdof>: Missing number of dof." << endl;
                  break;
               }
               ret = _cmd->get(_nb_dof);
               if (!ret)
                  *_ofh << "    nbdof " << _nb_dof << endl;
               break;

            case 9:
               if (_verb>1)
                  cout << "Saving mesh in file ..." << endl;
               if (_saved) {
                  cout << "Warning: You are trying to delete an existing mesh.\n";
                  cout << "retype command 'save' to confirm." << endl; 
                  *_ofl << "In rita>mesh>cube>save>: Trying to delete a created mesh." << endl;
                  _saved = false;
                  break;
               }
               if (_theMesh!=nullptr)
                  delete _theMesh, _theMesh = nullptr;
               _mesh_file = "rita-cube.m";
               if (_cmd->getNbArgs()>0)
                  _cmd->get(_mesh_file);
               _theMesh = new OFELI::Mesh(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,cxmin,cxmax,cymin,
                                          cymax,czmin,czmax,HEXAHEDRON,size_t(_nb_dof));
               _theMesh->put(_mesh_file);
               _data->mesh_name.push_back("M"+to_string(_data->theMesh.size()));
               _data->theMesh.push_back(_theMesh);
               _generator = 3;
               _generated = true;
               cout << "Mesh of cube complete and saved in file " << _mesh_file << endl;
               if (!ret)
                  *_ofh << "    save " << _mesh_file << endl;
               _saved = true;
               break;

            case 10:
            case 11:
               if (_verb>1)
                  cout << "Getting back to higher level ..." << endl;
               if (!_saved) {
                  if (_theMesh!=nullptr)
                     delete _theMesh, _theMesh = nullptr;
                  _theMesh = new OFELI::Mesh(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,cxmin,
                                             cxmax,cymin,cymax,czmin,czmax,HEXAHEDRON,size_t(_nb_dof));
                  _data->mesh_name.push_back("M"+to_string(_data->theMesh.size()));
                  _data->theMesh.push_back(_theMesh);
                  _saved = true;
               }
               *_ofh << "    end" << endl;
               if (_verb && !_saved)
                  cout << "Mesh 'cube' complete." << endl;
               _ret = 90;
               return;

            case 12:
            case 13:
               _ret = 100;
               return;

            case 14:
               _ret = 200;
               return;

            case -2:
               break;

            default:
               cout << "Unknown Command: " << _cmd->token() << endl;
               cout << "Available commands: min, max, ne, codes, nbdof, save, end, <" << endl;
               cout << "Global commands: help, ?, set, quit, exit" << endl;
               *_ofl << "In rita>mesh>cube>: Unknown command: " << _cmd->token() << endl;
               break;
         }
      }
   }
   _ret = 0;
   return;
}


void mesh::setCode()
{
   int nb=0, c=0, np=0, nc=0, ns=0, nv=0;
   int c_ok=0, points_ok=0, curves_ok=0, surfaces_ok=0, volumes_ok=0;
   vector<int> points, curves, surfaces, volumes;
   _ret = 0;
   if (_generator>0 && _generator<4) {
      cout << "Error: This keyword is not allowed for generated mesh." << endl;
      *_ofl << "In rita>mesh>code>: Keyword not allowed for generated mesh" << endl;
      return;
   }
   _generator = 10;
   if (_verb>1) {
      cout << "Default codes for generated nodes:    0" << endl;
      cout << "Default codes for generated elements: 1" << endl;
   }

   const vector<string> kw = {"help","?","set","value","points","curves","surfaces","volumes",
                              "end","<","quit","exit","EXIT"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case 3:
            c = _cmd->int_token(0);
            if (c<0)
               c = MY_RANDOM - c;
            c_ok++;
            break;

         case 4:
            np = nb;
            for (int j=0; j<np; ++j)
               points.push_back(_cmd->int_token(j));
            points_ok++;
            break;

         case 5:
            nc = nb;
            for (int j=0; j<nc; ++j)
               curves.push_back(_cmd->int_token(j));
            curves_ok++;
            break;

         case 6:
            ns = nb;
            for (int j=0; j<ns; ++j)
               surfaces.push_back(_cmd->int_token(j));
            surfaces_ok++;
            break;

         case 7:
            nv = nb;
            for (int j=0; j<nv; ++j)
               volumes.push_back(_cmd->int_token(j));
            volumes_ok++;
            break;

         default:
            cout << "Error: Unknown argument: " << _kw[n] << endl;
            *_ofl << "In rita>mesh>code>: Unknown argument: " << _kw[n] << endl;
            _ret = 1;
	    return;
      }
   }
   if (!c_ok) {
      cout << "Error: Missing code value." << endl;
      *_ofl << "In rita>mesh>code>: Missing code value." << endl;
      _ret = 1;
      return;
   }
   if (points_ok>1 || curves_ok>1 || surfaces_ok>1 || volumes_ok>1) {
      cout << "Error: Each entity must be given once only." << endl;
      *_ofl << "In rita>mesh>code>:  Each entity must be given once only." << endl;
      _ret = 1;
      return;
   }
   if (np+nc+ns+nv==0) {
      cout << "Error: At least one entity must be given." << endl;
      *_ofl << "In rita>mesh>code>: At least one entity must be given." << endl;
      _ret = 1;
      return;
   }
   *_ofh << "  code  value=" << c;
   if (points_ok) {
      *_ofh << "  points=";
      for (int j=0; j<np; ++j)
         _Pcode[c].l.push_back(points[j]), _Pcode[c].nb++;
      for (int j=0; j<np-1; ++j)
         *_ofh << points[j] << ",";
      *_ofh << points[np-1] << endl;
   }
   if (curves_ok) {
      *_ofh << "  curves=";
      for (int j=0; j<nc; ++j)
         _Ccode[c].l.push_back(curves[j]), _Ccode[c].nb++;
      for (int j=0; j<nc-1; ++j)
         *_ofh << curves[j] << ",";
      *_ofh << curves[nc-1] << endl;
   }
   if (surfaces_ok) {
      *_ofh << "  surfaces=";
      for (int j=0; j<ns; ++j)
         _Scode[c].l.push_back(surfaces[j]), _Scode[c].nb++;
      for (int j=0; j<ns-1; ++j)
         *_ofh << surfaces[j] << ",";
      *_ofh << surfaces[ns-1] << endl;
   }
   if (volumes_ok) {
      *_ofh << "  volumes=";
      for (int j=0; j<nv; ++j)
         _Vcode[c].l.push_back(volumes[j]), _Vcode[c].nb++;
      for (int j=0; j<nv-1; j++)
         *_ofh << volumes[j] << ",";
      *_ofh << volumes[nv-1] << endl;
   }
   _ret = 0;
   return;
}


void mesh::setPoint()
{
   _ret = 0;
   int nb=0;
   bool n_ok=false, x_ok=false, y_ok=false, z_ok=false, h_ok=false, del_ok=false;
   static int it = 0;
   if (it==0) {
      _point.n = 1;
      _point.x = _point.y = _point.z = 0.;
      _point.h = 0.1;
   }
   it++;
   const vector<string> kw = {"help","?","set","label","n","coord","size","end","<","quit","exit","EXIT"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case 3:
         case 4:
            _point.n = _cmd->int_token(0);
            n_ok = true;
            break;

         case 5:
            _point.x = _point.y = _point.z = 0.;
            if (nb==1) {
               _point.x = _cmd->double_token(0);
               x_ok = true;
            }
            else if (nb==2) {
              _point.x = _cmd->double_token(0);
              _point.y = _cmd->double_token(1);
               x_ok = y_ok = true;
            }
            else {
              _point.x = _cmd->double_token(0);
              _point.y = _cmd->double_token(1);
              _point.z = _cmd->double_token(2);
              x_ok = y_ok = z_ok = true;
            }
            break;

         case 6:
            _point.h = _cmd->double_token(0);
            h_ok = true;
            break;

         case 7:
            del_ok = true;
            break;

         default:
            cout << "Error: Unknown argument: " << _kw[n] << endl;
            *_ofl << "In rita>mesh>point>: Unknown argument: " << _kw[n] << endl;
            _ret = 1;
	    return;
      }
   }
   if (del_ok) {
      if (!n_ok) {
         cout << "Error: Point label must be given to be deleted." << endl;
         *_ofl << "In rita>mesh>point>: No point label given." << endl;
         _ret = 1;
         return;
      }
      _points.erase(_point.n);
      *_ofh << "  point  label="<< _point.n << "  delete" << endl;
      _ret = 0;
      return;
   }
   if (!n_ok) {
      _point.n++;
      if (_verb>1)
         cout << "Assumed point label: " << _point.n << endl;
   }
   if (_point.n<=0) {
      cout << "Error: Label must be positive." << endl;
      *_ofl << "In rita>mesh>point>: Label must be positive." << endl;
      _ret = 1;
      return;
   }
   if (!x_ok && _verb>1) {
      cout << "Warning: No x-coordinate given. Assumed x-coordinate: " << _point.x << endl;
      *_ofl << "In rita>mesh>point>: No x-coordinate given. Assumed x-coordinate: " << _point.x << endl;
   }
   if (!y_ok && _verb>1) {
      cout << "Warning: No y-coordinate given. Assumed y-coordinate: " << _point.y << endl;
      *_ofl << "In rita>mesh>point>: No y-coordinate given. Assumed y-coordinate: " << _point.y << endl;
   }
   if (!z_ok && _verb>1) {
      cout << "Warning: No z-coordinate given. Assumed z-coordinate: " << _point.z << endl;
      *_ofl << "In rita>mesh>point>: No z-coordinate given. Assumed z-coordinate: " << _point.z << endl;
   }
   if (!h_ok && _verb) {
      cout << "Warning: No mesh size given. Assumed mesh size: " << _point.h << endl;
      *_ofl << "In rita>mesh>point>: No mesh size given. Assumed mesh size: " << _point.h << endl;
   }
   if (_generator>0 && _generator<=3) {
      if (_generator==1) {
         cout << "Error: A 1-D mesh has already been generated.\n"
                 "       Retype command to delete older mesh and generate new one" << endl;
      }
      else if (_generator==2) {
         cout << "Error: A 2-D mesh has already been generated.\n"
                 "       Retype command to delete older mesh and generate new one" << endl;
      }
      else if (_generator==3) {
         cout << "Error: A 3-D mesh has already been generated.\n"
                 "       Retype command to delete older mesh and generate new one" << endl;
      }
      *_ofl << "In rita>mesh>point>: Mesh already generated. Needs retyping command to confirm.";
      _generator = 10;
   }
   _points[_point.n] = _point;
   *_ofh << "  point  label=" << _point.n << "  x=" << _point.x << "  y=" << _point.y
         << "  z=" << _point.z << "  size=" << _point.h << endl;
#ifndef USE_GMSH
   _theDomain->insertVertex(_point.x,_point.y,_point.z,_point.h,0);
#endif
   _ret = 0;
   return;
}


void mesh::setCurve()
{
   static int nn = 0;
   int nb=0, n1=0, n2=0, n3=0;
   bool n_ok=false, line_ok=false, circle_ok=false, del_ok=false;
   _ret = 0;
   if (_generator>0 && _generator<=3) {
      if (_generator==1) {
         cout << "Error: A 1-D mesh has already been generated.\n"
                 "       Retype command to delete older mesh and generate new one" << endl;
      }
      else if (_generator==2) {
         cout << "Error: A 2-D mesh has already been generated.\n"
                 "       Retype command to delete older mesh and generate new one" << endl;
      }
      else if (_generator==3) {
         cout << "Error: A 3-D mesh has already been generated.\n"
                 "       Retype command to delete older mesh and generate new one" << endl;
      }
      *_ofl << "In rita>mesh>line>: Mesh already generated. Needs retyping command to confirm.";
      _generator = 0;
      return;
   }
   _generator = 10;
   const vector<string> kw = {"help","?","set","label","n","line","circle","del$ete","end","<","quit","exit","EXIT"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case 3:
         case 4:
            nn = _cmd->int_token(0);
            n_ok = true;
            break;

         case 5:
            if (nb!=2) {
               cout << "Error: Illegal number of arguments for a line" << endl;
               *_ofl << "In rita>mesh>curve>: Illegal number of arguments for a line" << endl;
               _ret = 1;
               return;
            }
            n1 = _cmd->int_token(0);
            n2 = _cmd->int_token(1);
            line_ok = true;
            break;

         case 6:
            if (nb!=3) {
               cout << "Error: Illegal number of arguments for a circle" << endl;
               *_ofl << "In rita>mesh>curve>: Illegal number of arguments for a circle" << endl;
               _ret = 1;
               return;
            }
            n1 = _cmd->int_token(0);
            n2 = _cmd->int_token(1);
            n3 = _cmd->int_token(2);
            circle_ok = true;
            break;

         case 7:
            del_ok = true;
            break;

         default:
            cout << "Error: Unknown argument: " << _kw[n] << endl;
            *_ofl << "In rita>mesh>curve>: Unknown argument: " << _kw[n] << endl;
            _ret = 1;
	    return;
      }
   }

   if (del_ok) {
      if (!n_ok) {
         cout << "Error: Curve label must be given to be deleted." << endl;
         *_ofl << "In rita>mesh>curve>: No curve label given." << endl;
         _ret = 1;
         return;
      }
      _curve.erase(nn);
      *_ofh << "  curve  label="<< nn << "  delete" << endl;
      _ret = 0;
      return;
   }
   if (nn<=0) {
      cout << "Label must be positive." << endl;
      *_ofl << "In rita>mesh>curve>: Label must be positive." << endl;
      _ret = 1;
      return;
   }
   if (!n_ok) {
      nn++;
      if (_verb>1)
         cout << "Assumed curve label: " << nn << endl;
   }
   if (!line_ok && !circle_ok) {
      cout << "Error: A line or a circle must be defined." << endl;
      *_ofl << "In rita>mesh>curve>: A line or a circle must be defined." << endl;
      _ret = 1;
      return;
   }
   if (line_ok && circle_ok) {
      cout << "Error: Curve cannot be defined as a line and a circle simultaneously." << endl;
      *_ofl << "In rita>mesh>curve>: Curve cannot be defined as a line and a circle simultaneously." << endl;
      _ret = 1;
      return;
   }
   _curve[nn].nb = nn;
   *_ofh << "  curve label=" << nn;
   if (line_ok) {
      if (_points.find(n1)==_points.end()) {
         cout << "Error: Undefined point: " << n1 << endl;
         *_ofl << "In rita>mesh>curve>: Undefined end point: " << n1 << endl;
         _ret = 1;
         return;
      }
      if (_points.find(n2)==_points.end()) {
         cout << "Error: Undefined point: " << n2 << endl;
         *_ofl << "In rita>mesh>curve>: Undefined end point: " << n2 << endl;
         _ret = 1;
         return;
      }
      _curve[nn].l.clear(); _curve[nn].nb = 0;
      _curve[nn].l.push_back(n1); _curve[nn].nb++;
      _curve[nn].l.push_back(n2); _curve[nn].nb++;
      _curve[nn].type = 1;
      *_ofh << "  line=" << n1 << "," << n2 << endl;
   }
   if (circle_ok) {
      if (_points.find(n1)==_points.end()) {
         cout << "Error: Undefined point: " << n1 << endl;
         *_ofl << "In rita>mesh>curve>: Undefined point: " << n1 << endl;
         _ret = 1;
         return;
      }
      if (_points.find(n2)==_points.end()) {
         cout << "Error: Undefined point: " << n2 << endl;
         *_ofl << "In rita>mesh>curve>: Undefined point: " << n2 << endl;
         _ret = 1;
         return;
      }
      if (_points.find(n3)==_points.end()) {
         cout << "Error: Undefined point: " << n3 << endl;
         *_ofl << "In rita>mesh>curve>: Undefined point: " << n3 << endl;
         _ret = 1;
         return;
      }
      _curve[nn].l.clear();
      _curve[nn].l.push_back(n1);
      _curve[nn].l.push_back(n2);
      _curve[nn].l.push_back(n3);
      _curve[nn].nb = 3;
      _curve[nn].type = 2;
      *_ofh << "  circle=" << n1 << "," << n2 << "," << n3 << endl;
   }
#ifndef USE_GMSH
   _theDomain->insertLine(_curve[nn].l[0],_curve[nn].l[1],0);
#endif
   _ret = 0;
   return;
}


void mesh::setContour()
{
   int nb=0, nn=0, s=1, n1=0, n2=0;
   vector<int> curv, surf;
   int n_ok=0, curves_ok=0, surfaces_ok=0, del_ok=0;
   _ret = 0;
   if (_generator>0 && _generator<=3) {
      if (_generator==1) {
         cout << "Error: A 1-D mesh has already been generated.\n"
                 "       Retype command to delete older mesh and generate new one" << endl;
      }
      else if (_generator==2) {
         cout << "Error: A 2-D mesh has already been generated.\n"
                 "       Retype command to delete older mesh and generate new one" << endl;
      }
      else if (_generator==3) {
         cout << "Error: A 3-D mesh has already been generated.\n"
                 "       Retype command to delete older mesh and generate new one" << endl;
      }
      *_ofl << "In rita>mesh>contour>: Mesh already generated.";
      _generator = 0;
      return;
   }
   _generator = 10;
   const vector<string> kw = {"help","?","set","label","n","curv$es","surf$aces","end","<","quit","exit","EXIT"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case 3:
         case 4:
            nn = _cmd->int_token(0);
            n_ok++;
            break;

         case 5:
            for (int i=0; i<nb; i++)
               curv.push_back(_cmd->int_token(i));
            curves_ok++;
            break;

         case 6:
            for (int j=0; j<nb; j++)
               surf.push_back(_cmd->int_token(j));
            surfaces_ok++;
            break;

         case 7:
            del_ok++;
            break;

         default:
            cout << "Error: Unknown argument: " << _kw[n] << endl;
            *_ofl << "In rita>mesh>contour>: Unknown argument: " << _kw[n] << endl;
            _ret = 1;
	    return;
      }
   }
   if (!n_ok) {
      cout << "Error: Missing contour label." << endl;
      cout << "In rita>mesh>contour>: Missing contour label." << endl;
      _ret = 1;
      return;
   }
   if (n_ok>1) {
      cout << "Error: Contour label cannot be given more than once." << endl;
      cout << "In rita>mesh>contour>: Contour label cannot be given more than once." << endl;
      _ret = 1;
      return;
   }
   if (del_ok) {
      if (del_ok>1) {
         cout << "Error: Keyword delete cannot be given more than once." << endl;
         cout << "In rita>mesh>contour>: Keyword delete cannot be given more than once." << endl;
         _ret = 1;
         return;
      }
      _Ccontour.erase(nn);
      _Scontour.erase(nn);
      *_ofh << "  contour  label="<< nn << "  delete" << endl;
      _ret = 0;
      return;
   }
   if (nn==0) {
      cout << "Contour label cannot be zero." << endl;
      *_ofl << "In rita>mesh>contour>: Contour label cannot be zero." << endl;
      _ret = 1;
      return;
   }
   if (!curves_ok && !surfaces_ok) {
      cout << "Contour must be defined either by curves or surfaces." << endl;
      *_ofl << "In rita>mesh>contour>: Contour must be defined either by curves or surfaces." << endl;
      _ret = 1;
      return;
   }
   if (curves_ok && surfaces_ok) {
      cout << "Contour cannot be defined by curves or surfaces simultaneously." << endl;
      *_ofl << "In rita>mesh>contour>: Contour cannot be defined by curves or surfaces simultaneously." << endl;
      _ret = 1;
      return;
   }
   s = 1;
   if (nn<0)
      s = -1, nn = -nn;
   if (curves_ok) {
      if (curves_ok>1) {
         cout << "Error: A curve contour has already been defined." << endl;
         *_ofl << "In rita>mesh>contour>: A curve contour already defined." << endl;
         _ret = 1;
         return;
      }
      _Ccontour[nn].nb = 0;
      for (int i=0; i<nb; ++i) {
         int l = curv[i];
         if (_curve.find(l)==_curve.end()) {
            cout << "Error: Undefined curve: " << l << endl;
            *_ofl << "In rita>mesh>contour>: Undefined curve: " << l << endl;
            _ret = 1;
            return;
         }
         _Ccontour[nn].l.push_back(l), _Ccontour[nn].nb++;
      }
      n1 = _curve[_Ccontour[nn].l[0]].l[0];
      if (s==-1)
         n1 = _curve[_Ccontour[nn].l[0]].l[1];
      n2 = _curve[_Ccontour[nn].l[_Ccontour[nn].nb-1]].l[1];
      if (s==-1)
         n2 = _curve[_Ccontour[nn].l[_Ccontour[nn].nb-1]].l[0];
      if (n1!=n2) {
         cout << "Error: Non closed contour." << endl;
         *_ofl << "In rita>mesh>contour>: Non closed contour." << endl;
         _ret = 1;
         return;
      }	 
      *_ofh << "  contour  label=" << nn << "  curves=";
      for (int i=0; i<nb-1; ++i)
         *_ofh << curv[i] << ",";
      *_ofh << curv[nb-1] << endl;
      _nb_Ccontour = _Ccontour.size();
   }
   if (surfaces_ok) {
      if (surfaces_ok>1) {
         cout << "Error: A surface contour has already been defined." << endl;
         *_ofl << "In rita>mesh>contour>: A surface contour already defined." << endl;
      }
      _Scontour[nn].nb = 0;
      for (int i=0; i<nb; ++i) {
         int l = surf[i];
         if (_surface.find(l)==_surface.end()) {
            cout << "Error: Undefined surface: " << l << endl;
            *_ofl << "In rita>mesh>contour>: Undefined surface: " << l << endl;
            _ret = 1;
            return;
         }
         _Scontour[nn].l.push_back(l); _Scontour[nn].nb++;
      }
      n1 = _curve[_Scontour[nn].l[0]].l[0];
      if (s==-1)
         n1 = _surface[_Scontour[nn].l[0]].l[1];
      n2 = _surface[_Scontour[nn].l[_Scontour[nn].nb-1]].l[1];
      if (s==-1)
         n2 = _surface[_Scontour[nn].l[_Scontour[nn].nb-1]].l[0];
      if (n1 != n2) {
         cout << "Error: Non closed contour." << endl;
         *_ofl << "In rita>mesh>contour>: Non closed contour." << endl;
         _ret = 1;
         return;
      }
      *_ofh << "  contour  label=" << nn << "  surfaces=";
      for (int i=0; i<nb-1; ++i)
         *_ofh << surf[i] << ",";
      *_ofh << curv[nb-1] << endl;
      _nb_Scontour = _Scontour.size();
   }
   _ret = 0;
   return;
}


void mesh::setSurface()
{
   static int nn = 1;
   int nb=0, n=_surface.size()+1;
   bool n_ok=false, cont_ok=false, del_ok=false;
   vector<int> cont;
   _ret = 0;
   if (_generator>0 && _generator<=3) {
      if (_generator==1) {
         cout << "Error: A 1-D mesh has already been generated.\n"
                 "       Retype command to delete older mesh and generate new one" << endl;
      }
      else if (_generator==2) {
         cout << "Error: A 2-D mesh has already been generated.\n"
                 "       Retype command to delete older mesh and generate new one" << endl;
      }
      else if (_generator==3) {
         cout << "Error: A 3-D mesh has already been generated.\n"
                 "       Retype command to delete older mesh and generate new one" << endl;
      }
      *_ofl << "In rita>mesh>surface>: Mesh already generated. Needs retyping command to confirm.";
      _generator = 0;
      return;
   }
   _generator = 10;
   const vector<string> kw = {"help","?","set","label","n","contours","end","<","quit","exit","EXIT"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case 3:
         case 4:
            nn = _cmd->int_token(0);
            n_ok = true;
            break;

         case 5:
            for (int j=0; j<nb; ++j)
               cont.push_back(_cmd->int_token(j));
            cont_ok = true;
            break;

         case 6:
            del_ok = true;
            break;

         default:
            cout << "Error: Unknown argument: " << _kw[n] << endl;
            *_ofl << "In rita>mesh>surface>: Unknown argument: " << _kw[n] << endl;
            _ret = 1;
	    return;
      }
   }
   if (del_ok) {
      if (!n_ok) {
         cout << "Error: Surface label must be given to be deleted." << endl;
         *_ofl << "In rita>mesh>surface>: No surface label given." << endl;
         _ret = 1;
         return;
      }
      _surface.erase(nn);
      *_ofh << "  surface  label="<< nn << "  delete" << endl;
      _ret = 0;
      return;
   }
   if (!n_ok) {
      nn++;
      if (_verb)
         cout << "Assumed surface label: " << nn << endl;
   }
   if (nn==0) {
      cout << "Error: Label 0 is forbidden." << endl;
      *_ofl << "In rita>mesh>surface>: Label 0 is forbidden" << endl;
      _ret = 1;
      return;
   }
   if (!cont_ok) {
      cout << "Surface contours must be given." << endl;
      *_ofl << "In rita>mesh>surface>: Surface contours must be given." << endl;
      _ret = 1;
      return;
   }
   _surface[nn].nb = 0;
   for (int i=0; i<nb; ++i) {
      int c = cont[i];
      if (_Ccontour.find(c)==_Ccontour.end()) {
         cout << "Error: Undefined curve contour: " << n << endl;
         *_ofl << "In rita>mesh>surface>: Undefined curve contour: " << c << endl;
         _ret = 1;
         return;
      }
      _surface[nn].l.push_back(c); _surface[nn].nb++;
   }
   *_ofh << "  surface  label=" << nn << "  contours=";
   for (int i=0; i<nb-1; ++i)
      *_ofh << cont[i] << ",";
   *_ofh << cont[nb-1] << endl;
   _nb_surface = _surface.size();
   _ret = 0;
   return;
}


void mesh::setVolume()
{
   _ret = 0;
   _generator = 10;
   _nb_volume = 0;
}


void mesh::setSubDomain()
{
   _ret = 0;
   if (_generator>0 && _generator<=3) {
      if (_generator==1) {
         cout << "Error: A 1-D mesh has already been generated.\n"
                 "       Retype command to delete older mesh and generate new one" << endl;
      }
      else if (_generator==2) {
         cout << "Error: A 2-D mesh has already been generated.\n"
                 "       Retype command to delete older mesh and generate new one" << endl;
      }
      else if (_generator==3) {
         cout << "Error: A 3-D mesh has already been generated.\n"
                 "       Retype command to delete older mesh and generate new one" << endl;
      }
      *_ofl << "In rita>mesh>vertex>: Mesh already generated. Needs retyping command to confirm.";
      _generator = 0;
      return;
   }
   _generator = 10;
   int ret=0;
   *_ofh << "   subdomain " << endl;   
   _sd.ln = 1, _sd.orientation = 1, _sd.code = 1;
   if (_verb) {
      cout << "Default line in subdomain: 1\n";
      cout << "Default orientation: 1\n";
      cout << "Default code: 10" << endl;
   }

   const vector<string> kw = {"help","?","set","line","orient$ation","quit","exit","EXIT"};
   while (1) {
      if (_cmd->readline("rita>mesh>subdomain> ")<0)
         continue;
      switch (_key=_cmd->getKW(kw)) {

         case 0:
         case 1:
            _cmd->setNbArg(0);
            cout << "\nAvailable Commands\n";
            cout << "line        : Label for a line in subdomain\n";
            cout << "orientation : Orientation\n";
            cout << "code        : Code to associate to subdomain\n";
            cout << "save        : Save subdomain\n";
            cout << "end or <    : go back to mesh menu" << endl;
            break;

         case 2:
            if (_verb)
               cout << "Setting configuration parameter(s) ..." << endl;
            setConfigure();
            break;

         case 3:
            if (_verb>1)
               cout << "Setting line ..." << endl;
            if (_cmd->setNbArg(1,"Give label of a line in subdomain.")) {
               *_ofl << "In rita>mesh>subdomain>line>: Missing label of vertex in subdomain." << endl;
               break;
            }
            ret = _cmd->get(_sd.ln);
            if (!ret)
               *_ofh << "    line " << _sd.ln << endl;
            break;

         case 4:
            if (_verb>1)
               cout << "Setting orientation ..." << endl;
            if (_cmd->setNbArg(1,"Give orientation of subdomain (1/-1).")) {
               *_ofl << "In rita>mesh>subdomain>orientation>: Missing orientation." << endl;
               break;

            }
            ret = _cmd->get(_sd.orientation);
            if (!ret)
               *_ofh << "    orientation " << _sd.orientation << endl;   
            break;

         case 5:
            if (_verb>1)
               cout << "Setting code ..." << endl;
            if (_cmd->setNbArg(1,"Give code to associate to subdomain")) {
              *_ofl << "In rita>mesh>subdomain>code>: Missing code to associate to subdomain." << endl;
               break;
            }
            ret = _cmd->get(_sd.code);
            if (!ret)
               *_ofh << "    code " << _sd.code << endl;   
            break;

         case 6:
            if (_verb>1)
               cout << "Saving subdomain ..." << endl;
            _cmd->setNbArg(0);
    //        _theDomain->insertSubDomain(_sd.ln,_sd.orientation,_sd.code);
            _nb_sub_domain++;
            _subdomains.push_back(_sd);
            cout << "Subdomain Added." << endl;
            cout << "Line " << _sd.ln << ", Orientation " << _sd.orientation
                 << ", Code " << _sd.code << endl;
            cout << "Total number of subdomains: " << _nb_sub_domain << endl;
            *_ofh << "    save" << endl;
            _saved = true;
            return;

         case 7:
         case 8:
            if (_verb>1)
               cout << "Getting back to mesh menu ..." << endl;
            _cmd->setNbArg(0);
            if (!_saved) {
               cout << "Warning: Subdomain not saved. You can type again subdomain "
                    << "and save generated subdomain." << endl;
               *_ofl << "In rita>mesh>subdomain>: Subdomain not saved." << endl;
            }
            *_ofh << "    end" << endl;
            _ret = 0;
            return;

         case 9:
         case 10:
            _ret = 100;
            return;

         case 11:
            _ret = 200;
            return;

         case -2:
         case -3:
         case -4:
            break;

         default:
            cout << "Unknown Command " << _cmd->token() << endl;
            cout << "Available commands: lines, save, end, <" << endl;
            cout << "Global commands:    help, ?, set, quit, exit" << endl;
            *_ofl << "In rita>mesh>subdomain>: Unknown command " << _cmd->token() << endl;
            break;
       }
   }
   _ret = 0;
   return;
}


void mesh::saveGeo(const string& file)
{
   _generator = 10;
   ofstream s(file.c_str());
   for (auto const& v: _points)
      s << "Point(" << v.first << ") = {" << v.second.x << ", " << v.second.y
        << ", " << v.second.z << ", " << v.second.h << "};" << endl;

   for (auto const& v: _Pcode) {
      s << "Physical Point(" << v.first << ") = {";
      for (int i=0; i<v.second.nb-1; ++i)
         s << v.second.l[i] << ", ";
      s << v.second.l[v.second.nb-1] << "};" << endl;
   }
                   
   for (auto const& v: _curve) {
      Entity curve = v.second;
      if (curve.type==1)
         s << "Line(" << v.first << ") = {" << curve.l[0] << ", " << curve.l[1] << "};" << endl;
      else if (curve.type==2)
         s << "Circle(" << v.first << ") = {" << curve.l[0] << ", " << curve.l[2] << ", " 
           << curve.l[1] << "};" << endl;
   }

   for (auto const& v: _Ccode) {
      Entity code = v.second;
      s << "Physical Curve(" << v.first << ") = {";
      for (int i=0; i<code.nb-1; ++i)
         s << code.l[i] << ", ";
      s << code.l[code.nb-1] << "};" << endl;
   }

   for (auto const& v: _Ccontour) {
      s << "Curve Loop(" << v.first << ") = {";
      for (int i=0; i<v.second.nb-1; ++i)
         s << v.second.l[i] << ", ";
      s << v.second.l[v.second.nb-1] << "};" << endl;
   }

   for (auto const& v: _surface) {
      s << "Plane Surface(" << v.first << ") = {";
      for (int i=0; i<v.second.nb-1; ++i)
         s << v.second.l[i] << ", ";
      s << v.second.l[v.second.nb-1] << "};" << endl;
   }

   for (auto const& v: _Scode) {
      s << "Physical Surface(" << v.first << ") = {";
      for (int i=0; i<v.second.nb-1; ++i)
         s << v.second.l[i] << ", ";
      s << v.second.l[v.second.nb-1] << "};" << endl;
   }

   for (auto const& v: _Vcode) {
      s << "Physical Volume(" << v.first << ") = {";
      for (int i=0; i<v.second.nb-1; ++i)
         s << v.second.l[i] << ", ";
      s << v.second.l[v.second.nb-1] << "};" << endl;
   }
   _geo = true;
}


void mesh::saveDomain(const string& file)
{
   ofstream id(file.c_str());
   id << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\" ?>\n<OFELI_File>\n<info>" << endl;
   id << "<title>Domain file created by rita</title>" << endl;
   id << "<date></date>\n<author></author>\n</info>" << endl;
   id << "<Project name=\"\">" << endl;
   id << "  <Domain dim=\"2\">" << endl;
   for (int i=0; i<_nb_point; ++i)
      id << "    <vertex>" << _points[i].x << " " << _points[i].y << " "
         << 0 << " " << _points[i].h << "</vertex>" << endl;
   for (int i=1; i<=_nb_curve; ++i)
      id << "    <line>" << _curve[i].l[0] << " " << _curve[i].l[1] << " "
         << 0 << "</line>" << endl;
/*   for (int i=0; i<_nb_circle; ++i)
      id << "    <circle>" << _circles[i].pt1 << " " << _circles[i].c << " " << _circles[i].pt2
         << " " << 0 << "</circle>" << endl;*/
   for (int i=0; i<_nb_sub_domain; ++i)
      id << "    <Subdomain>" << _subdomains[i].ln << " " << _subdomains[i].orientation << " "
         << _subdomains[i].code << "</Subdomain>" << endl;
   id << "  </Domain>\n</Project>\n</OFELI_File>" << endl;
}


void mesh::Generate()
{
   _ret = 0;
   if (_generator>0 && _generator<=3) {
      if (_generator==1) {
         cout << "A 1-D mesh has already been generated." << endl;
         *_ofl << "In rita>mesh>generate>: 1-D mesh already generated" << endl;
      }
      else if (_generator==2) {
         cout << "A 2-D mesh has already been generated." << endl;
         *_ofl << "In rita>mesh>generate>: 2-D mesh already generated" << endl;
      }
      else if (_generator==3) {
         cout << "A 3-D mesh has already been generated." << endl;
         *_ofl << "In rita>mesh>generate>: 3-D mesh already generated" << endl;
      }
      return;
   }
   _generator = 10;
   if (_theMesh!=nullptr)
      delete _theMesh, _theMesh = nullptr;
   _mesh_file = "rita.m";
   
// Save geo gmsh file and generate gmsh mesh
   if (!_geo)
      saveGeo("rita.geo");

#ifdef USE_GMSH
   if (_verb)
      cout << "Starting mesh generation using Gmsh ..." << endl;
   gmsh::initialize();
   gmsh::option::setNumber("General.Terminal",1);
   gmsh::model::add("rita");
   for (auto const& v: _points)
      gmsh::model::geo::addPoint(v.second.x,v.second.y,v.second.z,v.second.h,v.first);
   for (auto const& v: _curve) {
      if (v.second.type==1)
         gmsh::model::geo::addLine(v.second.l[0],v.second.l[1],v.first);
      else if (v.second.type==2)
         gmsh::model::geo::addCircleArc(v.second.l[0],v.second.l[2],v.second.l[1],v.first);
   }
   for (auto const& v: _Ccontour)
      gmsh::model::geo::addCurveLoop(v.second.l,v.first);
   for (auto const& v: _surface)
      gmsh::model::geo::addPlaneSurface(v.second.l,v.first);
   for (auto const& v: _Pcode)
      gmsh::model::addPhysicalGroup(0,v.second.l,v.first);
   for (auto const& v: _Ccode)
      gmsh::model::addPhysicalGroup(1,v.second.l,v.first);
   for (auto const& v: _Scode)
      gmsh::model::addPhysicalGroup(2,v.second.l,v.first);
   for (auto const& v: _Vcode)
      gmsh::model::addPhysicalGroup(3,v.second.l,v.first);
   gmsh::model::setPhysicalName(2,4,"Domain");
   gmsh::model::geo::synchronize();
   gmsh::model::mesh::generate(2);
   gmsh::write("rita.msh");
   gmsh::finalize();
   if (_verb)
      cout << "Gmsh mesh generation complete." << endl;
   _theMesh = new OFELI::Mesh("rita.msh",false,NODE_DOF,_nb_dof);
   _generated = true;
#else
   _theDomain->setDim(_dim);
   _theDomain->setNbDOF(size_t(_nb_dof));
   saveDomain("rita.dom");
   _theDomain->genMesh(_mesh_file);
   _theMesh = new OFELI::Mesh(_mesh_file,false,NODE_DOF,2);
   _generated = true;
   _generator = 4;
#endif
   *_ofh << "  generate" << endl;
   _data->mesh_name.push_back("M"+to_string(_data->theMesh.size()));
   _data->theMesh.push_back(_theMesh);
   _ret = 0;
   return;
}


void mesh::setNbDOF()
{
   int n = 0;
   if (_cmd->setNbArg(1,"Give number of degrees of freedom per node.")) {
      *_ofl << "In rita>mesh>nbdof>: Missing number of degrees of freedom." << endl;
      _ret = 1;
      return;
   }
   _ret = _cmd->get(n);
   if (!_ret) {
      *_ofh << "    nbdof " << n << endl;
      _nb_dof = n;
   }
}


void mesh::Plot()
{
   if (_theMesh==nullptr) {
      cout << "No mesh saved." << endl;
      *_ofl << "In rita>mesh>plot>: No mesh to plot." << endl;
      return;
   }
   _cmd->setNbArg(0);
   string file = "rita.msh";
   _theMesh->save(file);
   string com = "gmsh " + file;
   if (system(com.c_str())) {
      cout << "Error in system command." << endl;
      *_ofl << "In rita>mesh>plot>: Error in system command." << endl;
      _ret = 1;
      return;
   }
   remove(file.c_str());
   _ret = 0;
}


void mesh::Read()
{
   int ret=0;
   ifstream ip;
   string dom_file, file, bamg_file, geo_file, out_file, Cmd, msh_file;
   int mesh_ok=0, geo_ok=0, gmsh_ok=0;
   const vector<string> kw = {"help","?","set","mesh","geo","gmsh","end","<","quit","exit","EXIT"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg();
      switch (n) {

         case 3:
            file = _cmd->string_token();
            mesh_ok++;
            break;

         case 4:
            file = _cmd->string_token();
            geo_ok++;
            break;

         case 5:
            file = _cmd->string_token();
            gmsh_ok++;
            break;

         default:
            cout << "Error: Unknown argument: " << _kw[n] << endl;
            *_ofl << "In rita>mesh>read>: Unknown argument: " << _kw[n] << endl;
            _ret = 1;
            return;
      }
   }
   if (nb_args>0) {
      if (mesh_ok+geo_ok+gmsh_ok==0) {
         cout << "Error: One input file must be provided." << endl;
         *_ofl << "In rita>mesh>read>: One input file must be provided." << endl;
         _ret = 1;
         return;
      }
      if (mesh_ok+geo_ok+gmsh_ok>1) {
         cout << "Error: Only one input file must be provided." << endl;
         *_ofl << "In rita>mesh>read>: Only one input file must be provided." << endl;
         _ret = 1;
         return;
      }
      ip.open(file);
      if (ip.is_open())
         ip.close();
      else {
         cout << "Error: Unable to open file: " << file << endl;
         *_ofl << "In rita>mesh>read>: Unable to open file: " << file << endl;
         _ret = 1;
         return;
      }
      *_ofh << "  read";
      if (mesh_ok) {
         if (file.substr(file.find_last_of(".")+1)!="m") {
            cout << "Error: File extension must be \".m\"" << endl;
            *_ofl << "In rita>mesh>read>: File extension must be \".m\"";
            _ret = 1;
            return;
         }
         _theMesh = new OFELI::Mesh(file);
         if (_theMesh->getNbNodes()==0) {
            cout << "Error: Empty mesh." << endl;
            *_ofl << "In rita>mesh>read>: Empty mesh";
            _ret = 1;
            return;
         }
         _mesh_file = file;
         _nb_dof = _theMesh->getNbDOF() / _theMesh->getNbNodes();
         *_ofh << "  mesh=" << file;
         _data->mesh_name.push_back("M"+to_string(_data->theMesh.size()));
         _data->theMesh.push_back(_theMesh);
         _generated = true;
      }
      else if (geo_ok) {
         if (file.substr(file.find_last_of(".")+1)!="geo") {
            cout << "Error: File extension must be \".geo\"" << endl;
            *_ofl << "In rita>mesh>read>: File extension must be \".geo\"";
            _ret = 1;
            return;
         }
         msh_file = file.substr(0,file.rfind(".")) + "msh";
         Cmd = "gmsh -2 " + file + " -o " + msh_file;
         if (system(Cmd.c_str())) {
            cout << "Error in system command." << endl;
            *_ofl << "In rita>mesh>read>: Error in system command.";
            _ret = 1;
            return;
         }
         _theMesh = new OFELI::Mesh;
         _theMesh->get(file,GMSH);
         _data->mesh_name.push_back("M"+to_string(_data->theMesh.size()));
         _data->theMesh.push_back(_theMesh);
         _generated = true;
         *_ofh << "  geo=" << file;
      }
      else if (gmsh_ok) {
         if (file.substr(file.find_last_of(".")+1)!="msh") {
            cout << "Error: File extension must be \".msh\"" << endl;
            *_ofl << "In rita>mesh>read>: File extension must be \".msh\"";
            _ret = 1;
            return;
         }
         _generated = true;
         _theMesh = new OFELI::Mesh;
         _theMesh->get(file,GMSH);
         _data->mesh_name.push_back("M"+to_string(_data->theMesh.size()));
         _data->theMesh.push_back(_theMesh);
         *_ofh << "  gmsh=" << file;
      }
      *_ofh << endl;
   }
   else {
      _cmd->setNbArg(0);
      while (1) {
         if (_cmd->readline("rita>mesh>read> ")<0)
            continue;

         switch (_key=_cmd->getKW(kw)) {

            case 0:
            case 1:
               _cmd->setNbArg(0);
               cout << "\nAvailable Commands\n";
               cout << "mesh     : Read mesh in OFELI mesh file\n";
               cout << "geo      : Read mesh in OFELI mesh file\n";
               cout << "gmsh     : Read mesh in gmsh file" << endl;
               cout << "< or end : Go back to mesh menu" << endl;
               break;

            case 2:
               setConfigure();
               break;

            case 3:
               if (_cmd->setNbArg(1,"Give OFELI mesh file name.")) {
                  *_ofl << "In rita>mesh>read>mesh>: Missing Mesh file name." << endl;
                  break;
               }
               ret = _cmd->get(file);
               ip.open(file);
               if (ip.is_open()) {
                  ip.close();
                  _theMesh = new OFELI::Mesh(file);
                  *_ofh << "  read mesh " << file << endl;
                  if (_theMesh->getNbNodes()==0) {
                     cout << "Error: Empty mesh." << endl;
                     *_ofl << "In rita>mesh>read>mesh>: Empty mesh";
                     _ret = 1;
                     return;
                  }
                  _nb_dof = _theMesh->getNbDOF() / _theMesh->getNbNodes();
                  _data->mesh_name.push_back("M"+to_string(_data->theMesh.size()));
                  _data->theMesh.push_back(_theMesh);
                  _generated = true;
                  _ret = 90;
               }
               else {
                  cout << "Unable to open file: " << file << endl;
                  *_ofl << "In rita>mesh>read>mesh>: Unable to open file: " << file << endl;
                  _generated = false;
                  _ret = 0;
               }
               return;

            case 4:
               if (_cmd->setNbArg(1,"Give geo gmsh file name.")) {
                  *_ofl << "In rita>mesh>read>geo>: Missing geo file name." << endl;
                  break;
               }
               ret = _cmd->get(file);
               if (!ret) {
                  if (file.substr(file.find_last_of(".")+1)!=".geo") {
                     cout << "Error: File extension must be \".geo\"" << endl;
                     *_ofl << "In rita>mesh>read>geo>: File extension must be \".geo\"";
                     _ret = 1;
                     break;
                  }
                  msh_file = file.substr(0,file.rfind(".")) + "msh";
                  Cmd = "gmsh -2 " + file + " -o " + msh_file;
                  if (system(Cmd.c_str())) {
                     cout << "Error in system command." << endl;
                     *_ofl << "In rita>mesh>read>geo>: Error in system command.";
                     _ret = 1;
                     break;
                  }
                  _theMesh = new OFELI::Mesh;
                  _theMesh->get(msh_file,GMSH);
                  *_ofh << "  read geo " << file << endl;
                  _data->mesh_name.push_back("M"+to_string(_data->theMesh.size()));
                  _data->theMesh.push_back(_theMesh);
                  _geo = true;
                  _generated = false;
	       }
               _ret = 90;
               return;

            case 5:
               if (_cmd->setNbArg(1,"Give gmsh file name where to save mesh.")) {
                  *_ofl << "In rita>mesh>read>gmsh>: Missing gmsh file name." << endl;
                  break;
               }
               ret = _cmd->get(file);
               if (!ret) {
                  if (file.substr(file.find_last_of(".")+1)!=".msh") {
                     cout << "Error: File extension must be \".msh\"" << endl;
                     *_ofl << "In rita>mesh>read>gmsh>: File extension must be \".msh\"";
                     _ret = 1;
                     break;
                  }
                  _theMesh = new OFELI::Mesh;
                  _theMesh->get(file,GMSH);
                  _data->mesh_name.push_back("M"+to_string(_data->theMesh.size()));
                  _data->theMesh.push_back(_theMesh);
                  *_ofh << "  read gmsh " << file << endl;
	       }
               _generated = true;
               _ret = 90;
               return;

            case 6:
            case 7:
               _ret = 0;
               return;

            case 8:
            case 9:
               _ret = 100;
               return;

            case 10:
               _ret = 200;
               return;

            case -2:
            case -3:
            case -4:
               break;

            default:
               cout << "Unknown Command " << _cmd->token() << endl;
               cout << "Available commands: mesh, geo, gmsh" << endl;
               cout << "Global commands: help, ?, set, quit, exit" << endl;
               *_ofl << "In rita>mesh>read>: Unknown command " << _cmd->token() << endl;
               break;
         }
      }
   }
   _ret = 0;
   return;
}


void mesh::Save()
{
   string domain_f="rita.dom", geo_f="rita.geo", mesh_f="rita.m", gmsh_f="rita.msh";
   string t="", vtk_f="rita.vtk", gnuplot_f="rita-gpl.dat", matlab_f="rita-matlab.m";
   string tecplot_f="rita-tecplot.dat";
   int domain_ok=0, geo_ok=0, mesh_ok=0, gmsh_ok=0, vtk_ok=0, gnuplot_ok=0, matlab_ok=0, tecplot_ok=0;
   _ret = 0;
   const vector<string> kw = {"help","?","set","domain","geo","mesh","gmsh","vtk","gnuplot",
                              "matlab","tecplot","end","<","quit","exit","EXIT"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg();
      switch (n) {

         case 3:
#ifndef USE_GMSH
            if ((t=_cmd->string_token())!="")
               domain_f = t;
            domain_ok++;
#endif
            break;

         case 4:
            if ((t=_cmd->string_token())!="")
               geo_f = t;
            geo_ok++;
            break;

         case 5:
            if ((t=_cmd->string_token())!="")
               mesh_f = t;
            mesh_ok++;
            break;

         case 6:
            if ((t=_cmd->string_token())!="")
               gmsh_f = t;
            gmsh_ok++;
            break;

         case 7:
            if ((t=_cmd->string_token())!="")
               vtk_f = t;
            vtk_ok++;
            break;

         case 8:
            if ((t=_cmd->string_token())!="")
               gnuplot_f = t;
            gnuplot_ok++;
            break;

         case 9:
            if ((t=_cmd->string_token())!="")
               matlab_f = t;
            matlab_ok++;
            break;

         case 10:
            if ((t=_cmd->string_token())!="")
               tecplot_f = t;
            tecplot_ok++;
            break;

         default:
            cout << "Error: Unknown argument." << endl;
            *_ofl << "In rita>mesh>save>: Unknown argument." << endl;
            _ret = 1;
	    break;
      }
   }
   if (domain_ok+geo_ok+mesh_ok+gmsh_ok+vtk_ok+gnuplot_ok+matlab_ok+tecplot_ok==0) {
      cout << "Error: Nothing to save." << endl;
      cout << "In rita>mesh>save>: Nothing to save." << endl;
      _ret = 1;
      return;
   }
   if (mesh_ok+gmsh_ok+vtk_ok+gnuplot_ok+matlab_ok+tecplot_ok>0 && !_generated) {
      cout << "Error: No generated mesh to be saved." << endl;
      *_ofl << "In rita>mesh>save>: No generated mesh to be saved" << endl;
      _ret = 1;
      return;
   }
   *_ofh << "  save";
   if (domain_ok) {
      saveDomain(domain_f);
      *_ofh << "  domain=" << domain_f;
   }
   if (geo_ok) {
      saveGeo(geo_f);
      *_ofh << "  geo=" << geo_f;
   }
   if (mesh_ok) {
      _mesh_file = mesh_f;
      _theMesh->put(mesh_f);
      *_ofh << "  mesh=" << mesh_f;
   }
   if (gmsh_ok) {
      saveMesh(gmsh_f,*_theMesh,GMSH);
      *_ofh << "  gmsh=" << gmsh_f;
   }
   if (vtk_ok) {
      saveMesh(vtk_f,*_theMesh,VTK);
      *_ofh << "  vtk=" << vtk_f;
   }
   if (gnuplot_ok) {
      saveMesh(gnuplot_f,*_theMesh,GNUPLOT);
      *_ofh << "  gnuplot=" << gnuplot_f;
   }
   if (matlab_ok) {
      saveMesh(matlab_f,*_theMesh,MATLAB);
      *_ofh << "  matlab=" << matlab_f;
   }
   if (tecplot_ok) {
      saveMesh(tecplot_f,*_theMesh,TECPLOT);
      *_ofh << "  tecplot=" << tecplot_f;
   }
   *_ofh << endl;
   return;
}

 /*
void mesh::Save()
{
   int ret=0;
   string file;
   *_ofh << "  save" << endl;
   while (1) {
      if (_cmd->readline("rita>mesh>save> ")<0)
         continue;

      switch (_key=_cmd->getKW(_kw_save)) {

         case 0:
         case 1:
            _cmd->setNbArg(0);
            cout << "\nAvailable Commands\n";
#ifndef USE_GMSH
            cout << "domain  : Save Domain in a Domain format file\n";
#endif
            cout << "geo     : Save geometry in gmsh geo file\n";
            cout << "mesh    : Save mesh in OFELI mesh file\n";
            cout << "gmsh    : Save mesh in gmsh file\n";
            cout << "vtk     : Save mesh in vtk file\n";
            cout << "gnuplot : Save mesh in gnuplot file\n";
            cout << "matlab  : Save mesh in Matlab(R) file\n";
            cout << "tecplot : Save mesh in Tecplot(R) graphics file" << endl;
            break;

         case 2:
            setConfigure();
            break;

         case 3:
#ifndef USE_GMSH
            file = "rita.dom";
            if (_cmd->getNbArgs()==1)
               ret = _cmd->get(file);
            if (!ret)
               saveDomain(file);
            *_ofh << "    domain " << file << endl;
            _ret = 0;
#endif
            return;

         case 4:
            file = "rita.geo";
            if (_cmd->getNbArgs()==1)
               ret = _cmd->get(file);
            if (!ret)
               *_ofh << "    geo " << file << endl;
            saveGeo(file);
            _ret = 0;
            return;

         case 5:
            _mesh_file = "rita.m";
            if (_cmd->getNbArgs()==1)
               ret = _cmd->get(_mesh_file);
            if (_generated==false) {
               cout << "Error: No mesh generated to be saved." << endl;
               *_ofl << "In rita>mesh>save>mesh>: No mesh generated to be saved" << endl;
               return;
            }
            _theMesh->put(_mesh_file);
            if (!ret)
               *_ofh << "    mesh " << _mesh_file << endl;
            if (_verb)
               cout << "Mesh saved in file: " << _mesh_file << endl;
            _ret = 90;
            return;

         case 6:
            file = "rita.msh";
            if (_cmd->getNbArgs()==1)
               ret = _cmd->get(file);
            if (_generated==false) {
               cout << "Error: No mesh generated to be saved." << endl;
               *_ofl << "In rita>mesh>save>msh>: No mesh generated to be saved" << endl;
               return;
            }
            saveMesh(file,*_theMesh,VTK);
            if (!ret)
               *_ofh << "    gmsh " << file << endl;
            if (_verb)
               cout << "Mesh saved in file: " << file << endl;
            _ret = 90;
            return;

         case 7:
            if (_generated==false) {
               cout << "Error: No mesh generated." << endl;
               *_ofl << "In rita>mesh>save>vtk>: No mesh generated." << endl;
               _ret = 1;
               return;
            }
            file = "rita.vtk";
            if (_cmd->getNbArgs()==1)
               ret = _cmd->get(file);
            saveMesh(file,*_theMesh,VTK);
            if (!ret)
               *_ofh << "    vtk " << file << endl;
            _ret = 0;
            return;

         case 8:
            if (_generated==false) {
               cout << "Error: No mesh generated." << endl;
               *_ofl << "In rita>mesh>save>gnuplot>: No mesh generated." << endl;
               _ret = 1;
               return;
            }
            file = "rita-gpl.dat";
            if (_cmd->getNbArgs()==1)
               ret = _cmd->get(file);
            saveMesh(file,*_theMesh,GNUPLOT);
            if (!ret)
               *_ofh << "    gnuplot " << file << endl;
            _ret = 0;
            return;

         case 9:
            if (_generated==false) {
               cout << "Error: No mesh generated." << endl;
               *_ofl << "In rita>mesh>save>matlab>: No mesh generated." << endl;
               _ret = 1;
               return;
            }
            file = "rita-matlab.m";
            if (_cmd->getNbArgs()==1)
               ret = _cmd->get(file);
            saveMesh(file,*_theMesh,MATLAB);
            if (!ret)
               *_ofh << "    matlab " << file << endl;
            _ret = 0;
            return;

         case 10:
            if (_generated==false) {
               cout << "Error: No mesh generated." << endl;
               *_ofl << "In rita>mesh>save>tecplot>: No mesh generated." << endl;
               _ret = 1;
               return;
            }
            file = "rita-tecplot.dat";
            if (_cmd->getNbArgs()==1)
               ret = _cmd->get(file);
            saveMesh(file,*_theMesh,TECPLOT);
            if (!ret)
               *_ofh << "    tecplot " << file << endl;
            _ret = 0;
            return;

         case 11:
         case 12:
            _cmd->setNbArg(0);
            if (!_saved) {
               cout << "Warning: Nothing saved." << endl;
               *_ofl << "In rita>mesh>save>: Nothing saved." << endl;
            }
            *_ofh << "    end" << endl;
            _ret = 0;
            return;

         case 13:
         case 14:
            _ret = 100;
            return;

         case 15:
            _ret = 200;
            return;

         case -2:
         case -3:
         case -4:
            break;

         default:
            cout << "Unknown Command " << _cmd->token() << endl;
#ifdef USE_GMSH
            cout << "Available commands: geo, mesh, gmsh, vtk, gnuplot, matlab, tecplot, end, <" << endl;
#else
            cout << "Available commands: domain, geo, mesh, gmsh, vtk, gnuplot, matlab, tecplot, end, <" << endl;
#endif
            cout << "Global commands:    help, ?, set, quit, exit" << endl;
            *_ofl << "In rita>mesh>save>: Unknown command " << _cmd->token() << endl;
            break;
       }
   }

   if (_generator < 4) {
      cout << "Mesh already saved in file." << endl;
      *_ofl << "In rita>mesh>save>: Mesh already saved in file." << endl;
      _ret = 1;
      return;
   }
   if (_cmd->setNbArg(1,"Give name to file where to save."))
      return;
   if (_theMesh==nullptr) {
      cout << "Mesh was not properly created." << endl;
      *_ofl << "In rita>mesh>save>: mesh was not properly created" << endl;
      _ret = 1;
      return;
   }
   _cmd->get(file);
   _theMesh->save(file);
   _ret = 0;
   return;
}*/


void mesh::Clear()
{
   _cmd->setNbArg(0);
   if (_theMesh!=nullptr)
      delete _theMesh;
   _theMesh = nullptr;
   _saved = false;
   _dim = 1;
   _ret = 0;
/*   if (_theDomain)
      delete _theDomain;
   _theDomain = nullptr;*/
   _nb_point = _nb_curve = _nb_sub_domain = 0;
   _points.clear();
   _curve.clear();
   _surface.clear();
   _volume.clear();
   _subdomains.clear();
   *_ofh << "  clear" << endl;
   if (_verb)
      cout << "Mesh cleared." << endl;
}

} /* namespace RITA */
