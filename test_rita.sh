#!/bin/sh

# test_rita.sh 
# A script to test rita

echo "==============================================="
echo "Testing rita Tutorial ..."
echo "==============================================="

DIR="/usr/local"
DD=$DIR/share/rita
RITA=$DIR/bin/rita

echo "-----------------------------------------------"
echo "Test solution of algebraic equations (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   $RITA ${DD}/tutorial/ae/example1.rita
   $RITA ${DD}/tutorial/ae/example2.rita
   $RITA ${DD}/tutorial/ae/example3.rita
fi

echo "-----------------------------------------------------------"
echo "Test solution of ordinary differential equations (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   $RITA ${DD}/tutorial/ode/example1.rita
   $RITA ${DD}/tutorial/ode/example2.rita
fi

echo "----------------------------------------------------------"
echo "Test solution of partial differential equations (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   $RITA ${DD}/tutorial/pde/example1.rita
   $RITA ${DD}/tutorial/pde/example2.rita
   $RITA ${DD}/tutorial/pde/example3.rita
   $RITA ${DD}/tutorial/pde/example4.rita
   $RITA ${DD}/tutorial/pde/example5.rita
fi

echo "-------------------------------------------------"
echo "Test solution of optimization problems (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   $RITA ${DD}/tutorial/optim/example1.rita
   $RITA ${DD}/tutorial/optim/example2.rita
   $RITA ${DD}/tutorial/optim/example3.rita
   $RITA ${DD}/tutorial/optim/example4.rita
fi

echo "-------------------------------------"
echo "Test numerical integration (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   $RITA ${DD}/tutorial/integration/example1.rita
   $RITA ${DD}/tutorial/integration/example2.rita
fi

echo "-----------------------------------"
echo "Cleaning ..."
rm -f *.vtk *.pos *.sol *.msh *.geo rita-1d.m example3.m example5.m example1.dat ex2*
rm -f .rita.his .rita.log

echo "========================================================="
echo "rita Testing complete"
