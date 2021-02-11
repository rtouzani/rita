#!/bin/sh

# test_rita.sh 
# A script to test rita

echo "==============================================="
echo "Testing rita Tutorial ..."
echo "==============================================="


cd tutorial/ae

echo "-----------------------------------------------"
echo "Test solution of algebraic equations (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
../../src/rita example1.rita
../../src/rita example2.rita
../../src/rita example3.rita
fi

echo "-----------------------------------------------------------"
echo "Test solution of ordinary differential equations (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
cd ../ode
../../src/rita example1.rita
../../src/rita example2.rita
fi

echo "----------------------------------------------------------"
echo "Test solution of partial differential equations (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
cd ../pde
../../src/rita example1.rita
../../src/rita example2.rita
../../src/rita example3.rita
../../src/rita example4.rita
../../src/rita example5.rita
fi

echo "-------------------------------------------------"
echo "Test solution of optimization problems (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
cd ../optim
../../src/rita example1.rita
../../src/rita example2.rita
../../src/rita example3.rita
../../src/rita example4.rita
fi

echo "-------------------------------------"
echo "Test numerical integration (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
cd ../integration
../../src/rita example1.rita
../../src/rita example2.rita
fi

echo "-----------------------------------"
echo "Remove all created files (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
echo "Cleaning ..."
cd ..
make clean
fi

cd ../..
echo "========================================================="
echo "rita Testing complete"
