#include <stdlib.h>
#include <math.h>
#include <iostream>

int main()
{
   double RandomFloat(double a, double b);
   size_t n=500;
   double dx = 1./n, x=0., y=0.;
   for (size_t i=0; i<=n; ++i) {
      y = sin(3.1416*x) + RandomFloat(-0.1,0.1);
      if (i%25==0)
         std::cout << x << "  " << y << std::endl;
      x += dx;
   }
   return 0;
}


double RandomFloat(double a, double b)
{
   double random = ((double) rand()) / (double) RAND_MAX;
   return a + random*(b-a);
}
