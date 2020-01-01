/*NAME: Henry Unruh
 *CLASS: CS330, Numerical Computing
 *ASSIGNMENT: Estimating pi
 *SYNOPSIS: The goal of the assignment is to estimate the value of pi using four
 *different methods: the trapezoid, Simpson's (1/3)rd and (3/8)th's, and Boole's
 rules. 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//Pi inclusion statement, from PDF
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

//integral estimation, given from assignment
long double integral(long double x){
long double i;
i = 4.0/(1.0 + x*x);
return i;
}

//Function declarations
long double trapezoid();
long double simps13();
long double simps38();
long double booles();

//Trapezoid rule
long double trapezoid(long double x, long double y, int num)
{
long double est;
long double h;
long double c1;
long double c2;
long double tErr;
int i = 0;

est = 0;
h = (y-x)/num;
c1 = 0;
c2 = h;

  while(i < num){
  est = est + ((0.5*h)*(integral(c1) + integral(c2)));
   c1 = c2;
   c2 = (c2+h);
   i++;
  }

return tErr = (fabsl(M_PI - est));
}

//Simpson's (1/3)rd rule
long double simps13(long double x, long double y, int num)
{
long double est;
long double h;
long double c1;
long double c2;
long double c3;
long double s13Err;
int i = 0;

est = 0;
h = (y-x)/num;
c1 = 0;
c2 = h;
c3 = h+h;

  while(i < (num/2)){
  est = est + ((h/3)*(integral(c1) + (4*integral(c2))+integral(c3)));
   c1 = c1 + (h*2);
   c2 = c2 + (h*2);
   c3 = c3 + (h*2);
   i++;
  }

return s13Err = (fabsl(M_PI - est));
}

//Simpson's (3/8)th rule
long double simps38(long double x, long double y, int num)
{
long double est;
long double h;
long double c1;
long double c2;
long double c3;
long double c4;
long double s38Err;
int i = 0;

est = 0;
h = (y-x)/num;
c1 = 0;
c2 = h;
c3 = h*2;
c4 = h*3;

  while(i < num) {
  est = est + (3*h/8)*(integral(c1) + (3*integral(c2)) + (3*integral(c3)) + integral(c4)); 
   c1 = c1 + (h*4);
   c2 = c2 + (h*4);
   c3 = c3 + (h*4);
   c4 = c4 + (h*4);
   i = (i+3);
  }
return s38Err = (fabsl(M_PI - est));
}

//Boole's rule
long double booles(long double x, long double y, int num)
{
long double est;
long double h;
long double c1;
long double c2;
long double c3;
long double c4;
long double c5;
long double booleErr;
int i=0;

est = 0;
h = (y-x)/num;
c1 = 0;
c2 = h;
c3 = (h*2);
c4 = (h*3);
c5 = (h*4);

   while(i < num) {
   est = est + (2*h/45)*(7*integral(c1) + 32*integral(c2) + 12*integral(c3) + 32*integral(c4) + 7*integral(c5));
    c1 = c1 + (h*4);
    c2 = c2 + (h*4);
    c3 = c3 + (h*4);
    c4 = c4 + (h*4);
    c5 = c5 + (h*4);
    i = (i+4);
  }
return booleErr = (fabsl(M_PI - est));
}

int main() {  	
long double x;
long double y;
long double trapErr;
long double simps13Err;
long double simps38Err;
long double boolesErr;
int num = 12; //Value for num given from the PDF (12, 786433)

x=0;
y=1;

  while(num < 786433) { //Value for num given from the PDF (12, 786433)
   trapErr = trapezoid(x,y,num);
   simps13Err = simps13(x,y,num); 
   simps38Err = simps38(x,y,num);
   boolesErr = booles(x,y,num);
    printf("Num: %d | Trapezoid: %0.10Le | (1/3)rd: %0.10Le | (3/8)ths: %0.10Le | Booles: %0.10Le\n",num,trapErr,simps13Err,simps38Err,boolesErr);
    num = 2*num;
  }
  return 0;
}

