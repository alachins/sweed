////////////////////////////////////////////////////////////////////////////////
// File: exponential_integral_Ei.c                                            //
// Routine(s):                                                                //
//    Exponential_Integral_Ei                                                 //
//    xExponential_Integral_Ei                                                //
////////////////////////////////////////////////////////////////////////////////

/* #include <math.h>           // required for fabsl(), expl() and logl()         */
/* #include <float.h>          // required for LDBL_EPSILON, DBL_MAX */
/* #include <mpfr.h> */
/* #include <stdlib.h> */
/* #include <stdio.h> */
/* #include <string.h> */


/* #define ACC 1000 */

#include "SweeD.h"

//                         Internally Defined Routines                        //
double      Exponential_Integral_Ei( double x );
long double xExponential_Integral_Ei( long double x );

static long double Continued_Fraction_Ei( long double x );
static long double Power_Series_Ei( long double x );
static long double Argument_Addition_Series_Ei( long double x);


static void mpfr_Continued_Fraction_Ei( mpfr_t ei, mpfr_t x);
static void mpfr_Power_Series_Ei(mpfr_t ei, mpfr_t x);


//                         Internally Defined Constants                       //
static const long double epsilon = 10.0 * LDBL_EPSILON;


////////////////////////////////////////////////////////////////////////////////
// double Exponential_Integral_Ei( double x )                                 //
//                                                                            //
//  Description:                                                              //
//     The exponential integral Ei(x) is the integral with integrand          //
//                             exp(t) / t                                     //
//     where the integral extends from -inf to x.                             //
//     Note that there is a singularity at t = 0.  Therefore for x > 0, the   //
//     integral is defined to be the Cauchy principal value:                  //
//          lim { I[-inf, -eta] exp(-t) dt / t + I[eta, x] exp(-t) dt / t }   //
//     in which the limit is taken as eta > 0 approaches 0 and I[a,b]         //
//     denotes the integral from a to b.                                      //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
//     If x = 0.0, then Ei is -inf and -DBL_MAX is returned.                  //
//                                                                            //
//  Example:                                                                  //
//     double y, x;                                                           //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Exponential_Integral_Ei( x );                                      //
////////////////////////////////////////////////////////////////////////////////
double Exponential_Integral_Ei( double x )
{
   return (double) xExponential_Integral_Ei( (long double) x);
}


////////////////////////////////////////////////////////////////////////////////
// long double xExponential_Integral_Ei( long double x )                      //
//                                                                            //
//  Description:                                                              //
//     The exponential integral Ei(x) is the integral with integrand          //
//                             exp(t) / t                                     //
//     where the integral extends from -inf to x.                             //
//     Note that there is a singularity at t = 0.  Therefore for x > 0, the   //
//     integral is defined to be the Cauchy principal value:                  //
//          lim { I[-inf, -eta] exp(-t) dt / t + I[eta, x] exp(-t) dt / t }   //
//     in which the limit is taken as eta > 0 approaches 0 and I[a,b]         //
//     denotes the integral from a to b.                                      //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of the exponential integral Ei().         //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
//     If x = 0.0, then Ei is -inf and -DBL_MAX is returned.                  //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = xExponential_Integral_Ei( x );                                     //
////////////////////////////////////////////////////////////////////////////////

long double xExponential_Integral_Ei( long double x )
{
   if ( x < -5.0L ) return Continued_Fraction_Ei(x);
   if ( x == 0.0L ) return -DBL_MAX;
   if ( x < 6.8L )  return Power_Series_Ei(x);
   if ( x < 50.0L ) return Argument_Addition_Series_Ei(x);
   return Continued_Fraction_Ei(x);
}



void mpfr_Exponential_Integral_Ei( mpfr_t ei, mpfr_t x)
{

  if( mpfr_cmp_d(x, -5.0) < 0)
    mpfr_Continued_Fraction_Ei(ei, x );
  else if( mpfr_cmp_d(x, 6.8) < 0)
    mpfr_Power_Series_Ei(ei, x);
  else if( mpfr_cmp_d(x, 50.0) < 0)
    mpfr_eint(ei, x, GMP_RNDU);
  else if (mpfr_cmp_d(x, 50.0) >= 0)
    mpfr_Continued_Fraction_Ei(ei, x);
  else
    assert(0);
}



static void mpfr_Continued_Fraction_Ei( mpfr_t t, mpfr_t x )
{
  
  mpfr_t epsilon1;
  mpfr_t Am1, A0, Bm1, B0, a, b, Ap1, Bp1, j;
  mpfr_t tmp, tmp1, tmp2, tmp3, tmp4;
  mpfr_inits2(ACC, epsilon1, Am1, A0, Bm1, B0, a, b, Ap1, Bp1, tmp, tmp1, tmp2, tmp3, tmp4, j, (mpfr_ptr) 0);

  mpfr_set_ld(epsilon1, LDBL_EPSILON, GMP_RNDU);

  
  /* mpfr_set_ld(epsilon1, 1.0L, GMP_RNDU); */
  /* mpfr_mul_d(epsilon1, epsilon1, 0.5, GMP_RNDU); */
  /* mpfr_exp_t e = mpfr_get_emin(); */
  /* mpfr_set_exp(epsilon1, e); */
  

  /* mpfr_out_str(NULL, 10, 20, epsilon1, GMP_RNDU); */
  /* putchar('\n'); */
  
  mpfr_set_ld(Am1, 1.0L, GMP_RNDU);
  mpfr_set_ld(A0, 0.0L, GMP_RNDU);
  mpfr_set_ld(Bm1, 0.0L, GMP_RNDU);
  mpfr_set_ld(B0, 1.0L, GMP_RNDU);
  mpfr_set_ld(a, 0.0L, GMP_RNDU); //!
  mpfr_set_ld(b, 0.0L, GMP_RNDU);
  mpfr_set_ld(Ap1, 0.0L, GMP_RNDU);
  mpfr_set_ld(Bp1, 0.0L, GMP_RNDU);

  /* printf("*x: "); */
  /* mpfr_out_str(NULL, 10, 10, x, GMP_RNDU); */

  
  mpfr_exp(a, x, GMP_RNDU);
  
  mpfr_set(b, x, GMP_RNDU);
  mpfr_mul_si(b, b, -1, GMP_RNDU);
  mpfr_add_d(b, b, 1.0, GMP_RNDU);

  mpfr_set(Ap1, b, GMP_RNDU);
  mpfr_mul(Ap1, Ap1, A0, GMP_RNDU);
  mpfr_set(tmp, a, GMP_RNDU);
  mpfr_mul(tmp, tmp, Am1, GMP_RNDU);
  mpfr_add(Ap1, Ap1, tmp, GMP_RNDU);

  mpfr_set(Bp1, b, GMP_RNDU);
  mpfr_mul(Bp1, Bp1, B0, GMP_RNDU);
  mpfr_set(tmp, a, GMP_RNDU);
  mpfr_mul(tmp, tmp, Bm1, GMP_RNDU);
  mpfr_add(Bp1, Bp1, tmp, GMP_RNDU);

    
  mpfr_set_ld(j, 1.0L, GMP_RNDU);
    
  mpfr_set_ld(a, 1.0L, GMP_RNDU);
  /* a = 1.0L; */
  
  
  mpfr_set(tmp, Ap1, GMP_RNDU);
  mpfr_mul(tmp, tmp, B0, GMP_RNDU);
  mpfr_set(tmp1, A0, GMP_RNDU);
  mpfr_mul(tmp1, tmp1, Bp1, GMP_RNDU);
  mpfr_sub(tmp, tmp, tmp1, GMP_RNDU);
  mpfr_abs(tmp, tmp, GMP_RNDU);
  
  mpfr_abs(tmp1, tmp1, GMP_RNDU);
  mpfr_mul(tmp1, tmp1, epsilon1, GMP_RNDU);
  
  while (  mpfr_cmp(tmp, tmp1) > 0)
    {
  
      /* printf("** "); */
      /* mpfr_out_str(NULL, 10, 10, tmp, GMP_RNDU); */
      /* printf("\t"); */
      /* mpfr_out_str(NULL, 10, 10, tmp1, GMP_RNDU); */
      /* printf("\n"); */
      
      mpfr_set(tmp2, Bp1, GMP_RNDU);
      mpfr_abs(tmp2, tmp2, GMP_RNDU);
      
      if( mpfr_cmp_ld(tmp2, 1.0L) > 0){
	
	mpfr_set(tmp3, A0, GMP_RNDU);
	mpfr_div(tmp3, tmp3, Bp1, GMP_RNDU);
	mpfr_set(Am1, tmp3, GMP_RNDU);
	
	mpfr_set(tmp3, Ap1, GMP_RNDU);
	mpfr_div(tmp3, tmp3, Bp1, GMP_RNDU);
	mpfr_set(A0, tmp3, GMP_RNDU);
	
	mpfr_set(tmp3, B0, GMP_RNDU);
	mpfr_div(tmp3, tmp3, Bp1, GMP_RNDU);
	mpfr_set(Bm1, tmp3, GMP_RNDU);
	
	mpfr_set_ld(B0, 1.0L, GMP_RNDU);
      } 
      else {
	mpfr_set(Am1, A0, GMP_RNDU);
	mpfr_set(A0, Ap1, GMP_RNDU);
	mpfr_set(Bm1, B0, GMP_RNDU);
	mpfr_set(B0, Bp1, GMP_RNDU);
      }
      
      mpfr_set(tmp3, j, GMP_RNDU);
      mpfr_mul(tmp3, tmp3, tmp3, GMP_RNDU);
      mpfr_mul_si(tmp3, tmp3, -1, GMP_RNDU);
      mpfr_set(a, tmp3, GMP_RNDU);
      
      mpfr_add_d(b, b, 2.0, GMP_RNDU);
      
      mpfr_set(tmp3, b, GMP_RNDU);
      mpfr_mul(tmp3, tmp3, A0, GMP_RNDU);
      mpfr_set(Ap1, tmp3, GMP_RNDU);
      mpfr_set(tmp3, a, GMP_RNDU);
      mpfr_mul(tmp3, tmp3, Am1, GMP_RNDU);
      mpfr_add(Ap1, Ap1, tmp3, GMP_RNDU);
      
      mpfr_set(tmp3, b, GMP_RNDU);
      mpfr_mul(tmp3, tmp3, B0, GMP_RNDU);
      mpfr_set(Bp1, tmp3, GMP_RNDU);
      mpfr_set(tmp3, a, GMP_RNDU);
      mpfr_mul(tmp3, tmp3, Bm1, GMP_RNDU);
      mpfr_add(Bp1, Bp1, tmp3, GMP_RNDU);
      
      mpfr_add_d(j, j, 1.0, GMP_RNDU);

      
      mpfr_set(tmp, Ap1, GMP_RNDU);
      mpfr_mul(tmp, tmp, B0, GMP_RNDU);
      mpfr_set(tmp1, A0, GMP_RNDU);
      mpfr_mul(tmp1, tmp1, Bp1, GMP_RNDU);
      mpfr_sub(tmp, tmp, tmp1, GMP_RNDU);
      mpfr_abs(tmp, tmp, GMP_RNDU);
      
      mpfr_abs(tmp1, tmp1, GMP_RNDU);
      mpfr_mul(tmp1, tmp1, epsilon1, GMP_RNDU);
      
    }
    
  mpfr_set(t, Ap1, GMP_RNDU);
  mpfr_mul_si(t, t, -1, GMP_RNDU);
  mpfr_div(t, t, Bp1, GMP_RNDU);
  
  mpfr_clears(epsilon1, Am1, A0, Bm1, B0, a, b, Ap1, Bp1, tmp, tmp1, tmp2, tmp3, tmp4, j, (mpfr_ptr) 0);
  
}





////////////////////////////////////////////////////////////////////////////////
// static long double Continued_Fraction_Ei( long double x )                  //
//                                                                            //
//  Description:                                                              //
//     For x < -5 or x > 50, the continued fraction representation of Ei      //
//     converges fairly rapidly.                                              //
//                                                                            //
//     The continued fraction expansion of Ei(x) is:                          //
//        Ei(x) = -exp(x) { 1/(-x+1-) 1/(-x+3-) 4/(-x+5-) 9/(-x+7-) ... }.    //
//                                                                            //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
////////////////////////////////////////////////////////////////////////////////

static long double Continued_Fraction_Ei( long double x )
{
   long double Am1 = 1.0L;
   long double A0 = 0.0L;
   long double Bm1 = 0.0L;
   long double B0 = 1.0L;
   long double a = expl(x);
   long double b = -x + 1.0L;
   long double Ap1 = b * A0 + a * Am1;
   long double Bp1 = b * B0 + a * Bm1;
   int j = 1;

   a = 1.0L;
   while ( fabsl(Ap1 * B0 - A0 * Bp1) > epsilon * fabsl(A0 * Bp1) ) {

     printf("ORIG %Le\n", fabsl(Ap1 * B0 - A0 * Bp1));

      if ( fabsl(Bp1) > 1.0L) {
         Am1 = A0 / Bp1;
         A0 = Ap1 / Bp1;
         Bm1 = B0 / Bp1;
         B0 = 1.0L;
      } else {
         Am1 = A0;
         A0 = Ap1;
         Bm1 = B0;
         B0 = Bp1;
      }
      a = -j * j;
      b += 2.0L;
      Ap1 = b * A0 + a * Am1;
      Bp1 = b * B0 + a * Bm1;
      j += 1;
   }
   return (-Ap1 / Bp1);
}


////////////////////////////////////////////////////////////////////////////////
// static long double Power_Series_Ei( long double x )                        //
//                                                                            //
//  Description:                                                              //
//     For -5 < x < 6.8, the power series representation for                  //
//     (Ei(x) - gamma - ln|x|)/exp(x) is used, where gamma is Euler's gamma   //
//     constant.                                                              //
//     Note that for x = 0.0, Ei is -inf.  In which case -DBL_MAX is          //
//     returned.                                                              //
//                                                                            //
//     The power series expansion of (Ei(x) - gamma - ln|x|) / exp(x) is      //
//        - Sum(1 + 1/2 + ... + 1/j) (-x)^j / j!, where the Sum extends       //
//        from j = 1 to inf.                                                  //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
////////////////////////////////////////////////////////////////////////////////

static long double Power_Series_Ei( long double x )
{ 
   long double xn = -x;
   long double Sn = -x;
   long double Sm1 = 0.0L;
   long double hsum = 1.0L;
   long double g = 0.5772156649015328606065121L;
   long double y = 1.0L;
   long double factorial = 1.0L;
  
   if ( x == 0.0L ) return (long double) -DBL_MAX;
 
   while ( fabsl(Sn - Sm1) > epsilon * fabsl(Sm1) ) {
      Sm1 = Sn;
      y += 1.0L;
      xn *= (-x);
      factorial *= y;
      hsum += (1.0 / y);
      Sn += hsum * xn / factorial;
   }
   return (g + logl(fabsl(x)) - expl(x) * Sn);
}


static void mpfr_Power_Series_Ei(mpfr_t t, mpfr_t x)
{ 
  mpfr_t epsilon1;
  mpfr_t xn, Sn, Sm1, hsum, g, y, factorial;
  mpfr_t tmp, tmp1, tmp2, tmp3, tmp4;

  mpfr_inits2(ACC, epsilon1, xn, Sn, Sm1, hsum, g, y, factorial, tmp, tmp1, tmp2, tmp3, tmp4, (mpfr_ptr) 0);
  
  mpfr_set_ld(epsilon1, LDBL_EPSILON, GMP_RNDU);

  /* long double xn = -x; */
  mpfr_set(xn, x, GMP_RNDU);
  mpfr_mul_si(xn, xn, -1, GMP_RNDU);
  
   /* long double Sn = -x; */
  mpfr_set(Sn, xn, GMP_RNDU);

   /* long double Sm1 = 0.0L; */
  mpfr_set_ld(Sm1, 0.0L, GMP_RNDU);

  /* long double hsum = 1.0L; */
  mpfr_set_ld(hsum, 1.0L, GMP_RNDU);

  /*long double g = 0.5772156649015328606065121L; */
  //mpfr_set_ld(g, 0.5772156649015328606065121L, GMP_RNDU);
  mpfr_const_euler(g, GMP_RNDU);
  
  /*long double y = 1.0L;*/
  mpfr_set_ld(y, 1.0L, GMP_RNDU);

  /* long double factorial = 1.0L; */
  mpfr_set_ui(factorial, 1, GMP_RNDU);
  
    
  /* if ( x == 0.0L ) return (long double) -DBL_MAX; */
  assert( mpfr_cmp_ld(x, 0.0L) != 0);
    
  mpfr_set(tmp, Sn, GMP_RNDU);
  mpfr_sub(tmp, tmp, Sm1, GMP_RNDU);
  mpfr_abs(tmp, tmp, GMP_RNDU);

  mpfr_set(tmp1, Sm1, GMP_RNDU);
  mpfr_abs(tmp1, tmp1, GMP_RNDU);
  mpfr_mul(tmp1, tmp1, epsilon1, GMP_RNDU);
  
  while ( (mpfr_cmp(tmp, tmp1) > 0) ) {
    /*Sm1 = Sn;*/
    mpfr_set(Sm1, Sn, GMP_RNDU);
    /*      y += 1.0L;*/
    mpfr_add_ui(y, y, 1, GMP_RNDU);
    /*  xn *= (-x);  */
    mpfr_set(tmp, x, GMP_RNDU);
    mpfr_mul_si(tmp, tmp, -1, GMP_RNDU);
    mpfr_mul(xn, xn, tmp, GMP_RNDU);
    /*  factorial *= y; */
    mpfr_mul(factorial, factorial, y, GMP_RNDU);
    /*  hsum += (1.0 / y); */
    mpfr_set_ld(tmp, 1.0L, GMP_RNDU);
    mpfr_div(tmp, tmp, y, GMP_RNDU);
    mpfr_add(hsum, hsum, tmp, GMP_RNDU);
    /*Sn += hsum * xn / factorial;*/
    mpfr_set(tmp, hsum, GMP_RNDU);
    mpfr_mul(tmp, tmp, xn, GMP_RNDU);
    mpfr_div(tmp, tmp, factorial, GMP_RNDU);
    mpfr_add(Sn, Sn, tmp, GMP_RNDU);

    /* update the values */
    mpfr_set(tmp, Sn, GMP_RNDU);
    mpfr_sub(tmp, tmp, Sm1, GMP_RNDU);
    mpfr_abs(tmp, tmp, GMP_RNDU);
    
    mpfr_set(tmp1, Sm1, GMP_RNDU);
    mpfr_abs(tmp1, tmp1, GMP_RNDU);
    mpfr_mul(tmp1, tmp1, epsilon1, GMP_RNDU);
  
   }
  
  /*return (g + logl(fabsl(x)) - expl(x) * Sn);*/
  mpfr_set(tmp, x, GMP_RNDU);
  mpfr_abs(tmp, tmp, GMP_RNDU);
  mpfr_log(tmp, tmp, GMP_RNDU);
  mpfr_add(tmp, tmp, g, GMP_RNDU);
  
  mpfr_set(t, tmp, GMP_RNDU);
    
  mpfr_exp(tmp, x, GMP_RNDU);
  mpfr_mul(tmp, tmp, Sn, GMP_RNDU);

  mpfr_sub(t, t, tmp, GMP_RNDU);

  mpfr_clears(epsilon1, xn, Sn, Sm1, hsum, g, y, factorial, tmp, tmp1, tmp2, tmp3, tmp4, (mpfr_ptr) 0);

}




////////////////////////////////////////////////////////////////////////////////
// static long double Argument_Addition_Series_Ei(long double x)              //
//                                                                            //
//  Description:                                                              //
//     For 6.8 < x < 50.0, the argument addition series is used to calculate  //
//     Ei.                                                                    //
//                                                                            //
//     The argument addition series for Ei(x) is:                             //
//     Ei(x+dx) = Ei(x) + exp(x) Sum j! [exp(j) expj(-dx) - 1] / x^(j+1),     //
//     where the Sum extends from j = 0 to inf, |x| > |dx| and expj(y) is     //
//     the exponential polynomial expj(y) = Sum y^k / k!, the Sum extending   //
//     from k = 0 to k = j.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
////////////////////////////////////////////////////////////////////////////////
static long double Argument_Addition_Series_Ei(long double x)
{
   static long double ei[] = {
      1.915047433355013959531e2L,  4.403798995348382689974e2L,
      1.037878290717089587658e3L,  2.492228976241877759138e3L,
      6.071406374098611507965e3L,  1.495953266639752885229e4L,
      3.719768849068903560439e4L,  9.319251363396537129882e4L,
      2.349558524907683035782e5L,  5.955609986708370018502e5L,
      1.516637894042516884433e6L,  3.877904330597443502996e6L,
      9.950907251046844760026e6L,  2.561565266405658882048e7L,
      6.612718635548492136250e7L,  1.711446713003636684975e8L,
      4.439663698302712208698e8L,  1.154115391849182948287e9L,
      3.005950906525548689841e9L,  7.842940991898186370453e9L,
      2.049649711988081236484e10L, 5.364511859231469415605e10L,
      1.405991957584069047340e11L, 3.689732094072741970640e11L,
      9.694555759683939661662e11L, 2.550043566357786926147e12L,
      6.714640184076497558707e12L, 1.769803724411626854310e13L,
      4.669055014466159544500e13L, 1.232852079912097685431e14L,
      3.257988998672263996790e14L, 8.616388199965786544948e14L,
      2.280446200301902595341e15L, 6.039718263611241578359e15L,
      1.600664914324504111070e16L, 4.244796092136850759368e16L,
      1.126348290166966760275e17L, 2.990444718632336675058e17L,
      7.943916035704453771510e17L, 2.111342388647824195000e18L,
      5.614329680810343111535e18L, 1.493630213112993142255e19L,
      3.975442747903744836007e19L, 1.058563689713169096306e20L
   };
   int  k = (int) (x + 0.5);
   int  j = 0;
   long double xx = (long double) k;
   long double dx = x - xx;
   long double xxj = xx;
   long double edx = expl(dx);
   long double Sm = 1.0L;
   long double Sn = (edx - 1.0L) / xxj;
   long double term = DBL_MAX;
   long double factorial = 1.0L;
   long double dxj = 1.0L;

   while (fabsl(term) > epsilon * fabsl(Sn) ) {
      j++;
      factorial *= (long double) j;
      xxj *= xx;
      dxj *= (-dx);
      Sm += (dxj / factorial);
      term = ( factorial * (edx * Sm - 1.0L) ) / xxj;
      Sn += term;
   }
   
   return ei[k-7] + Sn * expl(xx); 
}
