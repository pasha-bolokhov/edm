#include <stdio.h>
#include <math.h>


double rho = - 5.535;   /* radial integral divided by "a" */

double Z   = 55.0;      /* for Cesium */
double alpha = 1./137.0;
double Ryd   = 1.097e+05;       /* cm^{-1}, Rydberg constant */

/*
 * Energy levels for Cesium
 * in cm^{-1}
 * "Atomic Energy Levels", C.E.Moore, v.3
 */
double  Cs_E_s12  = 0.0;        /* 6s1/2, the ground state */
double  Cs_E_p12  = 11178.24;   /* 6p1/2                   */
double  Cs_E_p32  = 11732.35;   /* 6p3/2                   */

double  eV_cm1    = 8.0655e+03; /* Conversion coef eV -> cm^{-1} */

double  Cs_IP     = 3.893;      /* Ionization potential, Volts */
double  Cs_IE     = 0.0;        /* Ionization energy    */


/*
 * Effective principal numbers
 */
double  nu_s12      = 0.0; 
double  nu_p12      = 0.0;
double  nu_p32      = 0.0;


double  ggamma   = 0.0;          /* s1/2 */
double  ggamma2  = 0.0;          /* p1/2, for a cross check */
double  ggamma1  = 0.0;          /* p3/2 */

double  K_LI     = 0.0;          /* Lorentz-invariant enhancement     */
double  K_appr   = 0.0;          /* LV enhancement, up to (Z\alpha)^2 */
double  K        = 0.0;          /* exact LV enhancement              */

double  theta    = 0.0;
double  eta      = 0.0;
double  zeta     = 0.0;
double  zeta1    = 0.0;
double  kappa    = 0.0;
double  kappa1   = 0.0;

double  tgamma1  = 0.0;
double  tgamma2  = 0.0;


int 
main(int argc, char *argv[])
{
  printf("------------------------------------------------------------\n");
  printf(" Calculates EDM enhancement due to the CPT-odd LV operator.\n");

  /* Calculate the Effective Principal Numbers first */
  Cs_IE = Cs_IP * eV_cm1;              /* Ionization energy */
  nu_s12  = sqrt( Ryd / ( Cs_IE   -  Cs_E_s12 ) );
  nu_p12  = sqrt( Ryd / ( Cs_IE   -  Cs_E_p12 ) );
  nu_p32  = sqrt( Ryd / ( Cs_IE   -  Cs_E_p32 ) );
  printf("------------------------------------------------------------\n");
  printf(" Effective Principal Numbers for Cs\n");
  printf(" ionization energy = %g cm^{-1}\n", Cs_IE);
  printf(" nu_s12 = %g\n", nu_s12);
  printf(" nu_p12 = %g\n", nu_p12);
  printf(" nu_p32 = %g\n", nu_p32);

  /* Calculate gamma's */
  ggamma  = sqrt ( pow(( 0.5 + 0.5 ), 2.0)  -  pow( Z*alpha, 2.0 ) );
  ggamma2 = sqrt ( pow(( 0.5 + 0.5 ), 2.0)  -  pow( Z*alpha, 2.0 ) );
  ggamma1 = sqrt ( pow(( 1.5 + 0.5 ), 2.0)  -  pow( Z*alpha, 2.0 ) );
  printf("------------------------------------------------------------\n");
  printf(" gamma(s1/2) = %g, gamma(p1/2) = %g, gamma(p3/2) = %g\n",
	 ggamma, ggamma2, ggamma1);
  
  K_LI = - (16./3.) *
    pow(Z,3.0) * pow(alpha,2.0) * rho * 
    (Ryd / ( Cs_E_p12 - Cs_E_s12 )) /
      ( ggamma * ( 4.0 * pow(ggamma,2.0)  -  1.0 ) 
        * pow( nu_s12 * nu_p12, 1.5 ) );
  printf(" Lorentz-invariant enhancement factor = %g\n", K_LI);

  /* Lorentz-Noninvariant factor only up to (Z\alpha)^2 */
  K_appr = 
  (7./30.) *
    pow(Z,3.0) * pow(alpha, 2.0) * 2.0 * Ryd * rho /
    ( pow ( nu_s12 * nu_p32, 1.5 ) * ( - ( Cs_E_p32 - Cs_E_s12 )) );
  printf("------------------------------------------------------------\n");
  printf(" LV enhancement up to leading (Z\\alpha)^2 = %g\n", K_appr);

  /* Exact Lorentz-violating factor */
  printf("------------------------------------------------------------\n");
  theta = ggamma1 + ggamma;
  eta   = ggamma1 - ggamma;
  kappa  = -1.0;
  kappa1 = -2.0;
  zeta  =  ggamma  + kappa;
  zeta1  = ggamma1  + kappa1;

  tgamma1 = 2.082532;    /* Gamma( 2. + eta ) */
  tgamma2 = 1.027;       /* Gamma( 2. - eta ) */
  
  K  =

    2.0   *        /* overall factor, from PT */
    ( - 2. / 3. ) *
    Z * 2.0 * Ryd *         /*  alpha/a  = 2 * Ryd   */
    rho /                   /*  R1/a     = rho       */
    (  pow ( nu_s12 * nu_p32, 1.5 ) * ( - ( Cs_E_p32 - Cs_E_s12 ))  );

  K  = K /
    (  ( theta - 1.0 )  * 
/*
tgamma( 2. + eta ) * tgamma( 2. - eta ) 
*/
       tgamma1 * tgamma2 );

  K  = K *
    (  
     ( 1 - eta * eta )   
        - 
     ( zeta1 * ( 1.0 - eta )   +   zeta * ( 1.0 + eta ) ) / theta 
        +
     2.0 * 
     (  zeta * zeta1  -  (1./5.) * pow( Z* alpha, 2.0 ) ) /
                         ( theta * ( theta + 1.0 ) )
     );
  printf( " 2 + eta = %g,  2 - eta = %g\n", 2.0 + eta, 2.0 - eta);
  printf("------------------------------------------------------------\n");
  printf(" exact LV enhancement factor for Cs 6p3/2 - 6s1/2 mixing = %g\n",
	 K);

  {
    double a = tgamma(1.0);
    double b = tgamma(2.0);
    double c = exp(gamma(3.0));
    double d = exp(gamma(4.0));
    
    printf("%e %e %e %e\n", a,b,c,d);
    printf(" 2! = %g, 3! = %g, 4! = %g, 5! = %g, "
	   "Gamma(1/2) = sqrt(pi) = %g, should be %g\n", 
	   tgamma(3.0), 	 tgamma(4.0), 	 tgamma(5.0), 	 tgamma(6.0), 
	   tgamma(0.5), sqrt(3.141592653589793238));
  }
  
}

