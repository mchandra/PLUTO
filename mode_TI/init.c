/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sepy 10, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rdt dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
//  double delrho;
//  long iseed;
//  int i;
//  iseed=1;
//  v[RHO] = g_inputParam[n_e]*1.1718*( 1. + g_inputParam[amp] \
           *sin(10.0*2.0*CONST_PI*x1/4.e5)*sin(10.0*2.0*CONST_PI*x2/4.e5));
  v[RHO] = 0.0;
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  #if HAVE_ENERGY
//    v[PRS] =  g_inputParam[n_e]*1.1718*CONST_kB*0.78*1.16e7/(0.6167*CONST_mp)/1.0e10;
  v[PRS] = 0.0;
  #endif
  v[TRC] = 0.0;
 
  #if PHYSICS == MHD || PHYSICS == RMHD
//   v[BX1] = g_inputParam[B_in]*1.0/sqrt(2.0)/sqrt(4.0*CONST_PI*CONST_mp*1.0e10);
//   v[BX2] = g_inputParam[B_in]*1.0/sqrt(2.0)/sqrt(4.0*CONST_PI*CONST_mp*1.0e10);
   v[BX1] = 0.0;
   v[BX2] = 0.0;
   v[BX3] = 0.0;
  #endif
}

void Init_ICM(Data *d, Grid *grid)
/* Initializes the density, velocity and pressure by calculating the pressure
 * using Newton-Raphson. 
 * Written by Deovrat Prasad on 6th April, 2015.
 ************************************************************************** */
{
     int i, j, k;
     double *x1, *x2, *x3;
     x1 = grid[IDIR].x;
     x2 = grid[JDIR].x;
     x3 = grid[KDIR].x;
     TOT_LOOP(k,j,i){
     d->Vc[RHO][k][j][i] =  g_inputParam[n_e]*1.1718;
     D_EXPAND(
     d->Vc[VX1][k][j][i]=0.0;  ,
     d->Vc[VX2][k][j][i]=0.0;  ,
     d->Vc[VX3][k][j][i]=0.0;
     )
     #if HAVE_ENERGY
     d->Vc[PRS][k][j][i] = g_inputParam[n_e]*1.1718*CONST_kB*g_inputParam[TkeV]*1.16e7/(0.6167*CONST_mp)/1.0e10;
     #endif
     }
/* Introduction of Large Amplitude Random Pertubation in 
 * density at small scales.
 */
     int n,n1, l, l1;
     long iseed;
     double kx, ky;
     double phi, aklm, xs, ys, pert;
     double delrho[NX3_TOT][NX2_TOT][NX1_TOT];
     iseed = 1;
     memset(delrho, 0, sizeof delrho);
     for(n = 4; n<=20; n++){
     for(n1 = -n; n1<=n; n1+=2*n){
        kx = 2.0*CONST_PI*n1/(2.*g_inputParam[RMAX_1]);
        for(l = 4; l<=20; l++){
        for(l1 = -l; l1<=l; l1+=2*l){
           ky = 2.0*CONST_PI*l1/(2.*g_inputParam[RMAX_1]);
           pert=ran2(&iseed);
           phi = 2.0*CONST_PI*pert;
           pert=ran2(&iseed);
           aklm = 0.15*g_inputParam[amp]*(0.5-pert)/sqrt(pow(n,2)+pow(l,2));

           TOT_LOOP(k,j,i){
              delrho[k][j][i] += aklm*d->Vc[RHO][k][j][i]*cos(phi + kx*x1[i] + ky*x2[j] );
           }
        }}
    }}
    TOT_LOOP(k,j,i){
       d->Vc[RHO][k][j][i] += delrho[k][j][i];
       #if PHYSICS == MHD 
       D_EXPAND(
       d->Vc[BX1][k][j][i] = g_inputParam[B_in]*1.0/sqrt(2.0)/sqrt(4.0*CONST_PI*CONST_mp*1.0e10);,
       d->Vc[BX2][k][j][i] = g_inputParam[B_in]*1.0/sqrt(2.0)/sqrt(4.0*CONST_PI*CONST_mp*1.0e10);,
       d->Vc[BX3][k][j][i] = 0.0;
       )
       #ifdef STAGGERED_MHD
       D_EXPAND(
        d->Vs[BX1s][k][j][i] = d->Vc[BX1][k][j][i]; ,
        d->Vs[BX2s][k][j][i] = d->Vc[BX2][k][j][i]; ,
        d->Vs[BX3s][k][j][i] = d->Vc[BX3][k][j][i]; 
       )
       #endif
       #endif
    }
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{
  int k, j, i;
  double dvol=0.0, dmass=0.0, mass=0.0, vol=0.0, momx1=0.0, momx2=0.0, KE1=0.0, KE2=0.0, TE=0.0;
  double *dx1, *dx2, *dx3;
  double sendarray[48],recvarray[48];
  FILE *fdiag;

  dx1 = grid[IDIR].dx; dx2 = grid[JDIR].dx; dx3 = grid[KDIR].dx;

  DOM_LOOP(k,j,i){
    dvol = dx1[i]*dx2[j];
//*dx3[k]; 
    dmass = dvol*d->Vc[RHO][k][j][i];
    vol += dvol;
//    print("%20.10e %20.10e\n", dvol, dmass);
    mass += dmass;
    momx1 += d->Vc[VX1][k][j][i]*dmass;
    momx2 += d->Vc[VX2][k][j][i]*dmass;
    TE += d->Vc[PRS][k][j][i]*dvol/(g_gamma-1.0);
    KE1 += 0.5*d->Vc[VX1][k][j][i]*d->Vc[VX1][k][j][i]*dmass;
    KE2 += 0.5*d->Vc[VX2][k][j][i]*d->Vc[VX2][k][j][i]*dmass;
  }
    if (g_stepNumber==0){
      fdiag = fopen ("Diagnostics.out", "w");
      fprintf(fdiag,"# [1]time [2]g_dt [3]mass [4]TE [5]KE1 [6]KE2 [7]KE3 [8]MOM1 [9]MOM2 [10]MOM3\n");
    }
    else fdiag = fopen ("Diagnostics.out", "a");

  #ifdef PARALLEL
   sendarray[0]=mass; sendarray[1]=vol; sendarray[2]=momx1; sendarray[3]=momx2;
   sendarray[4]=KE1; sendarray[5]=KE2; sendarray[6]=TE;

   MPI_Reduce (sendarray, recvarray, 7, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);

   mass=recvarray[0]; vol=recvarray[1]; momx1=recvarray[2]; momx2=recvarray[3];
   KE1=recvarray[4]; KE2=recvarray[5]; TE=recvarray[6];

   if (prank==0){
  #endif
      fprintf(fdiag,"%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",g_time,
      g_dt, mass, TE, KE1, KE2, 0.0, momx1, momx2, 0.0);
  fclose(fdiag);
  #ifdef PARALLEL
   }
  #endif
}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;

  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;

  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){};
  }

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X1_END){  /* -- X1_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_END){  /* -- X2_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_END){  /* -- X3_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif

 #define IM1 2147483563
 #define IM2 2147483399
 #define AM (1.0/IM1)
 #define IMM1 (IM1-1)
 #define IA1 40014
 #define IA2 40692
 #define IQ1 53668
 #define IQ2 52774
 #define IR1 12211
 #define IR2 3791
 #define NTAB 32
 #define NDIV (1+IMM1/NTAB)
 #define EPS 1.2e-7
 #define RNMX (1.0-EPS)

float ran2(long *idum)

/*!
 *  Random number generator taken from Numerical Recipes in C.
 *  Returns a uniform deviate between 0.0 and 1.0.
 *  Call with idum a negative integer to initialize; 
 *  thereafter do not alter idum between successive deviates in a sequence.
 *  RNMX should approximate the largest floating value that is less than 1.
 ************************************************************************* */
{
    int j ;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0){
       if (-(*idum) < 1) *idum =1;
       else *idum = -(*idum);
       idum2 = (*idum);
       for (j=NTAB + 7; j>=0; j--) {
           k=(*idum)/IQ1;
           *idum=IA1*(*idum-k*IQ1) - k*IR1;
           if (*idum < 0) *idum += IM1;
           if (j < NTAB) iv[j] = *idum;
       }
       iy = iv[0];
    }
    k=(*idum)/IQ1;
    *idum =IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 <0) idum2 += IM2;
    j = iy/NDIV;
    iy=iv[j]-idum2;
    iv[j]= *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

