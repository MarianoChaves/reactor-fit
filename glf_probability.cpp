/* GLoBESfit -- GLoBES fitting tools
*  (C) 2019-2020 The GLoBESfit Team
*
* GLoBESfit is mainly intended for academic purposes. Proper
* credit must be given if you use GLoBESfit or parts of it. Please
* read the section 'Credit' in the README file.
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* (C) 2019, 2020 Patrick Huber, J. M. Berryman */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <gsl/gsl_sf_expint.h>
#include <globes/globes.h>   /* GLoBES library */

#include "glf_types.h"
#include "glf_probability.h"

#define GLF_LE_KM_GEV 1.0/(GLB_EV_TO_KM_FACTOR*4.0*1E9)



double fT(double E)
{
  return (3.77735 - 0.380584*E + 0.0189118*E*E - 0.000333625*E*E*E);
	//return(1./(2.12-22.81/pow(E,1./3.))+35.36/sqrt(E)-11.72/E);
}

/***************************************************************************
 *                        H E L P E R   F U N C T I O N S                  *
 ***************************************************************************/

/* Minimum of two numbers */
static inline double min(double x, double y)
{
  if (x < y)
    return x;
  else
    return y;
}

/* Square of real number */
static inline double square(double x)
{
  return x*x;
}


/**********************************
*  INPUT FOR PROBABILITY ENGINES  *
***********************************/

static double delta_atm; /* Delta m_{31}^2 */
static double theta_13;
static double theta_12;
static double delta_solar; /* Delta m_{21}^2 */
static double real;
static double imag;

#define me 0.511
#define Delta 1.29

#define gT 0.987
#define gA 1.27
#define gS 1.02

/**************************************************************************
*   This function is the result of                                        *
*       \int_La^Lb dL \sin^2(q*L)/L^2 / \int_La^Lb dL 1/L^2		  *
***************************************************************************/

/***************************************************************************
 * Store oscillation parameters in internal data structures.               *
 * For more sophisticated probability engines, this would be the right     *
 * place to pre-compute the mixing matrix and parts of the Hamiltonian in  *
 * order to speed up the calls to the actual probability matrix function.  *
 ***************************************************************************/

int glf_get_oscillation_parameters(glb_params p, void *user_data)
{
  theta_13 = glbGetOscParams(p, GLF_THETA_13);
  imag =  glbGetOscParams(p, GLF_IMAG);
  delta_atm =  glbGetOscParams(p, GLF_DELTA_ATM);
  theta_12 =  glbGetOscParams(p, GLF_THETA_12);
  delta_solar =  glbGetOscParams(p, GLF_DELTA_SOLAR);
  real =  glbGetOscParams(p, GLF_REAL);
  return 0;
}

/***************************************************************************
 * Write oscillation parameters from internal data structures into p.      *
 ***************************************************************************/

int glf_set_oscillation_parameters(glb_params p, void *user_data)
{
  glbSetOscParams(p,theta_13,GLF_THETA_13);  
  glbSetOscParams(p,imag,GLF_IMAG);
  glbSetOscParams(p,delta_atm,GLF_DELTA_ATM);
  glbSetOscParams(p,theta_12,GLF_THETA_12);
  glbSetOscParams(p,delta_solar,GLF_DELTA_SOLAR);
  glbSetOscParams(p,real,GLF_REAL);
  return 0;
}

/*********************************
*  DEFINING PROBABILITY ENGINES  *
**********************************/



int glf_nsi_probability_pre(
double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{


  int i, j, k;
  glf_distance_data *distances = (glf_distance_data *) user_data;
  size_t dlength= distances->length;
  double (*data)[2]= distances->data;
  double (*data2)[2]= distances->data2;

  double numerator = 0.0;
  double numerator2= 0.0;

  std::complex<double> et, eu;
  et={ real, imag };
  eu={ 0, 0 };

  double th12 = theta_12;
  double th13 = theta_13;

  double dm21 = delta_solar;
  double dm31 = delta_atm;

  double c12 = cos(th12);
  double s12 = sin(th12);

  double c13=cos(th13);
  double s13=sin(th13);

  double DD21 = GLF_LE_KM_GEV*dm21/E; /* Recall: E in GeV! */
  double DD31 = GLF_LE_KM_GEV*dm31/E;

  double lambda[3] = {0.0, DD21, DD31};

  //if(distances->name==12 & floor(E*100000)==590 & dm21 > 0.000088)
  //{printf("\n DM = %lf, Energy = %lf, lambda = %lf, DM/E=%lf\n", dm21*10000, E*1000, lambda[1], GLF_LE_KM_GEV*dm21/E);};

  double px, pxx, dx, dxx;
   //SCALAR
  
    dx=-gS/(1+3*gA*gA)*me/(E*1000-Delta);
    dxx=gS*gS/(1+3*gA*gA);
    px=0;
    pxx=gS*gS/(3*gA*gA);
/*  
  //TENSOR
  px=-gT/gA*me/fT(E*1000);
  pxx=gT*gT/(gA*gA);
  dx=3*gA*gT/(1+3*gA*gA)*me/(E*1000-Delta);
  dxx=3*gT*gT/(1+3*gA*gA);
*/
/* //LEFT
    px=1;
    pxx=1;
    dx=1;
    dxx=1;
*/


  double ue[3], um[3], ut[3];
  std::complex<double> v[3]; /* The elements of the "e" row of the PMNS matrix */
  ue[0] = c12*c13;
  ue[1] = c13*s12;
  ue[2] = s13;
  um[0] = -s12;
	um[1] = c12;
	um[2] = 0;
  ut[0] = -c12*s13;
	ut[1] = -s12*s13;
	ut[2] = c13;

  v[0] = eu*um[0]+et*ut[0];
  v[1] = eu*um[1]+et*ut[1];
  v[2] = eu*um[2]+et*ut[2];

  double A[3][3];
  double B[3][3];
  for(k=0; k<3; k++){
    for(int l=0; l<3; l++){
 //Full 
      std::complex<double> Ap = ue[k]*ue[l] + px*std::conj(v[k])*ue[l] + px*ue[k]*v[l] + pxx*std::conj(v[k])*v[l];
      std::complex<double> Ad = ue[k]*ue[l] + dx*v[k]*ue[l] + dx*ue[k]*std::conj(v[l]) + dxx*v[k]*std::conj(v[l]);
      A[k][l] = double( std::real(Ap*Ad) );
      B[k][l] = double( std::imag(Ap*Ad) );

/* //Linear 
      A[k][l] = ue[k]*ue[l] * ue[l]*ue[k] 
                + ue[k]*ue[l]*( px*v[k]*ue[l] + px*ue[k]*v[l])
                + ue[k]*ue[l]*( dx*v[k]*ue[l] + dx*ue[k]*v[l]);
      A[k][l] += ( px*v[k]*ue[l] + px*ue[k]*v[l])*( dx*v[k]*ue[l] + dx*ue[k]*v[l]);
*/
    }
  }

  double logq; /* NB: working with log10 of q! */
  double factor = 0.0;
  double factor2 = 0.0;

  double step = (data[1][0]-data[0][0]);
  double step2 = (data2[1][0]-data2[0][0]);
  int pos = 0;
  int pos2 = 0;

  /* Set all probabilities to zero initially */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      P[i][j] = 0.0;

  //P[0][0] = 1.0;
  for(i=0;i<3;i++)
  {
    P[0][0] += A[i][i]+B[i][i];
  }

  int test = 0;
  for (j=0; j<3; j++){
    for (k=j+1; k<3; k++){
      logq = log10(fabs(lambda[k]-lambda[j])); 
  
      if (logq <= data[0][0]){
	/* For logq too small, quadratic oscillations */
        //factor = data[0][1]*pow(10.0, 2.0*(logq-data[0][0])); // for sin2(x)
        //factor = 2*data[0][1]*pow(10.0, 2.0*(logq-data[0][0])); // for sin2(x)
        factor = 1-data[0][1]*pow(10.0, 2.0*(logq-data[0][0])); // for cos(2x)
        factor2 = 2*data2[0][1]*pow(10.0, 2.0*(logq-data2[0][0])); // for sin(2x)
        test = 1;
      }
      else if (logq >= data[dlength-1][0]){
 	/* For logq too large, oscillations average out */
        //factor = 0.5; // for sin2(x)
        factor = factor2 = 0.0; // for cos(2x) and sin(2x)
        test = 2;
      }
      else{
        pos = floor( (logq-data[0][0] )/step );
        pos2 = floor( (logq-data2[0][0])/step2 );
        numerator = (data[pos+1][1]-data[pos][1]);
        numerator2 = (data2[pos2+1][1]-data2[pos2][1]);
        factor = data[pos][1] + (logq-data[pos][0])*numerator/step;
        factor2 = data2[pos2][1] + (logq-data2[pos2][0])*numerator2/step2;
        test = 3;
      }
      
      P[0][0] += 2*A[k][j]*factor+2*B[k][j]*factor2;
      //P[0][0] += 2*A[k][j]*(1-2*factor);
      //P[0][0] += -4.0*ue[k]*ue[k]*ue[j]*ue[j] * factor;
    //if(distances->name==12 & k==1 & j==0 & floor(E*100000)==590){printf("\n DM = %f, Energy = %f, factor = %f, lambda = %f, DM/E=%f\n",dm21,E*1000, factor, lambda[k],DD21);};
    }

  }

  
  return 0;
}

int glf_nsi_probability_kamland(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{
	double Ret, Reu, Iet, Ieu;
	Reu = real;
	Ieu = imag;

	Ret = 0;
	Iet = 0;

/*	//SCALAR
	double pxL=0.0;
	double pxx=gS*gS/(3*gA*gA);
	double dxL=-gS/(1+3*gA*gA)*me/(E-Delta);
	double dxx=gS*gS/(1+3*gA*gA);
*/
	//TENSOR
	double pxL=-gT/gA*me/fT(E);
	double pxx=gT*gT/(gA*gA);
	double dxL=3*gA*gT/(1+3*gA*gA)*me/(E-Delta);
	double dxx=3*gT*gT/(1+3*gA*gA);	

/*	//LEFT
	double pxL=1;
	double pxx=1;
	double dxL=1;
	double dxx=1;
*/

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

	std::complex<double> minusX = {-1,0};
	std::complex<double> I = sqrt(minusX);
	int k,l;

  for (int i=0; i < 3; i++)
    {for (int j=0; j < 3; j++)
      {P[i][j] = 0.0;}}


	std::complex<double> Ue[3];
	std::complex<double> Ueconj[3];
	std::complex<double> Um[3];
	std::complex<double> Umconj[3];
	std::complex<double> Ut[3];
	std::complex<double> Utconj[3];

	std::complex<double> UeNSI[3];
	std::complex<double> UeNSIconj[3];

	double Dm2[3][3];
	std::complex<double> nsi[3][3];
	
	UeNSI[0]=UeNSI[1]=UeNSI[2]=UeNSIconj[0]=UeNSIconj[1]=UeNSIconj[2]=0.0;

	Dm2[1][0]=delta_solar; Dm2[0][1]=-delta_solar; // dm21
	Dm2[2][0]=delta_atm; Dm2[0][2]=-delta_atm; // dm31

	Dm2[2][1]=Dm2[2][0]-Dm2[1][0]; // dm32=dm31-dm21
	Dm2[1][2]=-Dm2[2][1];

	Dm2[0][0]=Dm2[1][1]=Dm2[2][2]=0.0;

  double t12 = theta_12;
  double t13 = theta_13;

	Ue[0] = cos(t12)*cos(t13);
	Ue[1] = sin(t12)*cos(t13);
	Ue[2] = sin(t13);
	Ueconj[0] = cos(t12)*cos(t13);
	Ueconj[1] = sin(t12)*cos(t13);
	Ueconj[2] = sin(t13);

	Um[0] = -sin(t12);
	Um[1] = cos(t12);
	Um[2] = 0;
	Umconj[0] = -sin(t12);
	Umconj[1] = cos(t12);
	Umconj[2] = 0;

	Ut[0] = -cos(t12)*sin(t13);
	Ut[1] = -sin(t12)*sin(t13);
	Ut[2] = cos(t13);
	Utconj[0] = -cos(t12)*sin(t13);
	Utconj[1] = -sin(t12)*sin(t13);
	Utconj[2] = cos(t13);

	nsi[0][0]={0.0,0.0};
	nsi[1][1]={0.0,0.0};
	nsi[2][2]={0.0,0.0};

	nsi[0][1]={Reu,Ieu};
	nsi[1][0]={Reu,-Ieu};

	nsi[0][2]={Ret,Iet};
	nsi[2][0]={Ret,-Iet};

	nsi[1][2]={0.0,0.0};
	nsi[2][1]={0.0,0.0};


	UeNSI[0]=nsi[0][1]*Um[0]+nsi[0][2]*Ut[0];
	UeNSI[1]=nsi[0][1]*Um[1]+nsi[0][2]*Ut[1];
	UeNSI[2]=nsi[0][1]*Um[2]+nsi[0][2]*Ut[2];

	UeNSIconj[0]=conj(UeNSI[0]);
	UeNSIconj[1]=conj(UeNSI[1]);
	UeNSIconj[2]=conj(UeNSI[2]);

	double baseline[20]={160,179,191,138,214,146,88,349,345,295,431,401,561,755,830,783,712,986,735,709};
	double power[20]={24.3,13.7,10.2,4.5,10.6,4.9,1.6,14.2,13.2,3.3,6.5,3.8,6.0,10.1,5.3,3.3,11.5,17.4,9.2,8.2};

	std::complex<double> Prob=0;
	double powerall=0;
	for (int item=0; item<20; item+=1)
	{
	for(k=0;k<3;k++)
	{
	for(l=0;l<3;l++)
	{	
		double L=baseline[item];
		double Po=(1/181.7)*power[item];
		//Prob+=0.4698*Po/L/L*(exp(-I*Dm2[k][l]*L/E*(2.*1.267))*(Ueconj[k]*Ue[l]+pxL*UeNSIconj[k]*Ue[l]+pxL*Ueconj[k]*UeNSI[l]+pxx*UeNSIconj[k]*UeNSI[l])*(Ueconj[l]*Ue[k]+dxL*UeNSI[k]*Ueconj[l]+dxL*UeNSIconj[l]*Ue[k]+dxx*UeNSI[k]*UeNSIconj[l]));
		Prob+=(exp(-I*Dm2[k][l]*L/E*(2.*1.267))*(Ueconj[k]*Ue[l]+pxL*UeNSIconj[k]*Ue[l]+pxL*Ueconj[k]*UeNSI[l]+pxx*UeNSIconj[k]*UeNSI[l])*(Ueconj[l]*Ue[k]+dxL*UeNSI[k]*Ueconj[l]+dxL*UeNSIconj[l]*Ue[k]+dxx*UeNSI[k]*UeNSIconj[l]));
		
	};
	};};
	

	P[0][0]=double( std::real(Prob))/2;


  return 0;
}


