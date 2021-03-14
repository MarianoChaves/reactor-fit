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

/* (C) 2005, 2007, 2019, 2020 Patrick Huber, J.M. Berryman */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <argp.h>
#include <vector>
#include <iostream>
#include <fstream>

#include <globes/globes.h>   /* GLoBES library */

#include "glf_types.h"
#include "glf_spectrum_chi.h"
#include "glf_precomputed_probabilities.h"
#include "glf_probability.h"
#include "nu_pre.h"
#include "Binned.h"
#include "emcee.h"
#include "experiments.h"

#define VERSION "1.0"
#define NUMP 6
double glf_systematic[4];


/* Program documentation. */
static char doc[] ="Data fitting for SBL experiments with GLoBES";

double sin2(double x){return sin(2*x);};
double cos2(double x){return cos(2*x);};



/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

int main(int argc, char *argv[])
{
  
  int i, s, s1, s2, s3, s4, s5;
  /* Initialize libglobes */
  glbInit(argv[0]);
  FILE *probabi = NULL;
  probabi = fopen("prob.csv","w");
  fprintf(probabi,"x1,x2,F1,F2\n");
  glbSetInitialStep(0.01);
  glbSelectMinimizer(GLB_MIN_NESTED_POWELL);

#ifndef  DEBUG

/*************************************************
*                                                *
*        DEFINING THE OSCILLATION ENGINES        *
*                                                *
**************************************************/
Pre pre_sin2AD1(&sin2);
Pre pre_sin2AD2(&sin2);
Pre pre_sin2AD3(&sin2);
Pre pre_sin2AD4(&sin2);
Pre pre_sin2AD5(&sin2);
Pre pre_sin2AD6(&sin2);
Pre pre_sin2AD7(&sin2);
Pre pre_sin2AD8(&sin2);

Pre pre_cos2AD1(&cos2);
Pre pre_cos2AD2(&cos2);
Pre pre_cos2AD3(&cos2);
Pre pre_cos2AD4(&cos2);
Pre pre_cos2AD5(&cos2);
Pre pre_cos2AD6(&cos2);
Pre pre_cos2AD7(&cos2);
Pre pre_cos2AD8(&cos2);

Pre pre_sin2ND_RN(&sin2);
Pre pre_sin2FD_RN(&sin2);

Pre pre_cos2ND_RN(&cos2);
Pre pre_cos2FD_RN(&cos2);

Pre pre_sin2ND_DC(&sin2);
Pre pre_sin2FD_DC(&sin2);

Pre pre_cos2ND_DC(&cos2);
Pre pre_cos2FD_DC(&cos2);

Pre pre_sin2_KL(&sin2);
Pre pre_cos2_KL(&cos2);

  


  /*************
  *  DAYA BAY  *
  **************/

Binned Q(0.1, 1000, 2000, "log");
std::vector<double> qf = Q.energy;

pre_cos2AD1.setQ(qf);
pre_cos2AD1.setGeo(&(geo::DB_AD1));
pre_sin2AD1.setQ(qf);
pre_sin2AD1.setGeo(&(geo::DB_AD1));
for(int i =0; i<=2000; i++){
  fprintf(probabi,"%f,%f,%f,%f\n",log10(qf[i]), DB_EH1_AD1_s.data[i][0], DB_EH1_AD1_s.data[i][1],pre_sin2AD1.F[i]);
  DB_EH1_AD1_s.data[i][0]=log10(qf[i]);
  DB_EH1_AD1_s.data[i][1]=pre_cos2AD1.F[i];
  DB_EH1_AD1_s.data2[i][0]=log10(qf[i]);
  DB_EH1_AD1_s.data2[i][1]=pre_sin2AD1.F[i];
}
glbDefineOscEngine(NUMP, &glf_nsi_probability_pre,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "DayaBay_EH1_AD1", (void *) &DB_EH1_AD1_s);


pre_cos2AD2.setQ(qf);
pre_cos2AD2.setGeo(&(geo::DB_AD2));
pre_sin2AD2.setQ(qf);
pre_sin2AD2.setGeo(&(geo::DB_AD2));
for(int i =0; i<=2000; i++){
  DB_EH1_AD2_s.data[i][0]=log10(qf[i]);
  DB_EH1_AD2_s.data[i][1]=pre_cos2AD2.F[i];
  DB_EH1_AD2_s.data2[i][0]=log10(qf[i]);
  DB_EH1_AD2_s.data2[i][1]=pre_sin2AD2.F[i];
}
glbDefineOscEngine(NUMP, &glf_nsi_probability_pre,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "DayaBay_EH1_AD2", (void *) &DB_EH1_AD2_s);


pre_cos2AD3.setQ(qf);
pre_cos2AD3.setGeo(&(geo::DB_AD3));
pre_sin2AD3.setQ(qf);
pre_sin2AD3.setGeo(&(geo::DB_AD3));
for(int i =0; i<=2000; i++){
  DB_EH2_AD3_s.data[i][0]=log10(qf[i]);
  DB_EH2_AD3_s.data[i][1]=pre_cos2AD3.F[i];
  DB_EH2_AD3_s.data2[i][0]=log10(qf[i]);
  DB_EH2_AD3_s.data2[i][1]=pre_sin2AD3.F[i];
}

glbDefineOscEngine(NUMP, &glf_nsi_probability_pre,
		   	&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "DayaBay_EH2_AD3", (void *) &DB_EH2_AD3_s);

pre_cos2AD8.setQ(qf);
pre_cos2AD8.setGeo(&(geo::DB_AD8));
pre_sin2AD8.setQ(qf);
pre_sin2AD8.setGeo(&(geo::DB_AD8));
for(int i =0; i<=2000; i++){
  DB_EH2_AD8_s.data[i][0]=log10(qf[i]);
  DB_EH2_AD8_s.data[i][1]=pre_cos2AD8.F[i];
  DB_EH2_AD8_s.data2[i][0]=log10(qf[i]);
  DB_EH2_AD8_s.data2[i][1]=pre_sin2AD8.F[i];
}
glbDefineOscEngine(NUMP, &glf_nsi_probability_pre,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "DayaBay_EH2_AD8", (void *) &DB_EH2_AD8_s);

pre_cos2AD4.setQ(qf);
pre_cos2AD4.setGeo(&(geo::DB_AD4));
pre_sin2AD4.setQ(qf);
pre_sin2AD4.setGeo(&(geo::DB_AD4));
for(int i =0; i<=2000; i++){
  DB_EH3_AD4_s.data[i][0]=log10(qf[i]);
  DB_EH3_AD4_s.data[i][1]=pre_cos2AD4.F[i];
  DB_EH3_AD4_s.data2[i][0]=log10(qf[i]);
  DB_EH3_AD4_s.data2[i][1]=pre_sin2AD4.F[i];
}
glbDefineOscEngine(NUMP, &glf_nsi_probability_pre,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "DayaBay_EH3_AD4", (void *) &DB_EH3_AD4_s);


pre_cos2AD5.setQ(qf);
pre_cos2AD5.setGeo(&(geo::DB_AD5));
pre_sin2AD5.setQ(qf);
pre_sin2AD5.setGeo(&(geo::DB_AD5));
for(int i =0; i<=2000; i++){
  DB_EH3_AD5_s.data[i][0]=log10(qf[i]);
  DB_EH3_AD5_s.data[i][1]=pre_cos2AD5.F[i];
  DB_EH3_AD5_s.data2[i][0]=log10(qf[i]);
  DB_EH3_AD5_s.data2[i][1]=pre_sin2AD5.F[i];
}
glbDefineOscEngine(NUMP, &glf_nsi_probability_pre,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "DayaBay_EH3_AD5", (void *) &DB_EH3_AD5_s);


pre_cos2AD6.setQ(qf);
pre_cos2AD6.setGeo(&(geo::DB_AD6));
pre_sin2AD6.setQ(qf);
pre_sin2AD6.setGeo(&(geo::DB_AD6));
for(int i =0; i<=2000; i++){
  DB_EH3_AD6_s.data[i][0]=log10(qf[i]);
  DB_EH3_AD6_s.data[i][1]=pre_cos2AD6.F[i];
  DB_EH3_AD6_s.data2[i][0]=log10(qf[i]);
  DB_EH3_AD6_s.data2[i][1]=pre_sin2AD6.F[i];
}
glbDefineOscEngine(NUMP, &glf_nsi_probability_pre,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "DayaBay_EH3_AD6", (void *) &DB_EH3_AD6_s);


pre_cos2AD7.setQ(qf);
pre_cos2AD7.setGeo(&(geo::DB_AD7));
pre_sin2AD7.setQ(qf);
pre_sin2AD7.setGeo(&(geo::DB_AD7));
for(int i =0; i<=2000; i++){
  DB_EH3_AD7_s.data[i][0]=log10(qf[i]);
  DB_EH3_AD7_s.data[i][1]=pre_cos2AD7.F[i];
  DB_EH3_AD7_s.data2[i][0]=log10(qf[i]);
  DB_EH3_AD7_s.data2[i][1]=pre_sin2AD7.F[i];
}
glbDefineOscEngine(NUMP, &glf_nsi_probability_pre,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "DayaBay_EH3_AD7", (void *) &DB_EH3_AD7_s);


  /*****************
  *  DOUBLE CHOOZ  *
  ******************/

Binned Q_DC(0.1, pow(10.0, 3), 2000, "log");
std::vector<double> qf_DC = Q_DC.energy;
pre_cos2ND_DC.setQ(qf_DC);
pre_cos2ND_DC.setGeo(&(geo::DC_ND));
pre_sin2ND_DC.setQ(qf_DC);
pre_sin2ND_DC.setGeo(&(geo::DC_ND));
for(int i =0; i<=2000; i++){
  DC_ND_s.data[i][0]=log10(qf_DC[i]);
  DC_ND_s.data[i][1]=pre_cos2ND_DC.F[i];
  DC_ND_s.data2[i][0]=log10(qf_DC[i]);
  DC_ND_s.data2[i][1]=pre_sin2ND_DC.F[i];
}

glbDefineOscEngine(NUMP, &glf_nsi_probability_pre,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "DC_ND", (void *) &DC_ND_s);


Binned Q2_DC(0.1, pow(10.0, 2.5), 1750, "log");
std::vector<double> qf2_DC = Q2_DC.energy;
pre_cos2FD_DC.setQ(qf2_DC);
pre_cos2FD_DC.setGeo(&(geo::DC_FD));
pre_sin2FD_DC.setQ(qf2_DC);
pre_sin2FD_DC.setGeo(&(geo::DC_FD));
for(int i =0; i<=1750; i++){
  DC_FD_s.data[i][0]=log10(qf2_DC[i]);
  DC_FD_s.data[i][1]=pre_cos2FD_DC.F[i];
  DC_FD_s.data2[i][0]=log10(qf2_DC[i]);
  DC_FD_s.data2[i][1]=pre_sin2FD_DC.F[i];

}

glbDefineOscEngine(NUMP, &glf_nsi_probability_pre,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "DC_FD", (void *) &DC_FD_s);


  /*********
  *  RENO  *
  **********/

Binned Q_RN(0.1, pow(10.0, 2.5), 1750, "log");
std::vector<double> qf_RN = Q_RN.energy;
pre_cos2ND_RN.setQ(qf_RN);
pre_cos2ND_RN.setGeo(&(geo::RENO_ND));
pre_sin2ND_RN.setQ(qf_RN);
pre_sin2ND_RN.setGeo(&(geo::RENO_ND));
for(int i =0; i<=1750; i++){
  RENO_ND_s.data[i][0]=log10(qf_RN[i]);
  RENO_ND_s.data[i][1]=pre_cos2ND_RN.F[i];
  RENO_ND_s.data2[i][0]=log10(qf_RN[i]);
  RENO_ND_s.data2[i][1]=pre_sin2ND_RN.F[i];
}
  
glbDefineOscEngine(NUMP, &glf_nsi_probability_pre,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "RENO_ND", (void *) &RENO_ND_s);

pre_cos2FD_RN.setQ(qf_RN);
pre_cos2FD_RN.setGeo(&(geo::RENO_FD));
pre_sin2FD_RN.setQ(qf_RN);
pre_sin2FD_RN.setGeo(&(geo::RENO_FD));
for(int i =0; i<=1750; i++){
  //fprintf(probabi,"%f,%f,%f,%f\n",log10(qf_RN[i]), RENO_FD_s.data[i][0], RENO_FD_s.data[i][1],pre_sin2FD_RN.F[i]);
  RENO_FD_s.data[i][0]=log10(qf_RN[i]);
  RENO_FD_s.data[i][1]=pre_cos2FD_RN.F[i];
  RENO_FD_s.data2[i][0]=log10(qf_RN[i]);
  RENO_FD_s.data2[i][1]=pre_sin2FD_RN.F[i];

}

glbDefineOscEngine(NUMP, &glf_nsi_probability_pre,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "RENO_FD", (void *) &RENO_FD_s);

  /************
  *  KAMLAND  *
  ************/



Binned Q_KL(0.1/1000, 1, 3500, "log");
std::vector<double> qf_KL = Q_KL.energy;
pre_cos2_KL.setQ(qf_KL);
pre_cos2_KL.setGeo(&(geo::KAMLAND));
pre_sin2_KL.setQ(qf_KL);
pre_sin2_KL.setGeo(&(geo::KAMLAND));
for(int i =0; i<=3500; i++){
  //fprintf(probabi,"%f,%f,%f,%f\n",log10(qf_KL[i]), RENO_FD_s.data[i][0], RENO_FD_s.data[i][1],pre_sin2_KL.F[i]);
  KAMLAND_s.data[i][0]=log10(qf_KL[i]);
  KAMLAND_s.data[i][1]=pre_cos2_KL.F[i];
  KAMLAND_s.data2[i][0]=log10(qf_KL[i]);
  KAMLAND_s.data2[i][1]=pre_sin2_KL.F[i];
}

glbDefineOscEngine(NUMP, &glf_nsi_probability_pre,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "KAMLAND", (void *) &KAMLAND_s);

printf("\n\n*******The probabilities have already been calculated******\n\n");
#endif

/***********************************************
*                                              *
*        DEFINING CHI-SQUARED FUNCTIONS        *
*                                              *
************************************************/

  int YesNo[10];

  glbDefineChiFunction(&glf_spectra_chi, 0, "spectra-chisq", (void *)YesNo);

/*************************************
*                                    *
*        THE MEAT OF THE CODE        *
*                                    *
**************************************/


/* Loading in experiment files... */

  s=glbInitExperiment("glb/spectra-DB-Near.glb",&glb_experiment_list[0],
	&glb_num_of_exps);  
  s2=glbInitExperiment("glb/spectra-DB-Far.glb",&glb_experiment_list[0],
	&glb_num_of_exps);
  s3=glbInitExperiment("glb/spectra-DC.glb",&glb_experiment_list[0],
	&glb_num_of_exps);
  s4=glbInitExperiment("glb/spectra-RN.glb",&glb_experiment_list[0],
	&glb_num_of_exps);
  s5=glbInitExperiment("glb/spectra-KL.glb",&glb_experiment_list[0],
	&glb_num_of_exps);

  /* Testing for failure */

  if(s<-1) {fprintf(stderr,"%s: FATAL: Unrecoverable parse error\n",argv[0]);exit(1);}
  if(s2<-1) {fprintf(stderr,"%s: FATAL: Unrecoverable parse error\n",argv[0]);exit(1);}
  if(s3<-1) {fprintf(stderr,"%s: FATAL: Unrecoverable parse error\n",argv[0]);exit(1);}
  if(s4<-1) {fprintf(stderr,"%s: FATAL: Unrecoverable parse error\n",argv[0]);exit(1);}
  if(s5<-1) {fprintf(stderr,"%s: FATAL: Unrecoverable parse error\n",argv[0]);exit(1);}
  
  /* Initialize parameter vector(s) */
  glb_params test_values = glbAllocParams();
  glb_params starting_values = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  glb_params minimum = glbAllocParams();

  glb_projection proj=glbAllocProjection();

  glbSetOscParams(test_values, asin(sqrt(0.084))/2.0, GLF_THETA_13);
  glbSetOscParams(test_values, 2.525e-3, GLF_DELTA_ATM);
  glbSetOscParams(test_values, 0.590, GLF_THETA_12);
  glbSetOscParams(test_values, 7.55e-5, GLF_DELTA_SOLAR);
  glbSetOscParams(test_values, 0.0, GLF_REAL);
  glbSetOscParams(test_values, 0.0, GLF_IMAG);

  //                        TH_13 0,   REAL 1,    IMAG 2,    D_ATM 3  , TH_12 4,  D_SOLAR 5
  glbDefineProjection(proj, GLB_FIXED, GLB_FIXED, GLB_FIXED, GLB_FIXED, GLB_FIXED, GLB_FIXED);
  glbSetDensityProjectionFlag(proj,GLB_FIXED,GLB_ALL);

  //glbSetOscParams(input_errors, 2.525e-3*0.1, GLF_THETA_13);
  glbSetOscParams(input_errors, asin(1), GLF_THETA_13);
  glbSetOscParams(input_errors, 0.0, GLF_DELTA_ATM);
  glbSetOscParams(input_errors, 0.0, GLF_THETA_12);
  glbSetOscParams(input_errors, 0.0, GLF_DELTA_SOLAR);
  glbSetOscParams(input_errors, 0.0, GLF_REAL);
  glbSetOscParams(input_errors, 0.0, GLF_IMAG);
  glbSetDensityParams(input_errors,0,GLB_ALL);
  glbSetInputErrors(input_errors);

  glbCopyParams(test_values,starting_values);
  glbSetStartingValues(starting_values);
  glbSetProjection(proj);

  /* The simulated data are computed */
  glbSetOscillationParameters(test_values);
  glbSetRates();

/*************************************************************************

  Initializing energy-scale uncertainty systematics...                 
  Also initializing the minimization over these nuisance parameters... 

*************************************************************************/



/* 
   We define this object now so that we can read out the ending values of
  the nuisance parameters...
*/

  double *EndNuisances;

/***************************************************************

  This block of code executes a scan over the parameters;

***************************************************************/
  FILE *output=NULL;
  output = fopen("output/output.csv","w");


//GRID SEARCH
  if(1==0){
      printf("\n\n*******Running parameter space******\n\n");

      fprintf(output,"th12,th13,dm31,dm21,real,imag,chi\n");

      /////////////////////////////////////////////////
      //////////////////////////////////////////////////
      double t12min, t12max, dt12;
      double t13min, t13max, dt13;
      double dreal, realmax, realmin;
      double dimag, imagmax, imagmin;
      double ddm, dmmax, dmmin;
      double ddm_s, dmmax_s, dmmin_s;

      int Nth12, i, Nth13, j, Nreal, k, Nimag, l, N;
      int Ndm, m, Ndm_s, m_s;

      t12min=log10(0.1);
      t12max=log10(10.);

      t13min=0.05;
      t13max=0.15;

      realmin=-1.0;
      realmax=+1.0;

      imagmin=-1.0;
      imagmax=+1.0;

      dmmax=3.5*(1.0e-3);
      dmmin=1.5*(1.0e-3);

      dmmax_s=10*(1.0e-5);
      dmmin_s=6*(1.0e-5);
      /////////////////////////////////////////////////
      /////////////////////////////////////////////////

      Nth12 = 20;
      Nth13 = 0;
      Nreal = 0;
      Nimag = 0;
      Ndm   = 0;
      Ndm_s = 20;

      dt12=(t12max-t12min)/Nth12;
      dt13=(t13max-t13min)/Nth13;
      dreal=(realmax-realmin)/Nreal;
      dimag=(imagmax-imagmin)/Nimag;
      ddm = (dmmax-dmmin)/Ndm;
      ddm_s = (dmmax_s-dmmin_s)/Ndm_s;

      double chi2;
      
      for(i = 0; i<=Nth12; i++){
        double t12test = t12min+dt12*i;
      for(j = 0; j<=Nth13; j++){
        double t13test=t13min+dt13*j;
      for(k = 0; k<=Nreal; k++){
        double realtest= realmin+dreal*k;
      for(l = 0; l<=Nimag; l++){
        double imagtest= imagmin+dimag*l;
      for(m = 0; m<=Ndm; m++){
        double dmtest = dmmin+ddm*m;
      for(m_s = 0; m_s<=Ndm_s; m_s++){
     
          double dmtest_s = dmmin_s+ddm_s*m_s;
          
          glbSetOscParams(test_values,atan(sqrt(pow(10.0,t12test))),GLF_THETA_12);
          //glbSetOscParams(test_values,asin(sqrt(t13test))/2,GLF_THETA_13);
          //glbSetOscParams(test_values,realtest,GLF_REAL);
          //glbSetOscParams(test_values,imagtest,GLF_IMAG);
          //glbSetOscParams(test_values,dmtest,GLF_DELTA_ATM);
          glbSetOscParams(test_values,dmtest_s,GLF_DELTA_SOLAR);
          
          std::vector<double> test_values{ asin(sqrt(0.084))/2.0, 2.525e-3, atan(sqrt(pow(10.0,t12test))), dmtest_s, 0, 0, 0, 0, 0, 0};
          
          chi2=my_spectra_chi(test_values);

          fprintf(output,"%f,%f,%f,%f,%f,%f,%f\n", t12test, t13test, dmtest*1000, dmtest_s*100000, realtest, imagtest, chi2);
          EndNuisances = glbGetSysStartingValuesListPtr(0, 0, GLB_ON);
      
      }}}}}}
  }


//SAMPLING SEARCH DayaBay
  if(1==1){
    printf("\n\n*******Running parameter space******\n\n");

    srand(42);

    int nwalkers = 1000;
    int ndim = 3;
    int bournout = 200;
    int nsteps = 200;


    std::vector< std::vector<double> > init_pos;

//std::vector<double> test_values{ asin(sqrt(0.084))/2.0, 2.525e-3, atan(sqrt(pow(10.0,t12test))), dmtest_s, 0, 0, 0, 0, 0, 0};
    for(int k = 0; k < nwalkers; k++){
        double th13_test = double(rand()% 1000+2000)/100000;
        //double th13_test = double(rand()% 13000+2000)/100000;
        double dm_atm_test = double(rand()% 20000+15000)/10000;
        double real_test = double(rand()% 1000-500)/1000;
        double imag_test = double(rand()% 1000+1)/1000;
        double pull_1 = double(rand()% 1000+1)/1000;
        double pull_2 = double(rand()% 1000+1)/1000;
        double pull_3 = double(rand()% 1000+1)/1000;
        double pull_4 = double(rand()% 1000+1)/1000;
        std::vector<double> pos{th13_test, dm_atm_test, real_test};
        //std::vector<double> pos{th12_test, dm_sol_test, real_test, imag_test, pull_1, pull_2, pull_3, pull_4};
        init_pos.push_back(pos);

    }
    std::cout<<"\n\n****Initial proposal done!****\n\n";
    char file_name[256] = "ReSet_state.csv"; 
    //nu::Emcee my_sample(nwalkers, ndim, init_pos);
    nu::Emcee my_sample(nwalkers, ndim, file_name);
    //my_sample.load_state(file_name);

    my_sample.run_mcmc(my_spectra_chi, bournout);
    std::cout<<"\n\n****Burnout done!****\n\n";
    my_sample.save_state(file_name);

    //my_sample.reset();

    //my_sample.load_state(file_name);
    //my_sample.run_mcmc(my_spectra_chi, nsteps);
    //std::cout<<"\n\n****EMCEE done!****\n\n";

    char file_sample_name[256] = "samples.csv";
    char header[256] = "th13,dm31,real,walker\n";
    my_sample.save_chain_walker(file_sample_name,header);
  }

/***************************************************************

  This block we calculate the number of events to plot

***************************************************************/

  //DAYBAY
  if(1==0){ 
    
      double numerator, numerator2, denominator;
      int j;
      
      double DayaBayF1[4] = {0.5678, 0.0761, 0.30075, 0.05545};
      double DayaBayF2[4] = {0.56605, 0.07615, 0.3021, 0.0556};
      double DayaBayF3[4] = {0.5618, 0.0761, 0.30665, 0.0553};
      double DayaBayF8[4] = {0.56345, 0.076, 0.3052, 0.05555};
      
      double DayaBayF4[4] = {0.559, 0.076, 0.310, 0.055};
      double DayaBayF5[4] = {0.559, 0.076, 0.310, 0.055};
      double DayaBayF6[4] = {0.559, 0.076, 0.310, 0.055};
      double DayaBayF7[4] = {0.552, 0.076, 0.315, 0.057};

      double t12test= atan(sqrt(0.4));
      double t13test = asin(0.16375);
      double dm31_test = 2.35e-3;
      double dm21_test = 8e-5;
      double realtest=0.0;
      double imagtest=0.0;

      glbSetOscParams(test_values, t12test, GLF_THETA_12);
      glbSetOscParams(test_values, t13test, GLF_THETA_13);
      glbSetOscParams(test_values, realtest, GLF_REAL);
      glbSetOscParams(test_values, imagtest, GLF_IMAG);
      glbSetOscParams(test_values, dm31_test, GLF_DELTA_ATM);
      glbSetOscParams(test_values, dm21_test, GLF_DELTA_SOLAR);
      
      glbSetRates();
      double chi2=glbChiSys(test_values,GLB_ALL,GLB_ALL);
      FILE *outputf;
      outputf = fopen("output/events-db.csv","w");
      fprintf(outputf,"data21,data31\n");
      int exp = 0;
      double deltaDB[52] = {0.0};
      double DB12_fudge = 0.9933;
      double DB13_fudge = 0.9973;
      /************
      *  EH1 AD1  *
      *************/
      double *EH1_AD1_U235 = glbGetSignalFitRatePtr(0, 0);
      double *EH1_AD1_U238 = glbGetSignalFitRatePtr(0, 1);
      double *EH1_AD1_Pu239 = glbGetSignalFitRatePtr(0, 2);
      double *EH1_AD1_Pu241 = glbGetSignalFitRatePtr(0, 3);
      
      double *EH1_AD1_U235_0 = glbGetRuleRatePtr(0, 0);
      double *EH1_AD1_U238_0 = glbGetRuleRatePtr(0, 1);
      double *EH1_AD1_Pu239_0 = glbGetRuleRatePtr(0, 2);
      double *EH1_AD1_Pu241_0 = glbGetRuleRatePtr(0, 3);

      /************
      *  EH1 AD2  *
      *************/
      
      double *EH1_AD2_U235 = glbGetSignalFitRatePtr(1, 0);
      double *EH1_AD2_U238 = glbGetSignalFitRatePtr(1, 1);
      double *EH1_AD2_Pu239 = glbGetSignalFitRatePtr(1, 2);
      double *EH1_AD2_Pu241 = glbGetSignalFitRatePtr(1, 3);

      double *EH1_AD2_U235_0 = glbGetRuleRatePtr(1, 0);
      double *EH1_AD2_U238_0 = glbGetRuleRatePtr(1, 1);
      double *EH1_AD2_Pu239_0 = glbGetRuleRatePtr(1, 2);
      double *EH1_AD2_Pu241_0 = glbGetRuleRatePtr(1, 3);

      /************
      *  EH2 AD3  *
      *************/

      double *EH2_AD3_U235 = glbGetSignalFitRatePtr(2, 0);
      double *EH2_AD3_U238 = glbGetSignalFitRatePtr(2, 1);
      double *EH2_AD3_Pu239 = glbGetSignalFitRatePtr(2, 2);
      double *EH2_AD3_Pu241 = glbGetSignalFitRatePtr(2, 3);

      /************
      *  EH2 AD8  *
      *************/

      double *EH2_AD8_U235 = glbGetSignalFitRatePtr(3, 0);
      double *EH2_AD8_U238 = glbGetSignalFitRatePtr(3, 1);
      double *EH2_AD8_Pu239 = glbGetSignalFitRatePtr(3, 2);
      double *EH2_AD8_Pu241 = glbGetSignalFitRatePtr(3, 3);

      /************
      *  EH3 AD4  *
      *************/

      double *EH3_AD4_U235 = glbGetSignalFitRatePtr(4, 0);
      double *EH3_AD4_U238 = glbGetSignalFitRatePtr(4, 1);
      double *EH3_AD4_Pu239 = glbGetSignalFitRatePtr(4, 2);
      double *EH3_AD4_Pu241 = glbGetSignalFitRatePtr(4, 3);

      /************
      *  EH3 AD5  *
      *************/

      double *EH3_AD5_U235 = glbGetSignalFitRatePtr(5, 0);
      double *EH3_AD5_U238 = glbGetSignalFitRatePtr(5, 1);
      double *EH3_AD5_Pu239 = glbGetSignalFitRatePtr(5, 2);
      double *EH3_AD5_Pu241 = glbGetSignalFitRatePtr(5, 3);

      /************
      *  EH3 AD6  *
      *************/

      double *EH3_AD6_U235 = glbGetSignalFitRatePtr(6, 0);
      double *EH3_AD6_U238 = glbGetSignalFitRatePtr(6, 1);
      double *EH3_AD6_Pu239 = glbGetSignalFitRatePtr(6, 2);
      double *EH3_AD6_Pu241 = glbGetSignalFitRatePtr(6, 3);

      /************
      *  EH3 AD7  *
      *************/

      double *EH3_AD7_U235 = glbGetSignalFitRatePtr(7, 0);
      double *EH3_AD7_U238 = glbGetSignalFitRatePtr(7, 1);
      double *EH3_AD7_Pu239 = glbGetSignalFitRatePtr(7, 2);
      double *EH3_AD7_Pu241 = glbGetSignalFitRatePtr(7, 3);

      int index;

      /*
      REMINDER: We need to take the appropriate combinations of the Daya Bay bins
      because of the fine binning with which we have calculated their spectra!

      As a reminder, the Daya Bay spectra are broken into 227 bins (!!), each of width
      0.05 MeV. We restructure them to match the binning at Daya Bay, as follows:

      Bins 0-11: Daya Bay Bin 0 (0.7 - 1.3 MeV prompt)
      Bins 12 + 5*(i-1) + (0-4): Daya Bay Bin i ( = 1-24 )
                      (1.3 + 0.25*(i-1) + (0-0.25) MeV prompt)
      Bins 132-225: Daya Bay Bin 25 (7.3 - 12.0 MeV prompt)

      Below, we calculate the total number of events at each AD -- we initialize these here
      */

      double TotalAD1[26] = {0.0};
      double TotalAD2[26] = {0.0};
      double TotalAD3[26] = {0.0};
      double TotalAD8[26] = {0.0};

      double TotalAD4[26] = {0.0};
      double TotalAD5[26] = {0.0};
      double TotalAD6[26] = {0.0};
      double TotalAD7[26] = {0.0};

      for (i=0; i<12; i++){
          TotalAD1[0] += EH1_AD1_U235[i]*DayaBayF1[0];
          TotalAD1[0] += EH1_AD1_U238[i]*DayaBayF1[1];
          TotalAD1[0] += EH1_AD1_Pu239[i]*DayaBayF1[2];
          TotalAD1[0] += EH1_AD1_Pu241[i]*DayaBayF1[3];

          TotalAD2[0] += EH1_AD2_U235[i]*DayaBayF2[0];
          TotalAD2[0] += EH1_AD2_U238[i]*DayaBayF2[1];
          TotalAD2[0] += EH1_AD2_Pu239[i]*DayaBayF2[2];
          TotalAD2[0] += EH1_AD2_Pu241[i]*DayaBayF2[3];

          TotalAD3[0] += EH2_AD3_U235[i]*DayaBayF3[0];
          TotalAD3[0] += EH2_AD3_U238[i]*DayaBayF3[1];
          TotalAD3[0] += EH2_AD3_Pu239[i]*DayaBayF3[2];
          TotalAD3[0] += EH2_AD3_Pu241[i]*DayaBayF3[3];

          TotalAD8[0] += EH2_AD8_U235[i]*DayaBayF8[0];
          TotalAD8[0] += EH2_AD8_U238[i]*DayaBayF8[1];
          TotalAD8[0] += EH2_AD8_Pu239[i]*DayaBayF8[2];
          TotalAD8[0] += EH2_AD8_Pu241[i]*DayaBayF8[3];

          TotalAD4[0] += EH3_AD4_U235[i]*DayaBayF4[0];
          TotalAD4[0] += EH3_AD4_U238[i]*DayaBayF4[1];
          TotalAD4[0] += EH3_AD4_Pu239[i]*DayaBayF4[2];
          TotalAD4[0] += EH3_AD4_Pu241[i]*DayaBayF4[3];

          TotalAD5[0] += EH3_AD5_U235[i]*DayaBayF5[0];
          TotalAD5[0] += EH3_AD5_U238[i]*DayaBayF5[1];
          TotalAD5[0] += EH3_AD5_Pu239[i]*DayaBayF5[2];
          TotalAD5[0] += EH3_AD5_Pu241[i]*DayaBayF5[3];

          TotalAD6[0] += EH3_AD6_U235[i]*DayaBayF6[0];
          TotalAD6[0] += EH3_AD6_U238[i]*DayaBayF6[1];
          TotalAD6[0] += EH3_AD6_Pu239[i]*DayaBayF6[2];
          TotalAD6[0] += EH3_AD6_Pu241[i]*DayaBayF6[3];

          TotalAD7[0] += EH3_AD7_U235[i]*DayaBayF7[0];
          TotalAD7[0] += EH3_AD7_U238[i]*DayaBayF7[1];
          TotalAD7[0] += EH3_AD7_Pu239[i]*DayaBayF7[2];
          TotalAD7[0] += EH3_AD7_Pu241[i]*DayaBayF7[3];
      }

      for (j=1; j<25; j++){
        for (i=0; i<5; i++){
          TotalAD1[j] += EH1_AD1_U235[12+5*(j-1)+i]*DayaBayF1[0];
          TotalAD1[j] += EH1_AD1_U238[12+5*(j-1)+i]*DayaBayF1[1];
          TotalAD1[j] += EH1_AD1_Pu239[12+5*(j-1)+i]*DayaBayF1[2];
          TotalAD1[j] += EH1_AD1_Pu241[12+5*(j-1)+i]*DayaBayF1[3];

          TotalAD2[j] += EH1_AD2_U235[12+5*(j-1)+i]*DayaBayF2[0];
          TotalAD2[j] += EH1_AD2_U238[12+5*(j-1)+i]*DayaBayF2[1];
          TotalAD2[j] += EH1_AD2_Pu239[12+5*(j-1)+i]*DayaBayF2[2];
          TotalAD2[j] += EH1_AD2_Pu241[12+5*(j-1)+i]*DayaBayF2[3];

          TotalAD3[j] += EH2_AD3_U235[12+5*(j-1)+i]*DayaBayF3[0];
          TotalAD3[j] += EH2_AD3_U238[12+5*(j-1)+i]*DayaBayF3[1];
          TotalAD3[j] += EH2_AD3_Pu239[12+5*(j-1)+i]*DayaBayF3[2];
          TotalAD3[j] += EH2_AD3_Pu241[12+5*(j-1)+i]*DayaBayF3[3];

          TotalAD8[j] += EH2_AD8_U235[12+5*(j-1)+i]*DayaBayF8[0];
          TotalAD8[j] += EH2_AD8_U238[12+5*(j-1)+i]*DayaBayF8[1];
          TotalAD8[j] += EH2_AD8_Pu239[12+5*(j-1)+i]*DayaBayF8[2];
          TotalAD8[j] += EH2_AD8_Pu241[12+5*(j-1)+i]*DayaBayF8[3];

          TotalAD4[j] += EH3_AD4_U235[12+5*(j-1)+i]*DayaBayF4[0];
          TotalAD4[j] += EH3_AD4_U238[12+5*(j-1)+i]*DayaBayF4[1];
          TotalAD4[j] += EH3_AD4_Pu239[12+5*(j-1)+i]*DayaBayF4[2];
          TotalAD4[j] += EH3_AD4_Pu241[12+5*(j-1)+i]*DayaBayF4[3];

          TotalAD5[j] += EH3_AD5_U235[12+5*(j-1)+i]*DayaBayF5[0];
          TotalAD5[j] += EH3_AD5_U238[12+5*(j-1)+i]*DayaBayF5[1];
          TotalAD5[j] += EH3_AD5_Pu239[12+5*(j-1)+i]*DayaBayF5[2];
          TotalAD5[j] += EH3_AD5_Pu241[12+5*(j-1)+i]*DayaBayF5[3];

          TotalAD6[j] += EH3_AD6_U235[12+5*(j-1)+i]*DayaBayF6[0];
          TotalAD6[j] += EH3_AD6_U238[12+5*(j-1)+i]*DayaBayF6[1];
          TotalAD6[j] += EH3_AD6_Pu239[12+5*(j-1)+i]*DayaBayF6[2];
          TotalAD6[j] += EH3_AD6_Pu241[12+5*(j-1)+i]*DayaBayF6[3];

          TotalAD7[j] += EH3_AD7_U235[12+5*(j-1)+i]*DayaBayF7[0];
          TotalAD7[j] += EH3_AD7_U238[12+5*(j-1)+i]*DayaBayF7[1];
          TotalAD7[j] += EH3_AD7_Pu239[12+5*(j-1)+i]*DayaBayF7[2];
          TotalAD7[j] += EH3_AD7_Pu241[12+5*(j-1)+i]*DayaBayF7[3];
        }
      }

      for (i=132; i<226; i++){
          TotalAD1[25] += EH1_AD1_U235[i]*DayaBayF1[0];
          TotalAD1[25] += EH1_AD1_U238[i]*DayaBayF1[1];
          TotalAD1[25] += EH1_AD1_Pu239[i]*DayaBayF1[2];
          TotalAD1[25] += EH1_AD1_Pu241[i]*DayaBayF1[3];

          TotalAD2[25] += EH1_AD2_U235[i]*DayaBayF2[0];
          TotalAD2[25] += EH1_AD2_U238[i]*DayaBayF2[1];
          TotalAD2[25] += EH1_AD2_Pu239[i]*DayaBayF2[2];
          TotalAD2[25] += EH1_AD2_Pu241[i]*DayaBayF2[3];

          TotalAD3[25] += EH2_AD3_U235[i]*DayaBayF3[0];
          TotalAD3[25] += EH2_AD3_U238[i]*DayaBayF3[1];
          TotalAD3[25] += EH2_AD3_Pu239[i]*DayaBayF3[2];
          TotalAD3[25] += EH2_AD3_Pu241[i]*DayaBayF3[3];

          TotalAD8[25] += EH2_AD8_U235[i]*DayaBayF8[0];
          TotalAD8[25] += EH2_AD8_U238[i]*DayaBayF8[1];
          TotalAD8[25] += EH2_AD8_Pu239[i]*DayaBayF8[2];
          TotalAD8[25] += EH2_AD8_Pu241[i]*DayaBayF8[3];

          TotalAD4[25] += EH3_AD4_U235[i]*DayaBayF4[0];
          TotalAD4[25] += EH3_AD4_U238[i]*DayaBayF4[1];
          TotalAD4[25] += EH3_AD4_Pu239[i]*DayaBayF4[2];
          TotalAD4[25] += EH3_AD4_Pu241[i]*DayaBayF4[3];

          TotalAD5[25] += EH3_AD5_U235[i]*DayaBayF5[0];
          TotalAD5[25] += EH3_AD5_U238[i]*DayaBayF5[1];
          TotalAD5[25] += EH3_AD5_Pu239[i]*DayaBayF5[2];
          TotalAD5[25] += EH3_AD5_Pu241[i]*DayaBayF5[3];

          TotalAD6[25] += EH3_AD6_U235[i]*DayaBayF6[0];
          TotalAD6[25] += EH3_AD6_U238[i]*DayaBayF6[1];
          TotalAD6[25] += EH3_AD6_Pu239[i]*DayaBayF6[2];
          TotalAD6[25] += EH3_AD6_Pu241[i]*DayaBayF6[3];

          TotalAD7[25] += EH3_AD7_U235[i]*DayaBayF7[0];
          TotalAD7[25] += EH3_AD7_U238[i]*DayaBayF7[1];
          TotalAD7[25] += EH3_AD7_Pu239[i]*DayaBayF7[2];
          TotalAD7[25] += EH3_AD7_Pu241[i]*DayaBayF7[3];
    }

    /*
    Now return to business as usual: the ratio between EH2/EH1 and EH3/EH1
    */

      for (i=0; i<26; i++){
          denominator = TotalAD1[i] + TotalAD2[i];
          numerator = TotalAD3[i] + TotalAD8[i];
          numerator2 = TotalAD4[i] + TotalAD5[i] + TotalAD6[i] + TotalAD7[i];
          deltaDB[i] = DB12_fudge*numerator/denominator ;
          deltaDB[i+26] = DB13_fudge*numerator2/denominator ;
          fprintf(outputf,"%f,%f\n",deltaDB[i],deltaDB[i+26]);
      }
  }
    
  //KAMLAND
  if(1==0){ 
    printf("\n\n*******Calculating KAMLAND EVENTS******\n\n");
    double numerator, numerator2, denominator;
    int j;
    
    double Kamland_F[4] = {0.567, 0.078, 0.298, 0.057};

    double t12test= atan(sqrt(0.366));
    double t13test = asin(sqrt(0.023));
    double dm31_test = 2.35e-3;
    double dm21_test = 8e-5;
    double realtest=0.0;
    double imagtest=0.0;


    glbSetOscParams(test_values, t12test, GLF_THETA_12);
    glbSetOscParams(test_values, t13test, GLF_THETA_13);
    glbSetOscParams(test_values, realtest, GLF_REAL);
    glbSetOscParams(test_values, imagtest, GLF_IMAG);
    glbSetOscParams(test_values, dm31_test, GLF_DELTA_ATM);
    glbSetOscParams(test_values, dm21_test, GLF_DELTA_SOLAR);
    
    glbSetRates();
    double chi2=glbChiSys(test_values,GLB_ALL,GLB_ALL);
    EndNuisances = glbGetSysStartingValuesListPtr(0, 0, GLB_ON);

    FILE *outputf;
    //outputf = fopen("kamland_output/expected.csv","w");
    outputf = fopen("kamland_output/events_bfp.csv","w");
    fprintf(outputf,"data,background\n");
    int exp = 0;
    double deltaKamland[17] = {0.0};
    double emin, emax;

    /************
    *  KAMLAND  *
    *************/
    glbGetEminEmax(12, &emin, &emax);


    double DeltaE0 = (emax-emin)/17;
    double t0 = emin/DeltaE0;
    double *Kamland_Simu_aux;
    double Kamland_Simu[29] = {0.0};
    double Kamland_Background[29] = {0.0};
  
    int linit = 6;
    double Kamland_BG[17] = {117.7189997, 132.3789758, 92.90521391, 68.05917575, 31.88233391, 6.627548887, 1.676206456, 1.492317065, 5.673098239, 3.782065492, 0.654213836, 11.489414, 11.27458771, 1.420708301, 0.255498155, 0, 0};
    
    Kamland_Simu_aux = glbGetSignalFitRatePtr(exp+12, 0);
    for (int l=0; l<17; l++){
      Kamland_Background[l+linit]=Kamland_BG[l];
      Kamland_Simu[l+linit]=Kamland_Simu_aux[l];
    };


    double totaltest = 0;
    for (i=0; i<17; i++){
        double a = EndNuisances[0]*0.01;
        double b = EndNuisances[1]*0.01;

        double del = b*(i*t0+0.5)+i;
        int kd = floor(del);

        numerator = (1+b)*(1+a)*( (Kamland_Simu[linit+kd+1]-Kamland_Simu[linit+kd])*(del-kd)+ Kamland_Simu[linit+kd]);

        a = EndNuisances[2]*0.01;
        b = EndNuisances[3]*0.01;
        double background_rate = (1+b)*(1+a)*( (Kamland_Background[linit+kd+1]-Kamland_Background[linit+kd])*(del-kd)+ Kamland_Background[linit+kd]);
        totaltest+=numerator;

        fprintf(outputf,"%f, %f\n",numerator, background_rate);
    }
    printf("Total Rate: %f\n",totaltest);
}

  //RENO
  if(1==0){

    double numerator, numerator2, denominator;
    int j;


    double t13test = asin(sqrt(0.0896))/2;
    double t12test= atan(sqrt(0.436));
    double dm31_test = 2.75e-3;
    double dm21_test = 7.53e-5;
    double realtest=0.0;
    double imagtest=0.0;

    glbSetOscParams(test_values, t12test, GLF_THETA_12);
    glbSetOscParams(test_values, t13test, GLF_THETA_13);
    glbSetOscParams(test_values, realtest, GLF_REAL);
    glbSetOscParams(test_values, imagtest, GLF_IMAG);
    glbSetOscParams(test_values, dm31_test, GLF_DELTA_ATM);
    glbSetOscParams(test_values, dm21_test, GLF_DELTA_SOLAR);

    glbSetRates();
    double chi2=glbChiSys(test_values,GLB_ALL,GLB_ALL);
    FILE *outputf;
    outputf = fopen("output/events-rn.csv","w");
    fprintf(outputf,"data\n");
    int exp = 0;

    int index;

    /*********
    *  RENO  *
    **********/

    /*
    We introduce a fudge factor for RENO in order to reproduce their determinations
    of theta_13 and Dm31/Dmee
    */
    double RENO_fudge = 1.00884;
    double deltaRENO[25] = {0.0};
    double TotND, TotFD;
    double RENO_Ratio[25] = {0.12592, 0.120914, 0.11629, 0.12038, 0.1154, 0.117704, 0.115853, 0.116933, 0.117296, 0.115695, 0.118126, 0.118856, 0.120047, 0.12019, 0.117643, 0.118033, 0.118237, 0.119066, 0.121666, 0.119507, 0.115003, 0.122545, 0.119143, 0.121214, 0.124873};

    double RENO_ND_F[4] = {0.57264, 0.07309, 0.29911, 0.055161};
    double RENO_FD_F[4] = {0.57447, 0.07340, 0.29742, 0.054752};

    /******************
    *  NEAR DETECTOR  *
    *******************/

    double *RENO_ND_U235 = glbGetSignalFitRatePtr(10, 0);
    double *RENO_ND_U238 = glbGetSignalFitRatePtr(10, 1);
    double *RENO_ND_Pu239 = glbGetSignalFitRatePtr(10, 2);
    double *RENO_ND_Pu241 = glbGetSignalFitRatePtr(10, 3);

    /*****************
    *  FAR DETECTOR  *
    ******************/

    double *RENO_FD_U235 = glbGetSignalFitRatePtr(11, 0);
    double *RENO_FD_U238 = glbGetSignalFitRatePtr(11, 1);
    double *RENO_FD_Pu239 = glbGetSignalFitRatePtr(11, 2);
    double *RENO_FD_Pu241 = glbGetSignalFitRatePtr(11, 3);


    for (i=0; i<22; i++){ /* The first 22 bins are totally normal...  */
    TotND = RENO_ND_U235[i]* RENO_ND_F[0];
    TotND += RENO_ND_U238[i]* RENO_ND_F[1];
    TotND += RENO_ND_Pu239[i]* RENO_ND_F[2];
    TotND += RENO_ND_Pu241[i]* RENO_ND_F[3];

    TotFD = RENO_FD_U235[i]* RENO_FD_F[0];
    TotFD += RENO_FD_U238[i]* RENO_FD_F[1];
    TotFD += RENO_FD_Pu239[i]* RENO_FD_F[2];
    TotFD += RENO_FD_Pu241[i]* RENO_FD_F[3];

    deltaRENO[i] = TotFD/TotND;
    fprintf(outputf,"%f\n",deltaRENO[i]);
  }

  for (i=0; i<2; i++){ /* The next 2 bins are actually two bins combined...  */
    TotND = RENO_ND_U235[22+2*i]* RENO_ND_F[0];
    TotND += RENO_ND_U238[22+2*i]* RENO_ND_F[1];
    TotND += RENO_ND_Pu239[22+2*i]* RENO_ND_F[2];
    TotND += RENO_ND_Pu241[22+2*i]* RENO_ND_F[3];

    TotND += RENO_ND_U235[22+2*i+1]* RENO_ND_F[0];
    TotND += RENO_ND_U238[22+2*i+1]* RENO_ND_F[1];
    TotND += RENO_ND_Pu239[22+2*i+1]* RENO_ND_F[2];
    TotND += RENO_ND_Pu241[22+2*i+1]* RENO_ND_F[3];

    TotFD = RENO_FD_U235[22+2*i]* RENO_FD_F[0];
    TotFD += RENO_FD_U238[22+2*i]* RENO_FD_F[1];
    TotFD += RENO_FD_Pu239[22+2*i]* RENO_FD_F[2];
    TotFD += RENO_FD_Pu241[22+2*i]* RENO_FD_F[3];

    TotFD += RENO_FD_U235[22+2*i+1]* RENO_FD_F[0];
    TotFD += RENO_FD_U238[22+2*i+1]* RENO_FD_F[1];
    TotFD += RENO_FD_Pu239[22+2*i+1]* RENO_FD_F[2];
    TotFD += RENO_FD_Pu241[22+2*i+1]* RENO_FD_F[3];

    deltaRENO[i+22] = RENO_fudge*TotFD/TotND;
    fprintf(outputf,"%f\n",deltaRENO[i+22]);
  }

  TotND = 0.0;
  TotFD = 0.0;

  for (i=26; i<29; i++){ /* The last bin is the sum of three bins  */
    TotND += RENO_ND_U235[i]* RENO_ND_F[0];
    TotND += RENO_ND_U238[i]* RENO_ND_F[1];
    TotND += RENO_ND_Pu239[i]* RENO_ND_F[2];
    TotND += RENO_ND_Pu241[i]* RENO_ND_F[3];

    TotFD += RENO_FD_U235[i]* RENO_FD_F[0];
    TotFD += RENO_FD_U238[i]* RENO_FD_F[1];
    TotFD += RENO_FD_Pu239[i]* RENO_FD_F[2];
    TotFD += RENO_FD_Pu241[i]* RENO_FD_F[3];
  }
  deltaRENO[24] = RENO_fudge*TotFD/TotND;
  fprintf(outputf,"%f\n",deltaRENO[24]);
}	  

  //DOUBLE CHOOZ
  if(1==0){
                
    double numerator, numerator2, denominator;
    int j;
    
    double DC_Ratio[26] = {0.425358, 0.422108, 0.430931, 0.435575, 0.438825, 0.437896, 0.445791, 0.44254, 0.441147, 0.444397, 0.444862, 0.434181, 0.443933, 0.450899, 0.455542, 0.463436, 0.455078, 0.470402, 0.45322, 0.458793, 0.451363, 0.44997, 0.465294, 0.474581, 0.469938, 0.520089};

    double DCF[4] = {0.511, 0.087, 0.340, 0.062};

    double t13test = asin(sqrt(0.105))/2;
    double t12test= atan(sqrt(0.436));
    double dm31_test = 2.25e-3;
    double dm21_test = 7.53e-5;
    double realtest=0.0;
    double imagtest=0.0;

    glbSetOscParams(test_values, t12test, GLF_THETA_12);
    glbSetOscParams(test_values, t13test, GLF_THETA_13);
    glbSetOscParams(test_values, realtest, GLF_REAL);
    glbSetOscParams(test_values, imagtest, GLF_IMAG);
    glbSetOscParams(test_values, dm31_test, GLF_DELTA_ATM);
    glbSetOscParams(test_values, dm21_test, GLF_DELTA_SOLAR);
    
    glbSetRates();
    double chi2=glbChiSys(test_values,GLB_ALL,GLB_ALL);
    FILE *outputf;
    outputf = fopen("output/events_imSet_-08.csv","w");
    fprintf(outputf,"data\n");
    int exp = 0;
    double deltaDC[26] = {0.0};
    /*****************
    *  DOUBLE CHOOZ  *
    ******************/

    /*
      We introduce a fudge factor for Double Chooz in order to reproduce their 
      determinations of theta_13 and Dm31/Dmee
    */
    double DC_fudge = 1.0026;

      /******************
      *  NEAR DETECTOR  *
      *******************/

    double *DC_ND_U235 = glbGetSignalFitRatePtr(8, 0);
    double *DC_ND_U238 = glbGetSignalFitRatePtr(8, 1);
    double *DC_ND_Pu239 = glbGetSignalFitRatePtr(8, 2);
    double *DC_ND_Pu241 = glbGetSignalFitRatePtr(8, 3);
      
      /*****************
      *  FAR DETECTOR  *
      ******************/

    double *DC_FD_U235 = glbGetSignalFitRatePtr(9, 0);
    double *DC_FD_U238 = glbGetSignalFitRatePtr(9, 1);
    double *DC_FD_Pu239 = glbGetSignalFitRatePtr(9, 2);
    double *DC_FD_Pu241 = glbGetSignalFitRatePtr(9, 3);

    double TotND, TotFD;


    for (i=0; i<26; i++){
      TotND = DC_ND_U235[i]*DCF[0];
      TotND += DC_ND_U238[i]*DCF[1];
      TotND += DC_ND_Pu239[i]*DCF[2];
      TotND += DC_ND_Pu241[i]*DCF[3];

      TotFD = DC_FD_U235[i]*DCF[0];
      TotFD += DC_FD_U238[i]*DCF[1];
      TotFD += DC_FD_Pu239[i]*DCF[2];
      TotFD += DC_FD_Pu241[i]*DCF[3];

      deltaDC[i] = DC_fudge*TotFD/TotND;
      fprintf(outputf,"%f\n",deltaDC[i]);
    }

  }


  /* Destroy parameter vector(s) */
  glbFreeParams(test_values); 
  glbFreeParams(minimum); 
  fclose(output);
  exit(0);

}
