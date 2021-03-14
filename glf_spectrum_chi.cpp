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
#include <limits>
#include <math.h>
#include <string.h>
#include <gsl/gsl_sf_expint.h>
#include <globes/globes.h>   /* GLoBES library */

#include "glf_spectrum_chi.h"
#include "glf_spectrum_aux_1.h"
#include "glf_spectrum_aux_2.h"

extern double glf_systematic[4];

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

  /**************************
  *  OTHER RELEVANT INPUTS  *
  ***************************/

/*
  The order of the experiments in the corresponding .glb file is:

  0 - Daya Bay EH1 AD1 (special format to be compatible with NEOS)
  1 - Daya Bay EH1 AD2 (special format to be compatible with NEOS)
  2 - Daya Bay EH2 AD3 (special format to be compatible with NEOS)
  3 - Daya Bay EH2 AD8 (special format to be compatible with NEOS)

  4 - Daya Bay EH3 AD4
  5 - Daya Bay EH3 AD5
  6 - Daya Bay EH3 AD6
  7 - Daya Bay EH3 AD7

  8 - Double Chooz ND
  9 - Double Chooz FD

  10 - RENO ND
  11 - RENO FD


*/

/**************************************************
*                                                 *
*        THE CALCULATION OF THE CHI-SQUARED       *
*                                                 *
***************************************************/
double glf_spectra_chi_old(int exp, int rule, int n_params, double *x, double *errors,
                 void *user_data)
{ 

/* 
  Here is where I calculate the chi-squared following

  There will be a "delta" for each experiment; assume no correlated uncertainties

  "user_data" contains flags indicated which experiments are to be included in the fit

*/

  double deltaDB[52] = {0.0};
  double deltaDC[26] = {0.0};
  double deltaRENO[25] = {0.0};
  double deltaKamland[17] = {0.0};

  int i, j, k;
  double numerator, numerator2, denominator;

  int Nnuisance=4; /* The number of nuisance parameters */

  int CheckZeros;
  /* A quantity that we use to look for zero entries after shifting energy scales */

  /* Various important quantities used below */
  double emin, emax; 
  /* Get emin/emax for given experiment; used for energy-scale variation */

  /************
  *  Kamland  *
  ************/

  /* Get emin and emax for KAMLAND */
  glbGetEminEmax(exp+12, &emin, &emax);

  double DeltaE0 = (emax-emin)/17;
  double t0 = emin/DeltaE0;
  double *Kamland_Simu_aux;
  double Kamland_Simu[29] = {0.0};
  double Kamland_Background[29] = {0.0};
  
  int linit = 6;

  Kamland_Simu_aux = glbGetSignalFitRatePtr(exp+12, 0);
  for (int l=0; l<17; l++){
    Kamland_Background[l+linit]=Kamland_BG[l];
    Kamland_Simu[l+linit]=Kamland_Simu_aux[l];
  };

//  glbShiftEnergyScale(0.01*x[1]*1.0, glbGetSignalFitRatePtr(exp+12, 0),
//                      Kamland_Simu, 17, emin, emax);  
//    glbShiftEnergyScale(0.01*x[3]*1.0, Kamland_Background_aux,
//                      Kamland_Background, 17, emin, emax);





CheckZeros = 0;
double totaltest = 0;
for (i=0; i<17; i++){
  double a = x[0]*0.01;
  double b = x[1]*0.01;

  double del = b*(i*t0+0.5)+i;
  int kd = floor(del);
//  printf("%f\n", del-kd);


  numerator = (1+b)*(1+a)*( (Kamland_Simu[linit+kd+1]-Kamland_Simu[linit+kd])*(del-kd)+ Kamland_Simu[linit+kd]);

  a = x[2]*0.01;
  b = x[3]*0.01;
  double background_rate = (1+b)*(1+a)*( (Kamland_Background[linit+kd+1]-Kamland_Background[linit+kd])*(del-kd)+ Kamland_Background[linit+kd]);
 



  deltaKamland[i] = Kamland_SG[i]-background_rate - numerator;
}

  /*************
  *  DAYA BAY  *
  **************/

  /*
    We introduce a couple of fudge factors for Daya Bay in order to reproduce their 
    determinations of theta_13 and Dm31(Dmee)
  */
  double DB12_fudge = 0.9933;
  double DB13_fudge = 0.9973;

/*
  WATCH CAREFULLY HERE! The four near detectors at Daya Bay are also used to determine
  the normalization for NEOS, so they are formatted differently from the four far
  detectors! Namely, while the far detectors are binned in 0.2 MeV bins from 2.08-8.68
  MeV, the near detectors are binned in 0.1 MeV bins from 1.78-8.68 MeV bins. The finer
  resolution is driven by the binning used at NEOS; just combine the appropriate bins to
  get the numbers of events useful at Daya Bay!

  I've also included one additional bin at low energies for the far detectors, which I'm
  not including in the calculation, because it gives a different result for the lowest bin
  that I am including. I'm not sure why this is the case...

  I also need the three-neutrino numbers for the near detectors; these are to be used with
  the NEOS analysis.
*/
/*
    //
    //  EH3 AD1  *
    //
  double *EH1_AD1_U235 = glbGetSignalFitRatePtr(exp+0, 0);
  double *EH1_AD1_U238 = glbGetSignalFitRatePtr(exp+0, 1);
  double *EH1_AD1_Pu239 = glbGetSignalFitRatePtr(exp+0, 2);
  double *EH1_AD1_Pu241 = glbGetSignalFitRatePtr(exp+0, 3);

  double *EH1_AD1_U235_0 = glbGetRuleRatePtr(exp+0, 0);
  double *EH1_AD1_U238_0 = glbGetRuleRatePtr(exp+0, 1);
  double *EH1_AD1_Pu239_0 = glbGetRuleRatePtr(exp+0, 2);
  double *EH1_AD1_Pu241_0 = glbGetRuleRatePtr(exp+0, 3);

    //
    //  EH3 AD2  *
    //

  double *EH1_AD2_U235 = glbGetSignalFitRatePtr(exp+1, 0);
  double *EH1_AD2_U238 = glbGetSignalFitRatePtr(exp+1, 1);
  double *EH1_AD2_Pu239 = glbGetSignalFitRatePtr(exp+1, 2);
  double *EH1_AD2_Pu241 = glbGetSignalFitRatePtr(exp+1, 3);

  double *EH1_AD2_U235_0 = glbGetRuleRatePtr(exp+1, 0);
  double *EH1_AD2_U238_0 = glbGetRuleRatePtr(exp+1, 1);
  double *EH1_AD2_Pu239_0 = glbGetRuleRatePtr(exp+1, 2);
  double *EH1_AD2_Pu241_0 = glbGetRuleRatePtr(exp+1, 3);

    //
    //  EH3 AD3  *
    //

  double *EH2_AD3_U235 = glbGetSignalFitRatePtr(exp+2, 0);
  double *EH2_AD3_U238 = glbGetSignalFitRatePtr(exp+2, 1);
  double *EH2_AD3_Pu239 = glbGetSignalFitRatePtr(exp+2, 2);
  double *EH2_AD3_Pu241 = glbGetSignalFitRatePtr(exp+2, 3);
   
    //
    //  EH3 AD8  *
    //

  double *EH2_AD8_U235 = glbGetSignalFitRatePtr(exp+3, 0);
  double *EH2_AD8_U238 = glbGetSignalFitRatePtr(exp+3, 1);
  double *EH2_AD8_Pu239 = glbGetSignalFitRatePtr(exp+3, 2);
  double *EH2_AD8_Pu241 = glbGetSignalFitRatePtr(exp+3, 3);

    //
    //  EH3 AD5  *
    //

  double *EH3_AD4_U235 = glbGetSignalFitRatePtr(exp+4, 0);
  double *EH3_AD4_U238 = glbGetSignalFitRatePtr(exp+4, 1);
  double *EH3_AD4_Pu239 = glbGetSignalFitRatePtr(exp+4, 2);
  double *EH3_AD4_Pu241 = glbGetSignalFitRatePtr(exp+4, 3);
   
    //
    //  EH3 AD5  *
    //

  double *EH3_AD5_U235 = glbGetSignalFitRatePtr(exp+5, 0);
  double *EH3_AD5_U238 = glbGetSignalFitRatePtr(exp+5, 1);
  double *EH3_AD5_Pu239 = glbGetSignalFitRatePtr(exp+5, 2);
  double *EH3_AD5_Pu241 = glbGetSignalFitRatePtr(exp+5, 3);

    //
    //  EH3 AD6  *
    //

  double *EH3_AD6_U235 = glbGetSignalFitRatePtr(exp+6, 0);
  double *EH3_AD6_U238 = glbGetSignalFitRatePtr(exp+6, 1);
  double *EH3_AD6_Pu239 = glbGetSignalFitRatePtr(exp+6, 2);
  double *EH3_AD6_Pu241 = glbGetSignalFitRatePtr(exp+6, 3);
   
    //
    //  EH3 AD7  *
    //

  double *EH3_AD7_U235 = glbGetSignalFitRatePtr(exp+7, 0);
  double *EH3_AD7_U238 = glbGetSignalFitRatePtr(exp+7, 1);
  double *EH3_AD7_Pu239 = glbGetSignalFitRatePtr(exp+7, 2);
  double *EH3_AD7_Pu241 = glbGetSignalFitRatePtr(exp+7, 3);

  int index;


//  REMINDER: We need to take the appropriate combinations of the Daya Bay bins
//  because of the fine binning with which we have calculated their spectra!

//  As a reminder, the Daya Bay spectra are broken into 227 bins (!!), each of width
//  0.05 MeV. We restructure them to match the binning at Daya Bay, as follows:

//    Bins 0-11: Daya Bay Bin 0 (0.7 - 1.3 MeV prompt)
//    Bins 12 + 5*(i-1) + (0-4): Daya Bay Bin i ( = 1-24 )
//                                (1.3 + 0.25*(i-1) + (0-0.25) MeV prompt)
//    Bins 132-225: Daya Bay Bin 25 (7.3 - 12.0 MeV prompt)

//  Below, we calculate the total number of events at each AD -- we initialize these here


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


  //Now return to business as usual: the ratio between EH2/EH1 and EH3/EH1

  for (i=0; i<26; i++){
    denominator = TotalAD1[i] + TotalAD2[i];
    numerator = TotalAD3[i] + TotalAD8[i];
    numerator2 = TotalAD4[i] + TotalAD5[i] + TotalAD6[i] + TotalAD7[i];

    deltaDB[i] = DayaBayRatio21[i] - DB12_fudge*numerator/denominator ;
    deltaDB[i+26] = DayaBayRatio31[i] - DB13_fudge*numerator2/denominator ;
    
  }

  //
  //  DOUBLE CHOOZ  *
  //
  
  
  //  We introduce a fudge factor for Double Chooz in order to reproduce their 
  //  determinations of theta_13 and Dm31/Dmee
  
  double DC_fudge = 1.0026;

    //
    //  NEAR DETECTOR  *
    //

  double *DC_ND_U235 = glbGetSignalFitRatePtr(exp+8, 0);
  double *DC_ND_U238 = glbGetSignalFitRatePtr(exp+8, 1);
  double *DC_ND_Pu239 = glbGetSignalFitRatePtr(exp+8, 2);
  double *DC_ND_Pu241 = glbGetSignalFitRatePtr(exp+8, 3);
   
    //
    //  FAR DETECTOR  *
    //

  double *DC_FD_U235 = glbGetSignalFitRatePtr(exp+9, 0);
  double *DC_FD_U238 = glbGetSignalFitRatePtr(exp+9, 1);
  double *DC_FD_Pu239 = glbGetSignalFitRatePtr(exp+9, 2);
  double *DC_FD_Pu241 = glbGetSignalFitRatePtr(exp+9, 3);

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

    deltaDC[i] = DC_Ratio[i] - DC_fudge*TotFD/TotND;
  }
  
  //
  //  RENO  *
  //

  
  //  We introduce a fudge factor for RENO in order to reproduce their determinations
  //  of theta_13 and Dm31/Dmee
  
  double RENO_fudge = 1.00884;

    //
    //  NEAR DETECTOR  *
    //

  double *RENO_ND_U235 = glbGetSignalFitRatePtr(exp+10, 0);
  double *RENO_ND_U238 = glbGetSignalFitRatePtr(exp+10, 1);
  double *RENO_ND_Pu239 = glbGetSignalFitRatePtr(exp+10, 2);
  double *RENO_ND_Pu241 = glbGetSignalFitRatePtr(exp+10, 3);
   
    //
    //  FAR DETECTOR  *
    //

  double *RENO_FD_U235 = glbGetSignalFitRatePtr(exp+11, 0);
  double *RENO_FD_U238 = glbGetSignalFitRatePtr(exp+11, 1);
  double *RENO_FD_Pu239 = glbGetSignalFitRatePtr(exp+11, 2);
  double *RENO_FD_Pu241 = glbGetSignalFitRatePtr(exp+11, 3);

  for (i=0; i<22; i++){ // The first 22 bins are totally normal...  
    TotND = RENO_ND_U235[i]* RENO_ND_F[0];
    TotND += RENO_ND_U238[i]* RENO_ND_F[1];
    TotND += RENO_ND_Pu239[i]* RENO_ND_F[2];
    TotND += RENO_ND_Pu241[i]* RENO_ND_F[3];

    TotFD = RENO_FD_U235[i]* RENO_FD_F[0];
    TotFD += RENO_FD_U238[i]* RENO_FD_F[1];
    TotFD += RENO_FD_Pu239[i]* RENO_FD_F[2];
    TotFD += RENO_FD_Pu241[i]* RENO_FD_F[3];

    deltaRENO[i] = RENO_Ratio[i] - RENO_fudge*TotFD/TotND;
  }

  for (i=0; i<2; i++){ // The next 2 bins are actually two bins combined...  
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

    deltaRENO[i+22] = RENO_Ratio[i+22] - RENO_fudge*TotFD/TotND;
  }

  TotND = 0.0;
  TotFD = 0.0;

  for (i=26; i<29; i++){ // The last bin is the sum of three bins  
    TotND += RENO_ND_U235[i]* RENO_ND_F[0];
    TotND += RENO_ND_U238[i]* RENO_ND_F[1];
    TotND += RENO_ND_Pu239[i]* RENO_ND_F[2];
    TotND += RENO_ND_Pu241[i]* RENO_ND_F[3];

    TotFD += RENO_FD_U235[i]* RENO_FD_F[0];
    TotFD += RENO_FD_U238[i]* RENO_FD_F[1];
    TotFD += RENO_FD_Pu239[i]* RENO_FD_F[2];
    TotFD += RENO_FD_Pu241[i]* RENO_FD_F[3];
  }
  deltaRENO[24] = RENO_Ratio[24] - RENO_fudge*TotFD/TotND;
*/

  /********************************
  *  Putting the pieces together  *
  *********************************/

  /* Adding the contributions from each experiment */

double chi2 = 0.0;

// Kamland
for (i=0; i<17; i++){
    chi2 += deltaKamland[i]*deltaKamland[i]/Kamland_SG[i];// *KamlandCovariance[i][j];
}
/*
// Daya Bay piece (EH2/EH1 and EH3/EH1) 
for (i=0; i<52; i++){
  for (j=0; j<52; j++){
    chi2 += deltaDB[i]*deltaDB[j]*DayaBayCovariance[i][j];
  }
}

// Double Chooz
for (i=0; i<26; i++){
  for (j=0; j<26; j++){
    chi2 += deltaDC[i]* deltaDC[j]*DC_Covariance[i][j];
  }
}
//RENO 
for (i=0; i<25; i++){
  for (j=0; j<25; j++){
    chi2 += deltaRENO[i]*deltaRENO[j]*RENO_Covariance[i][j];
  }
}
*/

/*
  Adding nuisance parameters; helps insure convergence of minimization to include
  nuisance parameters for experiments that are not included in fit...
*/

  for (i=0; i<Nnuisance; i++){
    //if(fabs(x[i]) > errors[i])
      chi2 += square( (fabs(x[i]) ) / errors[i]);
    //if(fabs(x[i]) < errors[i]) printf("%f\t",x[i]);
  }
//  printf("\n");  

  /* Save the systematics parameters as starting values for the next step */
  //printf("%f\n",x[0]);
  for (i=0; i < Nnuisance; i++)
    glf_systematic[i] = x[i];

  return chi2;
}  

double glf_spectra_chi(int exp, int rule, int n_params, double *x, double *errors,
                 void *user_data){
                   return 0;
                 }

//#######################################################
//#######################################################

/*
#define GLF_THETA_13 0
#define GLF_REAL 1
#define GLF_IMAG 2
#define GLF_DELTA_ATM 3
#define GLF_THETA_12 4
#define GLF_DELTA_SOLAR 5
*/
double my_spectra_chi(std::vector<double> test_params)
{ 

  glb_params test_values = glbAllocParams();

  double th13_t = asin(sqrt(0.084))/2.0;
  double dm_atm_t = 2.525e-3;
  double th12_t = 0.590;
  double dm_sol_t = 7.55e-5;
  double real_t = 0.0;
  double imag_t = 0.0;

  /*GLOBAL*/
  /*
  th13_t=test_params[0];
  dm_atm_t=test_params[1];
  th12_t=test_params[2];
  dm_sol_t =test_params[3];
  real_t = test_params[4];
  imag_t = test_params[5];
  double x[4] = {test_params[6],test_params[7],test_params[8],test_params[9]};
  */

 /* DayaBay + DC+RN*/
  if((test_params[0]<0) || (test_params[0] >1)){
    double infinity = std::numeric_limits<double>::infinity();
    return -infinity;
  };
  if(test_params[1]<0 ){
    double infinity = std::numeric_limits<double>::infinity();
    return -infinity;
  };

  th13_t=asin(sqrt(test_params[0]));
  dm_atm_t=test_params[1]*(1.0e-3);
  real_t = test_params[2];
  double x[4] = {0};
  
  
  /*KAMLAND*/
  /*
  if((test_params[0]<-4) || (test_params[0] >10)){
    return -1000000000000;
  };

  if((test_params[1]<5) || (test_params[1] >10)){
    return -10000000000000;
  };

  th12_t= atan(sqrt(pow(10.0,test_params[0])));
  dm_sol_t =test_params[1]*(1.0e-5);
  //real_t = test_params[2];
  //imag_t = test_params[3];
  double x[4] = {test_params[2],test_params[3],test_params[4],test_params[5]};
  */
 
  glbSetOscParams(test_values, th13_t, 0); //GLF_THETA_13
  glbSetOscParams(test_values, dm_atm_t, 3); //GLF_DELTA_ATM
  glbSetOscParams(test_values, th12_t, 4); //GLF_THETA_12
  glbSetOscParams(test_values, dm_sol_t, 5); //GLF_DELTA_SOLAR
  glbSetOscParams(test_values, real_t, 1); //GLF_REAL
  glbSetOscParams(test_values, imag_t, 2); //GLF_IMAG

  
  double errors[4] = { 5.0, 2.0, 8.0, 2.0};
/* 
  Here is where I calculate the chi-squared following

  There will be a "delta" for each experiment; assume no correlated uncertainties

  "user_data" contains flags indicated which experiments are to be included in the fit

*/

  glbChiSys(test_values,GLB_ALL,GLB_ALL);

  double deltaDB[52] = {0.0};
  double deltaDC[26] = {0.0};
  double deltaRENO[25] = {0.0};
  double deltaKamland[17] = {0.0};

  int i, j, k;
  double numerator, numerator2, denominator;

  int Nnuisance=4; /* The number of nuisance parameters */

  int CheckZeros;
  /* A quantity that we use to look for zero entries after shifting energy scales */

  /* Various important quantities used below */
  double emin, emax; 
  /* Get emin/emax for given experiment; used for energy-scale variation */

  /************
  *  Kamland  *
  ************/

  /* Get emin and emax for KAMLAND */
  /*
  glbGetEminEmax(12, &emin, &emax);

  double DeltaE0 = (emax-emin)/17;
  double t0 = emin/DeltaE0;
  double *Kamland_Simu_aux;
  double Kamland_Simu[29] = {0.0};
  double Kamland_Background[29] = {0.0};
  
  int linit = 6;

  Kamland_Simu_aux = glbGetSignalFitRatePtr(12, 0);
  for (int l=0; l<17; l++){
    Kamland_Background[l+linit]=Kamland_BG[l];
    Kamland_Simu[l+linit]=Kamland_Simu_aux[l];
  };

//  glbShiftEnergyScale(0.01*x[1]*1.0, glbGetSignalFitRatePtr(12, 0),
//                      Kamland_Simu, 17, emin, emax);  
//    glbShiftEnergyScale(0.01*x[3]*1.0, Kamland_Background_aux,
//                      Kamland_Background, 17, emin, emax);





CheckZeros = 0;
double totaltest = 0;
for (i=0; i<17; i++){
  double a = x[0]*0.01;
  double b = x[1]*0.01;

  double del = b*(i*t0+0.5)+i;
  int kd = floor(del);
//  printf("%f\n", del-kd);


  numerator = (1+b)*(1+a)*( (Kamland_Simu[linit+kd+1]-Kamland_Simu[linit+kd])*(del-kd)+ Kamland_Simu[linit+kd]);

  a = x[2]*0.01;
  b = x[3]*0.01;
  double background_rate = (1+b)*(1+a)*( (Kamland_Background[linit+kd+1]-Kamland_Background[linit+kd])*(del-kd)+ Kamland_Background[linit+kd]);
 



  deltaKamland[i] = Kamland_SG[i]-background_rate - numerator;
}
*/
  /*************
  *  DAYA BAY  *
  **************/

  /*
    We introduce a couple of fudge factors for Daya Bay in order to reproduce their 
    determinations of theta_13 and Dm31(Dmee)
  */
  double DB12_fudge = 0.9933;
  double DB13_fudge = 0.9973;

/*
  WATCH CAREFULLY HERE! The four near detectors at Daya Bay are also used to determine
  the normalization for NEOS, so they are formatted differently from the four far
  detectors! Namely, while the far detectors are binned in 0.2 MeV bins from 2.08-8.68
  MeV, the near detectors are binned in 0.1 MeV bins from 1.78-8.68 MeV bins. The finer
  resolution is driven by the binning used at NEOS; just combine the appropriate bins to
  get the numbers of events useful at Daya Bay!

  I've also included one additional bin at low energies for the far detectors, which I'm
  not including in the calculation, because it gives a different result for the lowest bin
  that I am including. I'm not sure why this is the case...

  I also need the three-neutrino numbers for the near detectors; these are to be used with
  the NEOS analysis.
*/

    //
    //  EH3 AD1  *
    //
  double *EH1_AD1_U235 = glbGetSignalFitRatePtr(0, 0);
  double *EH1_AD1_U238 = glbGetSignalFitRatePtr(0, 1);
  double *EH1_AD1_Pu239 = glbGetSignalFitRatePtr(0, 2);
  double *EH1_AD1_Pu241 = glbGetSignalFitRatePtr(0, 3);

  double *EH1_AD1_U235_0 = glbGetRuleRatePtr(0, 0);
  double *EH1_AD1_U238_0 = glbGetRuleRatePtr(0, 1);
  double *EH1_AD1_Pu239_0 = glbGetRuleRatePtr(0, 2);
  double *EH1_AD1_Pu241_0 = glbGetRuleRatePtr(0, 3);

    //
    //  EH3 AD2  *
    //

  double *EH1_AD2_U235 = glbGetSignalFitRatePtr(1, 0);
  double *EH1_AD2_U238 = glbGetSignalFitRatePtr(1, 1);
  double *EH1_AD2_Pu239 = glbGetSignalFitRatePtr(1, 2);
  double *EH1_AD2_Pu241 = glbGetSignalFitRatePtr(1, 3);

  double *EH1_AD2_U235_0 = glbGetRuleRatePtr(1, 0);
  double *EH1_AD2_U238_0 = glbGetRuleRatePtr(1, 1);
  double *EH1_AD2_Pu239_0 = glbGetRuleRatePtr(1, 2);
  double *EH1_AD2_Pu241_0 = glbGetRuleRatePtr(1, 3);

    //
    //  EH3 AD3  *
    //

  double *EH2_AD3_U235 = glbGetSignalFitRatePtr(2, 0);
  double *EH2_AD3_U238 = glbGetSignalFitRatePtr(2, 1);
  double *EH2_AD3_Pu239 = glbGetSignalFitRatePtr(2, 2);
  double *EH2_AD3_Pu241 = glbGetSignalFitRatePtr(2, 3);
   
    //
    //  EH3 AD8  *
    //

  double *EH2_AD8_U235 = glbGetSignalFitRatePtr(3, 0);
  double *EH2_AD8_U238 = glbGetSignalFitRatePtr(3, 1);
  double *EH2_AD8_Pu239 = glbGetSignalFitRatePtr(3, 2);
  double *EH2_AD8_Pu241 = glbGetSignalFitRatePtr(3, 3);

    //
    //  EH3 AD5  *
    //

  double *EH3_AD4_U235 = glbGetSignalFitRatePtr(4, 0);
  double *EH3_AD4_U238 = glbGetSignalFitRatePtr(4, 1);
  double *EH3_AD4_Pu239 = glbGetSignalFitRatePtr(4, 2);
  double *EH3_AD4_Pu241 = glbGetSignalFitRatePtr(4, 3);
   
    //
    //  EH3 AD5  *
    //

  double *EH3_AD5_U235 = glbGetSignalFitRatePtr(5, 0);
  double *EH3_AD5_U238 = glbGetSignalFitRatePtr(5, 1);
  double *EH3_AD5_Pu239 = glbGetSignalFitRatePtr(5, 2);
  double *EH3_AD5_Pu241 = glbGetSignalFitRatePtr(5, 3);

    //
    //  EH3 AD6  *
    //

  double *EH3_AD6_U235 = glbGetSignalFitRatePtr(6, 0);
  double *EH3_AD6_U238 = glbGetSignalFitRatePtr(6, 1);
  double *EH3_AD6_Pu239 = glbGetSignalFitRatePtr(6, 2);
  double *EH3_AD6_Pu241 = glbGetSignalFitRatePtr(6, 3);
   
    //
    //  EH3 AD7  *
    //

  double *EH3_AD7_U235 = glbGetSignalFitRatePtr(7, 0);
  double *EH3_AD7_U238 = glbGetSignalFitRatePtr(7, 1);
  double *EH3_AD7_Pu239 = glbGetSignalFitRatePtr(7, 2);
  double *EH3_AD7_Pu241 = glbGetSignalFitRatePtr(7, 3);

  int index;


//  REMINDER: We need to take the appropriate combinations of the Daya Bay bins
//  because of the fine binning with which we have calculated their spectra!

//  As a reminder, the Daya Bay spectra are broken into 227 bins (!!), each of width
//  0.05 MeV. We restructure them to match the binning at Daya Bay, as follows:

//    Bins 0-11: Daya Bay Bin 0 (0.7 - 1.3 MeV prompt)
//    Bins 12 + 5*(i-1) + (0-4): Daya Bay Bin i ( = 1-24 )
//                                (1.3 + 0.25*(i-1) + (0-0.25) MeV prompt)
//    Bins 132-225: Daya Bay Bin 25 (7.3 - 12.0 MeV prompt)

//  Below, we calculate the total number of events at each AD -- we initialize these here


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


  //Now return to business as usual: the ratio between EH2/EH1 and EH3/EH1

  for (i=0; i<26; i++){
    denominator = TotalAD1[i] + TotalAD2[i];
    numerator = TotalAD3[i] + TotalAD8[i];
    numerator2 = TotalAD4[i] + TotalAD5[i] + TotalAD6[i] + TotalAD7[i];

    deltaDB[i] = DayaBayRatio21[i] - DB12_fudge*numerator/denominator ;
    deltaDB[i+26] = DayaBayRatio31[i] - DB13_fudge*numerator2/denominator ;
    
  }

  //
  //  DOUBLE CHOOZ  *
  //
  
  
  //  We introduce a fudge factor for Double Chooz in order to reproduce their 
  //  determinations of theta_13 and Dm31/Dmee
  
  double DC_fudge = 1.0026;

    //
    //  NEAR DETECTOR  *
    //

  double *DC_ND_U235 = glbGetSignalFitRatePtr(8, 0);
  double *DC_ND_U238 = glbGetSignalFitRatePtr(8, 1);
  double *DC_ND_Pu239 = glbGetSignalFitRatePtr(8, 2);
  double *DC_ND_Pu241 = glbGetSignalFitRatePtr(8, 3);
   
    //
    //  FAR DETECTOR  *
    //

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

    deltaDC[i] = DC_Ratio[i] - DC_fudge*TotFD/TotND;
  }
  
  //
  //  RENO  *
  //

  
  //  We introduce a fudge factor for RENO in order to reproduce their determinations
  //  of theta_13 and Dm31/Dmee
  
  double RENO_fudge = 1.00884;

    //
    //  NEAR DETECTOR  *
    //

  double *RENO_ND_U235 = glbGetSignalFitRatePtr(10, 0);
  double *RENO_ND_U238 = glbGetSignalFitRatePtr(10, 1);
  double *RENO_ND_Pu239 = glbGetSignalFitRatePtr(10, 2);
  double *RENO_ND_Pu241 = glbGetSignalFitRatePtr(10, 3);
   
    //
    //  FAR DETECTOR  *
    //

  double *RENO_FD_U235 = glbGetSignalFitRatePtr(11, 0);
  double *RENO_FD_U238 = glbGetSignalFitRatePtr(11, 1);
  double *RENO_FD_Pu239 = glbGetSignalFitRatePtr(11, 2);
  double *RENO_FD_Pu241 = glbGetSignalFitRatePtr(11, 3);

  for (i=0; i<22; i++){ // The first 22 bins are totally normal...  
    TotND = RENO_ND_U235[i]* RENO_ND_F[0];
    TotND += RENO_ND_U238[i]* RENO_ND_F[1];
    TotND += RENO_ND_Pu239[i]* RENO_ND_F[2];
    TotND += RENO_ND_Pu241[i]* RENO_ND_F[3];

    TotFD = RENO_FD_U235[i]* RENO_FD_F[0];
    TotFD += RENO_FD_U238[i]* RENO_FD_F[1];
    TotFD += RENO_FD_Pu239[i]* RENO_FD_F[2];
    TotFD += RENO_FD_Pu241[i]* RENO_FD_F[3];

    deltaRENO[i] = RENO_Ratio[i] - RENO_fudge*TotFD/TotND;
  }

  for (i=0; i<2; i++){ // The next 2 bins are actually two bins combined...  
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

    deltaRENO[i+22] = RENO_Ratio[i+22] - RENO_fudge*TotFD/TotND;
  }

  TotND = 0.0;
  TotFD = 0.0;

  for (i=26; i<29; i++){ // The last bin is the sum of three bins  
    TotND += RENO_ND_U235[i]* RENO_ND_F[0];
    TotND += RENO_ND_U238[i]* RENO_ND_F[1];
    TotND += RENO_ND_Pu239[i]* RENO_ND_F[2];
    TotND += RENO_ND_Pu241[i]* RENO_ND_F[3];

    TotFD += RENO_FD_U235[i]* RENO_FD_F[0];
    TotFD += RENO_FD_U238[i]* RENO_FD_F[1];
    TotFD += RENO_FD_Pu239[i]* RENO_FD_F[2];
    TotFD += RENO_FD_Pu241[i]* RENO_FD_F[3];
  }
  deltaRENO[24] = RENO_Ratio[24] - RENO_fudge*TotFD/TotND;


  /********************************
  *  Putting the pieces together  *
  *********************************/

  /* Adding the contributions from each experiment */

double chi2 = 0.0;

// Kamland
/*
for (i=0; i<17; i++){
    chi2 += deltaKamland[i]*deltaKamland[i]/Kamland_SG[i];// *KamlandCovariance[i][j];
}
*/
// Daya Bay piece (EH2/EH1 and EH3/EH1) 
for (i=0; i<52; i++){
  for (j=0; j<52; j++){
    chi2 += deltaDB[i]*deltaDB[j]*DayaBayCovariance[i][j];
  }
}

// Double Chooz
for (i=0; i<26; i++){
  for (j=0; j<26; j++){
    chi2 += deltaDC[i]* deltaDC[j]*DC_Covariance[i][j];
  }
}
//RENO 
for (i=0; i<25; i++){
  for (j=0; j<25; j++){
    chi2 += deltaRENO[i]*deltaRENO[j]*RENO_Covariance[i][j];
  }
}


/*
  Adding nuisance parameters; helps insure convergence of minimization to include
  nuisance parameters for experiments that are not included in fit...
*/
/*
  for (i=0; i<Nnuisance; i++){
    //if(fabs(x[i]) > errors[i])
      chi2 += square( (fabs(x[i]) ) / errors[i]);
    //if(fabs(x[i]) < errors[i]) printf("%f\t",x[i]);
  }
//  printf("\n");  
*/
  /* Save the systematics parameters as starting values for the next step */
  //printf("%f\n",x[0]);
  //for (i=0; i < Nnuisance; i++)
  //  glf_systematic[i] = x[i];

  //return chi2;
  return -0.5*chi2;
}  
