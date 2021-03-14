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

#ifndef GLF_SPECTRUM_AUX_2_H
#define GLF_SPECTRUM_AUX_2_H

/****************************************
*  INPUTS FOR THE SPECTRAL MEASUREMENTS *
*****************************************/

  /*************
  *  DAYA BAY  *
  **************/

/*
  The ratios of estimated number of events, N2/N1 and N3/N1. The binning is nonuniform:
	- The middle 24 bins are all 0.25 MeV from 1.3 to 7.3 MeV (prompt energy).
	- The first bin is 0.7-1.3 MeV.
	- The last bin is 7.3-12.0 MeV.
*/

static const double DayaBayRatio21[26] = {0.926932, 0.928954, 0.918739, 0.939571, 0.93593, 0.927909, 0.938499, 0.934976, 0.938284, 0.939889, 0.929927, 0.939565, 0.939061, 0.936539, 0.936932, 0.940709, 0.938614, 0.946836, 0.93802, 0.94908, 0.927143, 0.949198, 0.937487, 0.950238, 0.963763, 0.955083};

static const double DayaBayRatio31[26] = {0.288691, 0.287703, 0.276019, 0.277756, 0.275704, 0.273235, 0.271906, 0.271127, 0.274559, 0.272144, 0.271828, 0.275376, 0.27577, 0.276768, 0.275629, 0.280527, 0.279364, 0.281164, 0.283437, 0.281031, 0.268344, 0.27921, 0.276389, 0.293173, 0.271901, 0.273236};

/*
  These are arrays containing the fuel fractions for each Daya Bay AD, ordered:
  { U235, U238, Pu239, Pu241}

  The number here is the AD number (8 is in EH2!)

  These are from Table 9 of arXiv:1607.05378
*/

static const double DayaBayF1[4] = {0.5678, 0.0761, 0.30075, 0.05545};
static const double DayaBayF2[4] = {0.56605, 0.07615, 0.3021, 0.0556};
static const double DayaBayF3[4] = {0.5618, 0.0761, 0.30665, 0.0553};
static const double DayaBayF8[4] = {0.56345, 0.076, 0.3052, 0.05555};

static const double DayaBayF4[4] = {0.559, 0.076, 0.310, 0.055};
static const double DayaBayF5[4] = {0.559, 0.076, 0.310, 0.055};
static const double DayaBayF6[4] = {0.559, 0.076, 0.310, 0.055};
static const double DayaBayF7[4] = {0.552, 0.076, 0.315, 0.057};

  /*****************
  *  DOUBLE CHOOZ  *
  ******************/

static const double DC_Ratio[26] = {0.425358, 0.422108, 0.430931, 0.435575, 0.438825, 0.437896, 0.445791, 0.44254, 0.441147, 0.444397, 0.444862, 0.434181, 0.443933, 0.450899, 0.455542, 0.463436, 0.455078, 0.470402, 0.45322, 0.458793, 0.451363, 0.44997, 0.465294, 0.474581, 0.469938, 0.520089};

static const double DCF[4] = {0.511, 0.087, 0.340, 0.062};


  /*********
  *  RENO  *
  **********/

static const double RENO_Ratio[25] = {0.12592, 0.120914, 0.11629, 0.12038, 0.1154, 0.117704, 0.115853, 0.116933, 0.117296, 0.115695, 0.118126, 0.118856, 0.120047, 0.12019, 0.117643, 0.118033, 0.118237, 0.119066, 0.121666, 0.119507, 0.115003, 0.122545, 0.119143, 0.121214, 0.124873};

static const double RENO_ND_F[4] = {0.57264, 0.07309, 0.29911, 0.055161};
static const double RENO_FD_F[4] = {0.57447, 0.07340, 0.29742, 0.054752};

  /**************
  *   KAMLAND   *
  **************/

static const double Kamland_SG[17] = {156.7652499, 267.4114769, 292.5760472, 297.3054257, 216.3873188, 237.6310975, 245.8389359, 221.3605468, 233.8095146, 153.5731024, 99.97399211, 84.13111054, 57.15552292, 23.7697531, 12.94926282, 6.078511422, 1.705896572};
static const double Kamland_BG[17] = {117.7189997, 132.3789758, 92.90521391, 68.05917575, 31.88233391, 6.627548887, 1.676206456, 1.492317065, 5.673098239, 3.782065492, 0.654213836, 11.489414, 11.27458771, 1.420708301, 0.255498155, 0, 0};

static const double Kamland_ERROR[17] = {210, 284, 284, 284, 284, 284, 284, 284, 284, 200, 148, 148, 123, 50, 30, 30, 10};

static const double Kamland_F[4] = {0.567, 0.078, 0.298, 0.057};

#endif /* !GLF_SPECTRUM_AUX_2_H */
