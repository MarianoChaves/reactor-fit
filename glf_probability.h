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

#ifndef GLF_PROBABILITY_H
#define GLF_PROBABILITY_H




#define GLF_THETA_13 0
#define GLF_REAL 1
#define GLF_IMAG 2
#define GLF_DELTA_ATM 3
#define GLF_THETA_12 4
#define GLF_DELTA_SOLAR 5


/**************************************************************************/

/* Setting and getting parameters */

int glf_set_oscillation_parameters(glb_params p, void *user_data);
int glf_get_oscillation_parameters(glb_params p, void *user_data);

/**************************************************************************/

/* Generic probability matrices */


int glf_nsi_probability_pre(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
				 double filter_sigma, void *user_data);

int glf_nsi_probability_kamland(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data);
/**************************************************************************/

#endif /* !GLF_PROBABILITY_H */
