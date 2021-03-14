#include<glf_probability.h>
#include<glf_precomputed_probabilities.h>
#include<glf_pre.h>
#include<glf_types.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define GLF_LE_KM_GEV 1.0/(GLB_EV_TO_KM_FACTOR*4.0*1E9)

double bugey_15_cos[401][2];
double bugey_15_sin[401][2];

double glf_compute_cos(double (*func)(double), double *vect ,void *user_data)
{
  glf_distance_data *distances = (glf_distance_data *) user_data;

  size_t dlength= distances->length;
  double (*data)[2]= distances->data;

};