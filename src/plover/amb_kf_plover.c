#include <stdio.h>
#include <constants.h>
#include <math.h>
#include <amb_kf.h>
#include <plover/amb_kf_plover.h>
double simple_amb_measurement(double carrier, double code)
{
  return (carrier + (code / GPS_L1_LAMBDA_NO_VAC));
}
void incorporate_obs2(nkf_t * kf, double decor_obs[kf->obs_dim])
{
  for (int i = 0; i < kf->obs_dim; i++) {
    double k[kf->state_dim];
    double var3[kf->state_dim];
    for (int var4 = 0; var4 < kf->state_dim; var4++) {
      var3[var4] = kf->state_cov_D[var4];
    }
    incorporate_scalar_measurement(kf->state_dim, var0, var1, var2, var3, k);
  }
  return 0;
}
