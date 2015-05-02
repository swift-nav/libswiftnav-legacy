#ifndef PLOVER_GENERATED_AMB_KF_PLOVER_H
#define PLOVER_GENERATED_AMB_KF_PLOVER_H

#include <common.h>
#include <amb_kf.h>
double simple_amb_measurement(double carrier, double code);
void incorporate_obs2(nkf_t * kf, double decor_obs[kf->obs_dim]);


#endif /* PLOVER_GENERATED_AMB_KF_PLOVER_H */
