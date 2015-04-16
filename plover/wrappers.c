#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "extern_defs.h"
#include "generated.h"

#include "constants.h"
#include "track.h"


void wrap_pvt(double rx_state[],
              const u8 n_used,
              const navigation_measurement_t nav_meas[n_used],
              double correction[4],
              double G[n_used][4],
              double X[4][n_used]) {

  double sat_pos[n_used][3];
  double pseudo[n_used];

  for (int i = 0; i < n_used; i++) {
    for (int x = 0; x < 3; x++) {
      sat_pos[i][x] = nav_meas[i].sat_pos[x];
    }
    pseudo[i] = nav_meas[i].pseudorange;
  }

  pvt((int)n_used, sat_pos, pseudo, rx_state, correction, G, X);
  //pvt((int)n_used, sat_pos, pseudo, (double(*)[3])rx_state, correction, G, X);
}
