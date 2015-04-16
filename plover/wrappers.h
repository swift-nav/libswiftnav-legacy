void wrap_pvt(double rx_state[],
              const u8 n_used,
              const navigation_measurement_t nav_meas[n_used],
              double correction[4],
              double G[n_used][4],
              double X[4][n_used]);
