#include "plover/pvt.h"


static void rot_small (const double x, double * result);
void rot_small (const double x, double * result)
{
    result[3 * 0] = 1.0;
    result[3 * 0 + 1] = x;
    result[3 * 0 + 2] = 0.0;
    result[3 * 1] = -x;
    result[3 * 1 + 1] = 1.0;
    result[3 * 1 + 2] = 0.0;
    result[3 * 2] = 0.0;
    result[3 * 2 + 1] = 0.0;
    result[3 * 2 + 2] = 1.0;
}
double pvt (double * rx_state, const u8 n_used, const navigation_measurement_t * * nav_meas, double * omp, double * H)
{
    double los [n_used * 3];
    
    for (s32 idx = 0; idx < n_used; idx++) {
        s32 j = idx;
        double loc [3];
        double tau;
        double arg [3];
        
        for (s32 idx2 = 0; idx2 < 3; idx2++) {
            arg[idx2] = rx_state[idx2] - nav_meas[j]->sat_pos[idx2];
        }
        tau = norm(3, arg) / GPS_C;
        
        double xk_new [3];
        double loc2 [3 * 3];
        
        rot_small(GPS_OMEGAE_DOT * tau, loc2);
        for (s32 idx2 = 0; idx2 < 3; idx2++) {
            double sum = 0;
            
            for (s32 idx3 = 0; idx3 < 3; idx3++) {
                sum += loc2[3 * idx2 + idx3] * nav_meas[j]->sat_pos[idx3];
            }
            xk_new[idx2] = sum;
        }
        for (s32 idx2 = 0; idx2 < 3; idx2++) {
            loc[idx2] = xk_new[idx2] - rx_state[idx2];
        }
        for (s32 idx2 = 0; idx2 < 3; idx2++) {
            los[3 * idx + idx2] = loc[idx2];
        }
    }
    
    double G [n_used * 4];
    
    for (s32 idx = 0; idx < n_used; idx++) {
        s32 j = idx;
        double loc [4];
        double loc2 [3];
        double arg [3];
        
        for (s32 idx2 = 0; idx2 < 3; idx2++) {
            arg[idx2] = -los[3 * j + idx2];
        }
        normalize(3, arg, loc2);
        
        s32 loc3 [1];
        
        loc3[0] = 1;
        for (s32 idx2 = 0; idx2 < 3; idx2++) {
            loc[0 + idx2] = loc2[idx2];
        }
        for (s32 idx2 = 0; idx2 < 1; idx2++) {
            loc[3 + idx2] = loc3[idx2];
        }
        for (s32 idx2 = 0; idx2 < 4; idx2++) {
            G[4 * idx + idx2] = loc[idx2];
        }
    }
    for (s32 idx = 0; idx < n_used; idx++) {
        s32 i = idx;
        
        omp[idx] = nav_meas[i]->pseudorange - norm(3, &los[3 * i]);
    }
    
    double loc [(3 + 1) * (3 + 1)];
    
    for (s32 idx = 0; idx < 3 + 1; idx++) {
        for (s32 idx2 = 0; idx2 < 3 + 1; idx2++) {
            double sum = 0;
            
            for (s32 idx3 = 0; idx3 < n_used; idx3++) {
                sum += G[(3 + 1) * idx3 + idx] * G[(3 + 1) * idx3 + idx2];
            }
            loc[(3 + 1) * idx + idx2] = sum;
        }
    }
    matrix_inverse(3 + 1, loc, H);
    
    double X [4 * n_used];
    
    for (s32 idx = 0; idx < 4; idx++) {
        for (s32 idx2 = 0; idx2 < n_used; idx2++) {
            double sum = 0;
            
            for (s32 idx3 = 0; idx3 < 4; idx3++) {
                sum += H[4 * idx + idx3] * G[(3 + 1) * idx2 + idx3];
            }
            X[n_used * idx + idx2] = sum;
        }
    }
    
    double correction [4];
    
    for (s32 idx = 0; idx < 4; idx++) {
        double sum = 0;
        
        for (s32 idx2 = 0; idx2 < n_used; idx2++) {
            sum += X[n_used * idx + idx2] * omp[idx2];
        }
        correction[idx] = sum;
    }
    
    double correction_norm;
    double arg [3];
    
    for (s32 idx = 0; idx < 3; idx++) {
        arg[idx] = correction[idx];
    }
    correction_norm = norm(3, arg);
    rx_state[3] = 0;
    for (s32 idx = 0; idx < 4; idx++) {
        rx_state[idx] = correction[idx] + rx_state[idx];
    }
    if (1.0e-3 < correction_norm) {
        return -correction_norm;
    }
    
    double tempvX [n_used];
    
    for (s32 idx = 0; idx < n_used; idx++) {
        s32 j = idx;
        double loc2;
        double pdot_pred;
        double sum = 0;
        
        for (s32 idx2 = 0; idx2 < 3; idx2++) {
            sum += G[(3 + 1) * j + idx2] * nav_meas[j]->sat_vel[idx2];
        }
        pdot_pred = -sum;
        loc2 = nav_meas[j]->doppler * GPS_C / GPS_L1_HZ - pdot_pred;
        tempvX[idx] = loc2;
    }
    for (s32 idx = 0; idx < 4; idx++) {
        double sum = 0;
        
        for (s32 idx2 = 0; idx2 < n_used; idx2++) {
            sum += X[n_used * idx + idx2] * tempvX[idx2];
        }
        rx_state[4 + idx] = sum;
    }
    return correction_norm;
}

