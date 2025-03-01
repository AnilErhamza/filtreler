#include "madgwick_lib.h"

struct quaternion q_est = { 1, 0, 0, 0};      


struct quaternion madgwick_wom(float ax, float ay, float az, float gx, float gy, float gz){
    
    struct quaternion q_est_prev = q_est;
    struct quaternion q_est_dot = {0};            
    
    struct quaternion q_a = {0, ax, ay, az};    
    
    float F_g [3] = {0};                        
    float J_g [3][4] = {0};                     
    
    struct quaternion gradient = {0};
    
    struct quaternion q_w;                  
    q_w.q1 = 0;                              
    q_w.q2 = gx;
    q_w.q3 = gy;
    q_w.q4 = gz;
    
    quat_scalar(&q_w, 0.5);                  
    q_w = quat_mult(q_est_prev, q_w);        

    
    quat_Normalization(&q_a);              
    
    F_g[0] = 2*(q_est_prev.q2 * q_est_prev.q4 - q_est_prev.q1 * q_est_prev.q3) - q_a.q2;
    F_g[1] = 2*(q_est_prev.q1 * q_est_prev.q2 + q_est_prev.q3* q_est_prev.q4) - q_a.q3;
    F_g[2] = 2*(0.5 - q_est_prev.q2 * q_est_prev.q2 - q_est_prev.q3 * q_est_prev.q3) - q_a.q4;
    
    
    J_g[0][0] = -2 * q_est_prev.q3;
    J_g[0][1] =  2 * q_est_prev.q4;
    J_g[0][2] = -2 * q_est_prev.q1;
    J_g[0][3] =  2 * q_est_prev.q2;
    
    J_g[1][0] = 2 * q_est_prev.q2;
    J_g[1][1] = 2 * q_est_prev.q1;
    J_g[1][2] = 2 * q_est_prev.q4;
    J_g[1][3] = 2 * q_est_prev.q3;
    
    J_g[2][0] = 0;
    J_g[2][1] = -4 * q_est_prev.q2;
    J_g[2][2] = -4 * q_est_prev.q3;
    J_g[2][3] = 0;
    
    
    gradient.q1 = J_g[0][0] * F_g[0] + J_g[1][0] * F_g[1] + J_g[2][0] * F_g[2];
    gradient.q2 = J_g[0][1] * F_g[0] + J_g[1][1] * F_g[1] + J_g[2][1] * F_g[2];
    gradient.q3 = J_g[0][2] * F_g[0] + J_g[1][2] * F_g[1] + J_g[2][2] * F_g[2];
    gradient.q4 = J_g[0][3] * F_g[0] + J_g[1][3] * F_g[1] + J_g[2][3] * F_g[2];
    
    
    quat_Normalization(&gradient);
  
    
    quat_scalar(&gradient, BETA);             
    quat_sub(&q_est_dot, q_w, gradient);        
    quat_scalar(&q_est_dot, DELTA_T);
    quat_add(&q_est, q_est_prev, q_est_dot);    
    quat_Normalization(&q_est);                 

    struct quaternion q_est_out = q_est;

    return q_est_out;                                            
   
}


