#include "madgwick_lib.h"

struct quaternion q_est = { 1, 0, 0, 0};      


struct quaternion madgwick_wm(float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz){
    
    struct quaternion q_est_prev = q_est;
    struct quaternion q_est_dot = {0};            
    
    struct quaternion q_a = {0, ax, ay, az};
    struct quaternion q_m = {0, mx, my, mz};    
    
    float F_g [6] = {0};                        
    float J_g [6][4] = {0};                     
    
    struct quaternion gradient = {0};
    
    struct quaternion q_w;                  
    q_w.q1 = 0;                              
    q_w.q2 = gx;
    q_w.q3 = gy;
    q_w.q4 = gz;
    
    quat_scalar(&q_w, 0.5);                  
    q_w = quat_mult(q_est_prev, q_w);        

    quat_Normalization(&q_a); 
    quat_Normalization(&q_m);
    struct quaternion a = quat_conjugate(q_est_prev);
    struct quaternion c = quat_mult(q_m, a);
    struct quaternion h = quat_mult(q_est_prev,c);
    struct quaternion normalize = {0, h.q2, h.q3, 0};
    float normalized = quat_Norm(normalize);
    struct quaternion b = {0, normalized, 0, h.q4};        
    
    F_g[0] = 2*(q_est_prev.q2 * q_est_prev.q4 - q_est_prev.q1 * q_est_prev.q3) - q_a.q2;
    F_g[1] = 2*(q_est_prev.q1 * q_est_prev.q2 + q_est_prev.q3* q_est_prev.q4) - q_a.q3;
    F_g[2] = 2*(0.5 - q_est_prev.q2 * q_est_prev.q2 - q_est_prev.q3 * q_est_prev.q3) - q_a.q4;
    F_g[3] = 2*b.q2*(0.5 - q_est_prev.q3*q_est_prev.q3 - q_est_prev.q4*q_est_prev.q4) + 2*b.q4*(q_est_prev.q2*q_est_prev.q4 - q_est_prev.q1*q_est_prev.q3) - q_m.q2;
    F_g[4] = 2*b.q2*(q_est_prev.q2*q_est_prev.q3 - q_est_prev.q1*q_est_prev.q4) + 2*b.q4*(q_est_prev.q1*q_est_prev.q2 + q_est_prev.q3*q_est_prev.q4) - q_m.q3;
    F_g[5] = 2*b.q2*(q_est_prev.q1*q_est_prev.q3 + q_est_prev.q2*q_est_prev.q4) + 2*b.q4*(0.5 - q_est_prev.q2*q_est_prev.q2 - q_est_prev.q3*q_est_prev.q3) - q_m.q4;
    
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

    J_g[3][0] = -2*b.q4*q_est_prev.q3;
    J_g[3][1] = 2*b.q4*q_est_prev.q4;
    J_g[3][2] = -4*b.q2*q_est_prev.q3 - 2*b.q4*q_est_prev.q1;
    J_g[3][3] = -4*b.q2*q_est_prev.q4 + 2*b.q4*q_est_prev.q2;

    J_g[4][0] = -2*b.q2*q_est_prev.q4 + 2*b.q4*q_est_prev.q2;
    J_g[4][1] = 2*b.q2*q_est_prev.q3 + 2*b.q4*q_est_prev.q1;
    J_g[4][2] = 2*b.q2*q_est_prev.q2 + 2*b.q4*q_est_prev.q4;
    J_g[4][3] = -2*b.q2*q_est_prev.q1 + 2*b.q4*q_est_prev.q3;

    J_g[5][0] = 2*b.q2*q_est_prev.q3;
    J_g[5][1] = 2*b.q2*q_est_prev.q4 - 4*b.q4*q_est_prev.q2;
    J_g[5][2] = 2*b.q2*q_est_prev.q1 - 4*b.q4*q_est_prev.q3;
    J_g[5][3] = 2*b.q2*q_est_prev.q2;
    
    
    gradient.q1 = J_g[0][0] * F_g[0] + J_g[1][0] * F_g[1] + J_g[2][0] * F_g[2] + J_g[3][0] * F_g[3] + J_g[4][0] * F_g[4] + J_g[5][0] * F_g[5];
    gradient.q2 = J_g[0][1] * F_g[0] + J_g[1][1] * F_g[1] + J_g[2][1] * F_g[2] + J_g[3][1] * F_g[3] + J_g[4][1] * F_g[4] + J_g[5][1] * F_g[5];
    gradient.q3 = J_g[0][2] * F_g[0] + J_g[1][2] * F_g[1] + J_g[2][2] * F_g[2] + J_g[3][2] * F_g[3] + J_g[4][2] * F_g[4] + J_g[5][2] * F_g[5];
    gradient.q4 = J_g[0][3] * F_g[0] + J_g[1][3] * F_g[1] + J_g[2][3] * F_g[2] + J_g[3][3] * F_g[3] + J_g[4][3] * F_g[4] + J_g[5][3] * F_g[5];
    
    
    quat_Normalization(&gradient);
  
    
    quat_scalar(&gradient, BETA);             
    quat_sub(&q_est_dot, q_w, gradient);        
    quat_scalar(&q_est_dot, DELTA_T);
    quat_add(&q_est, q_est_prev, q_est_dot);    
    quat_Normalization(&q_est);                 

    struct quaternion q_est_out = q_est;

    return q_est_out;                                            
   
}
