#include <stdio.h>
#include "madgwick_lib.h"
#include "madgwick_wom.h"

int main(){
    struct quaternion q_test;
    float accel[3]={0.01904297, -0.05224609, 0.9780273};
    float gyro[3]={-0.9375,	-1.25,	0.875};

    struct euler euler;
    q_test= madgwick_wom(accel[0],accel[1],accel[2],gyro[0]* (PI/180),gyro[1]* (PI/180),gyro[2]* (PI/180));

    printf("%f,",q_test.q1);
    printf("%f,",q_test.q2);
    printf("%f,",q_test.q3);
    printf("%f",q_test.q4);

    euler = eulerAngles(q_test);

    printf("%f\n",euler.roll);
    printf("%f",euler.pitch);
    printf("%f",euler.yaw);


    return 0;


}


