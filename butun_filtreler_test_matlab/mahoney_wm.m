function Quaternion = mahoney_wm(Quaternion, Gyroscope, Accelerometer, Magnetometer)
q = Quaternion; % short name local variable for readability
SamplePeriod = 1/256;
Kp = 1;
Ki = 0.1; 
eInt = [0 0 0];
% Normalise accelerometer measurement
if(norm(Accelerometer) == 0), return; end   % handle NaN
Accelerometer = Accelerometer / norm(Accelerometer);    % normalise magnitude
 
% Normalise magnetometer measurement
if(norm(Magnetometer) == 0), return; end    % handle NaN
Magnetometer = Magnetometer / norm(Magnetometer);   % normalise magnitude
 
% Reference direction of Earth's magnetic feild
h = quaternProd(q, quaternProd([0 Magnetometer], quaternConj(q)));
b = [0 norm([h(2) h(3)]) 0 h(4)];
            
% Estimated direction of gravity and magnetic field
v = [2*(q(2)*q(4) - q(1)*q(3))
     2*(q(1)*q(2) + q(3)*q(4))
     q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2];
w = [2*b(2)*(0.5 - q(3)^2 - q(4)^2) + 2*b(4)*(q(2)*q(4) - q(1)*q(3))
     2*b(2)*(q(2)*q(3) - q(1)*q(4)) + 2*b(4)*(q(1)*q(2) + q(3)*q(4))
     2*b(2)*(q(1)*q(3) + q(2)*q(4)) + 2*b(4)*(0.5 - q(2)^2 - q(3)^2)]; 
 
% Error is sum of cross product between estimated direction and measured direction of fields
e = cross(Accelerometer, v) + cross(Magnetometer, w);
if(Ki > 0)
    eInt = eInt + e * SamplePeriod;   
else
    eInt = [0 0 0];
end
            
% Apply feedback terms
Gyroscope = Gyroscope + Kp * e + Ki * eInt;            
            
% Compute rate of change of quaternion
qDot = 0.5 * quaternProd(q, [0 Gyroscope(1) Gyroscope(2) Gyroscope(3)]);
 
% Integrate to yield quaternion
q = q + qDot * SamplePeriod;
Quaternion = q / norm(q); % normalise quaternion
 end