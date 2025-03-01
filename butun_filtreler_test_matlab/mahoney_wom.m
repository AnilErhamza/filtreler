function Quaternion = mahoney_wom(Quaternion, Gyroscope, Accelerometer)
q = Quaternion; 
SamplePeriod = 1/256;
Kp = 1;
Ki = 0; 
eInt = [0 0 0];

% Normalise accelerometer measurement
if(norm(Accelerometer) == 0), return; end   % handle NaN
Accelerometer = Accelerometer / norm(Accelerometer);	% normalise magnitude
 
% Estimated direction of gravity and magnetic flux
v = [2*(q(2)*q(4) - q(1)*q(3))
     2*(q(1)*q(2) + q(3)*q(4))
    q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2];
 
% Error is sum of cross product between estimated direction and measured direction of field
e = cross(Accelerometer, v); 
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