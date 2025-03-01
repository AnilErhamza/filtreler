function [Quaternion] = madgwick_wom(Quaternion, Gyroscope, Accelerometer)

SamplePeriod = 1/256;
Beta = 0.1;
q = Quaternion;
 % short name local variable for readability
% Normalise accelerometer measurement
if(norm(Accelerometer) == 0), return; end	% handle NaN
Accelerometer = Accelerometer / norm(Accelerometer);	% normalise magnitude

% Gradient decent algorithm corrective step
F = [2*(q(2)*q(4) - q(1)*q(3)) - Accelerometer(1)
     2*(q(1)*q(2) + q(3)*q(4)) - Accelerometer(2)
     2*(0.5 - q(2)^2 - q(3)^2) - Accelerometer(3)];
J = [-2*q(3),	2*q(4),    -2*q(1),	2*q(2)
      2*q(2),     2*q(1),     2*q(4),	2*q(3)
      0,         -4*q(2),    -4*q(3),	0    ];
step = (J'*F);
step = step / norm(step);	% normalise step magnitude

% Compute rate of change of quaternion
qDot = 0.5 * quaternProd(q, [0 Gyroscope(1) Gyroscope(2) Gyroscope(3)]) - Beta * step';

% Integrate to yield quaternion
q = q + qDot * SamplePeriod;
Quaternion = q / norm(q); % normalise quaternion




end

