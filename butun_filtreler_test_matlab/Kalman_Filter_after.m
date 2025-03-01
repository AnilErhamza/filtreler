function [psi theta phi] = Kalman_Filter_after(z)

persistent H Q R
persistent x P
persistent firstRun


if isempty(firstRun)

  H = eye(3);

  Q = 0.001*eye(3); %process noise covariance matrix
  R = 1*eye(3); %measurement noise covariance matrix

  x = [0 0 0]';  
  P = 1*eye(3);
 
  firstRun = 1;  
end
 A = eye(3);

% Kalman Filter algorithm
xp = A*x;
Pp = A*P*A' + Q;

K = Pp*H'*inv(H*Pp*H' + R);

x = xp + K*(z - H*xp);   
P = Pp - K*H*Pp;

 % Euler param. to yaw-pitch-roll
psi  =x(1); % yaw [rad]
theta=x(2); % pitch [rad]
phi  =x(3); % roll [rad]