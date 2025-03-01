function [m1 m2 m3] = Kalman_Filter_before_magneto(z)

persistent H Q R
persistent x P
persistent firstRun


if isempty(firstRun)

  
  H = eye(3);

  Q = 0.1*eye(3); 
  R = 200*eye(3);

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


m1  =x(1); 
m2=x(2); 
m3  =x(3); 