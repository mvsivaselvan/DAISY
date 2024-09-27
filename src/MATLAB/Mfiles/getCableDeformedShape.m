function P = getCableDeformedShape(d1,phi1,gamm1,d2,phi2,gamm2,Pmid, ...
                                   x01, RJ1, RE1, r1, Rb10, ...
                                   x02, RJ2, RE2, r2, Rb20)

% INPUTS 
% d1 = displacement of joint1 
% phi1 = exponential coordinates of rotation of joint
% gamm1 = distance between first and second control points
% d2, phi2, gamm2 = same of joint 2 (end)
% Pmid = displaced control points 3:N-2 (3*(N-4) matrix)
% x01 = reference position of joint 1 (3x1) vector
% RJ1 = rotation of joint 1 coordinate frame with respect to global
% RE1 = rotation of cable end with respect to joint 1 coordinate frame
% r1 = position of cable end relative to joint 1 in joint coordinate system
%     (end offset, 3x1 vector)
% x02, RJ2, RE2, r2 = same for joint 2
% OUTPUT
% P = displaced control points

qbar1 = [d1; phi1; gamm1];
qbar2 = [d2; phi2; gamm2];

eta = zeros(3,1);
rho1 = zeros(3,1);
rho2 = zeros(3,1);

q1 = CableBCTransinCoord(qbar1, ...
                         x01, RJ1, RE1, r1, Rb10, ...
                         eta, rho1,rho2);
q2 = CableBCTransinCoord(qbar2, ...
                         x02, RJ2, RE2, r2, Rb20, ...
                         eta, rho1,rho2);
P = [q1(1:3) q1(4:6) Pmid q2(4:6) q2(1:3)];
