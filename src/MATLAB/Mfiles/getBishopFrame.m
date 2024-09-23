function R = getBishopFrame(P, knots, d, xg)

% Computes the (twist-free) Bishop framing of a given curve
% INPUTS
% P = control points of curve arranged as 3XN matrix
% knots = knot vector for B-spline basis functions
% d = degree of B-spline basis functions used
% xg = collocation points (in arclength parametrization) where Bishop
%      frames are to be computed
% q0 = initial condition for integration (quaternion representation 
%      of R at s=0)
% OUTPUT
% R = cell array with length(xg) rotation matrices defining the Bishop
%     frame

tau0 = (P(:,2) - P(:,1));
normtau0 = sqrt(tau0'*tau0);
tau0 = tau0/normtau0;
theta = acos(tau0(1)); % tau0(1) = dot product of tau0 with [1;0;0]
if (abs(theta) < 1e-15) % tau0 is oriented along [1;0;0]
    q0 = [1; 0; 0; 0];
else
    rotaxis = [0; -tau0(3); tau0(2)]/sqrt(tau0(2)^2+tau0(3)^2);
    % i.e. unit vector along cross([1;0;0],tau0)
    q0 = [cos(theta/2); sin(theta/2)*rotaxis];
end

opt = odeset('RelTol',2.5e-14,'AbsTol',1e-16);
[~,Q] = ode45(@(x,q)Bishop(x,q,P,knots,d),[0; xg; 1],q0,opt);

ng = length(xg) + 2; % + 2 for the two end points 0 and 1
R = zeros(3,3*ng);

for n = 1:ng
    qn = Q(n,:)';
    R(:,3*(n-1)+1:3*(n-1)+3) = quat2rot(qn);
end
