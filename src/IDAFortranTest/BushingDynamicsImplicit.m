function [res, drdz, drdzp] = BushingDynamicsImplicit(t, z, zp, ...
                                m, II, kv, kr, cv, cr, h, g, thet0, ...
                                ux, uz, dtsample)
% Bushing dynamics including coupling between rocking and vertical
% INPUTS
% t = time
% z = state (Delta, theta, \dot{Delta}, \dot{theta})
% zp = zprime
% m = bushing mass
% II = bushing moment of inertia about pivot point
% kv = vertial stiffness of cover plate
% kr = rotational stiffness of cover plate
% cv = damping coefficient for vertical motion
% cr = damping coefficient for rotational motion
% h = height of bushing center of mass from pivot point
% g = acceleration due to gravity
% thet0 = misalignment from vertical (imperfection) in radians
% ux, uz = input ground acceleration in x and z directions (times series
%          vectors)
% dtsample = sampling time of input ground accelerations

Delt = z(1);
thet = z(2);
dDelt = z(3); % \dot{Delta}
dthet = z(4); % \dot{theta}
ddDelt = zp(3);
ddthet = zp(4);

c = cos(thet+thet0);
s = sin(thet+thet0);

% mass matrix
M = [m -m*h*s; -m*h*s II];

% input ground accelerations
n = floor(t/dtsample) + 1;
if (n > length(ux)-1)
    ux_ = 0;
    uz_ = 0;
else
    ux_ = ux(n) + (ux(n+1)-ux(n))/dtsample*(t - (n-1)*dtsample);
    uz_ = uz(n) + (uz(n+1)-uz(n))/dtsample*(t - (n-1)*dtsample);
end

% force vector
F = [-m*h*c*dthet^2+kv*Delt+m*g*(1+uz_)+cv*dDelt; ...
     kr*thet-m*g*(1+uz_)*h*s+m*g*ux_*h*c+cr*dthet];
 
res = zeros(4,1);
res(1:2) = zp(1:2) - z(3:4);
res(3:4) = M*zp(3:4) + F;

if nargout>1
    K = [kv m*h*(-c*ddthet+s*dthet^2); ...
         0  kr-m*h*(c*(g*(1+uz_)+ddDelt)+s*g*ux_)];
    C = [cv -2*m*h*c*dthet; 0 cr];
    drdz = [zeros(2) -eye(2); K C];
    drdzp = blkdiag(eye(2),M);
end
