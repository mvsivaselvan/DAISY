function zdot = BushingDynamics(t, z, ...
                                m, II, kv, kr, cv, cr, h, g, thet0, ...
                                ux, uz, dtsample)
% Bushing dynamics including coupling between rocking and vertical
% INPUTS
% t = time
% z = state (Delta, theta, \dot{Delta}, \dot{theta})
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

% mass matrix
M = [m -m*h*sin(thet+thet0); -m*h*sin(thet+thet0) II];

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
F = [-m*h*cos(thet+thet0)*dthet^2+kv*Delt+m*(g+uz_)+cv*dDelt; ...
     kr*thet-m*(g+uz_)*h*sin(thet+thet0)+m*ux_*h*cos(thet+thet0)+cr*dthet];
 
zdot = zeros(4,1);
zdot(1:2) = z(3:4);
zdot(3:4) = -M\F;
