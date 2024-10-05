function zdot = BushingRocking(t, z, ...
                               m, II, kr, cr, h, g, thet0, ...
                               ux, uz, dtsample)
% rocking dynamics without coupling with vertical
% INPUTS - same as in BushingDynamics

% input ground accelerations
n = floor(t/dtsample) + 1;
if (n > length(ux)-1)
    ux_ = 0;
    uz_ = 0;
else
    ux_ = ux(n) + (ux(n+1)-ux(n))/dtsample*(t - (n-1)*dtsample);
    uz_ = uz(n) + (uz(n+1)-uz(n))/dtsample*(t - (n-1)*dtsample);
end

zdot = zeros(2,1);
zdot(1) = z(2);
zdot(2) = (-kr*z(1)-cr*z(2)+m*g*(1+uz_)*h*sin(z(1)+thet0)...
           -m*g*h*ux_*cos(z(1)+thet0))/II;
