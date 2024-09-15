function R = getBishopFrame(P, knots, d, xg, q0)

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

if (nargin < 5) % q0 not provided
    tau0 = (P(:,2) - P(:,1));
    tau0 = tau0/norm(tau0);
    theta = acos(tau0(1)); % tau0(1) = dot product of tau0 with [1;0;0]
    if (abs(theta) < 1e-15) % tau0 is oriented along [1;0;0]
        q0 = [1; 0; 0; 0];
    else
        rotaxis = [0; -tau0(3); tau0(2)]/sqrt(tau0(2)^2+tau0(3)^2); 
                    % i.e. unit vector along cross([1;0;0],tau0)
        q0 = [cos(theta/2); sin(theta/2)*rotaxis];
    end
end
    

k = 0.1; % feedback gain to prevent drift of q from unit norm

opt = odeset('RelTol',2.5e-14,'AbsTol',1e-16);
[~,Q] = ode45(@Bishop,xg,q0,opt);

ng = length(xg);
R = zeros(3,3*ng);

for n = 1:ng
    qn = Q(n,:)';
    R(:,3*(n-1)+1:3*(n-1)+3) = quat2rot(qn);
end

    function qdot = Bishop(x, q)
        colmat = spcol(knots, d+1, [x; x; x]);
        % B = colmat(1,:)'; % not used
        Bp = colmat(2,:)';
        Bpp = colmat(3,:)';
        
        xip = P*Bp;
        xipp = P*Bpp;
        
        nxip = norm(xip);
        tau = xip/nxip;
        
        PP = eye(3) - tau*tau';
        taup = (PP*xipp)/nxip;
        
        omega = cross(tau, taup);
        Omega = [0 -omega'; omega hat(omega)];
        
        qdot = 0.5*Omega*q + k/2/norm(q)^2*(1 - q'*q)*q;        
    end

    function RR = quat2rot(qq)
        RR = 1/norm(qq)^2 * ...
           [qq(1)^2+qq(2)^2-qq(3)^2-qq(4)^2 2*(qq(2)*qq(3)-qq(1)*qq(4)) ...
                                           2*(qq(2)*qq(4)+qq(1)*qq(3)); ...
            2*(qq(2)*qq(3)+qq(1)*qq(4)) qq(1)^2+qq(3)^2-qq(2)^2-qq(4)^2 ...
                                           2*(qq(3)*qq(4)-qq(1)*qq(2)); ...
            2*(qq(2)*qq(4)-qq(1)*qq(3)) 2*(qq(3)*qq(4)+qq(1)*qq(2)) ...
                                         qq(1)^2+qq(4)^2-qq(2)^2-qq(3)^2]; 
    end

end