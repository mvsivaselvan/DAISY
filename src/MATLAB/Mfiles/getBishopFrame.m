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

I3 = [1 0 0; 0 1 0; 0 0 0];

opt = odeset('RelTol',2.5e-14,'AbsTol',1e-16);
[~,Q] = ode45(@Bishop,xg,q0,opt);

ng = length(xg);
R = zeros(3,3*ng);

for n = 1:ng
    qn = Q(n,:)';
    R(:,3*(n-1)+1:3*(n-1)+3) = quat2rot(qn);
end

    function qdot = Bishop(x, q)
        if (coder.target('MATLAB')) % running in MATLAB
           colmat = spcol(knots, d+1, [x;x;x]);
        else % generating code
           colmat = spcolC(knots, length(knots), d+1, x, 1, 3);
        end
        % B = colmat(1,:)'; % not used
        Bp = colmat(2,:)';
        Bpp = colmat(3,:)';

        xip = P*Bp;
        xipp = P*Bpp;

        nxip = sqrt(xip'*xip);
        % if (nxip == 0) % This situation will never heppen, it is to avoid
        %                % nan and inf checks in code generation
        %     qdot = [0; 0; 0; 0];
        %     return
        % end

        tau = xip/nxip;

        PP = I3 - tau*tau';
        taup = (PP*xipp)/nxip;

        omega = mycross(tau, taup);
        Omega = [0 -omega'; omega hat(omega)];

        k = 0.1; % feedback gain to prevent drift of q from unit norm
        
        qTq = q'*q;
        qdot = 0.5*Omega*q + k/2/qTq*(1 - qTq)*q;        
    end

    function RR = quat2rot(qq)
        normqq = sqrt(qq'*qq);
        RR = 1/normqq^2 * ...
           [qq(1)^2+qq(2)^2-qq(3)^2-qq(4)^2 2*(qq(2)*qq(3)-qq(1)*qq(4)) ...
                                           2*(qq(2)*qq(4)+qq(1)*qq(3)); ...
            2*(qq(2)*qq(3)+qq(1)*qq(4)) qq(1)^2+qq(3)^2-qq(2)^2-qq(4)^2 ...
                                           2*(qq(3)*qq(4)-qq(1)*qq(2)); ...
            2*(qq(2)*qq(4)-qq(1)*qq(3)) 2*(qq(3)*qq(4)+qq(1)*qq(2)) ...
                                         qq(1)^2+qq(4)^2-qq(2)^2-qq(3)^2]; 
    end

end