function [F, K, C, M, B] = RigidBodyForce(d, phi, dd, phid, ddd, phidd, ...
                                       m, II, KT, CT, KR, CR, rr, ...
                                       RJ, u)

% Inputs
% d = displacement of joint in joint coordinate system
% phi = exponential coordinates of rotation in joint coordinate system
% dd = \dot{d}
% phid = \dot{phi}
% ddd = \ddot(d}
% phidd = \ddot(phi}
% m = mass
% II = moment of inertia about joint
% KT = translational stiffness matrix (3x3)
% CT = translational damping matrix (3x3)
% KR = rotational stiffness matrix (3x3)
% CR = rotational damping matrix (3x3)
% rr = position vector of center of mass relative to joint in joint
%      coordinate system
% RJ = oritentation of joint coordinate system relative to global
% u = force per unit mass (acceleration due to gravity and ground) 
% Outputs
% F = joint force (6x1 - force and moment)
% K = stiffness matrix (derivative of force w.r.t. (d,phi))
% C = derivative of force w.r.t. (ddot,phidot)
% M = mass matrix (derivative of force w.r.t. (dddot,phiddot))
% B = derivative of force w.r.t. input

if nargout == 1
    a = rothelper(phi, 2);
else
    a = rothelper(phi, 3);
end

I3 = eye(3);
hphi = hat(phi);
hphi2 = hphi^2;
R = I3 + a(1)*hphi + a(2)*hphi2;
dexp_ = I3 - a(2)*hphi + a(9)*hphi2;

wtJoint = m*(RJ'*u); % Weight in the joint frame
wtBody = R'*wtJoint; % weight in the body frame

w = dexp_*phid;
D1w_phid = D1w(phid);
wd = dexp_*phidd + D1w_phid*phid;

wxr = mycross(w,rr);
abar = mycross(wd,rr) + mycross(w,wxr);
mubar =m*mycross(rr,R'*ddd) + II*wd + mycross(w,II*w) - mycross(rr,wtBody);
wmubar = dexp_*mubar;

F = zeros(6,1);
F(1:3) = m*ddd + m*R*abar + KT*d + CT*dd - wtJoint;
F(4:6) = R*wmubar  + KR*phi + CR*phid;

if nargout > 1 % need derivatives K, C and M
    rhat = hat(rr);
    what = hat(w);
    
    tmp1 = D1sqw(phid,phid) + D1w(phidd);
    tmp2 = hat(wxr) + what*rhat;
    tmp3 = D1w_phid + D2D1w(phid);
    tmp4 = -hat(II*w) + what*II;
    
    dabardphi = -rhat*tmp1 - tmp2*D1w_phid;
    dabardphid = -rhat*tmp3 - tmp2*dexp_;
    dabardphidd = -rhat*dexp_;
    
    dFdphi = m*R*(-hat(abar)*dexp_ + dabardphi);
    dFdphid = m*R*dabardphid;
    dFdphidd = m*R*dabardphidd;
    
    RTddd = R'*ddd;
    dmubardphi = m*(RTddd*rr' - (rr'*RTddd)*I3)*dexp_ ...
               + II*tmp1 + tmp4*D1w_phid ...
               + ((rr'*wtBody)*I3 - wtBody*rr')*dexp_;
    dmubardphid = II*tmp3 + tmp4*dexp_;
    dmubardphidd = II*dexp_;
    
    dmudddd = m*dexp_'*rhat*R';
    dmudphi = R*(-hat(wmubar)*dexp_ + D1w(mubar) + dexp_*dmubardphi) + KR;
    dmudphid = dexp_'*dmubardphid + CR;
    dmudphidd = dexp_'*dmubardphidd;
    
    K = [KT dFdphi; zeros(3) dmudphi];
    C = [CT dFdphid; zeros(3) dmudphid];
    M = [m*I3 dFdphidd; dmudddd dmudphidd];
end

if nargout > 4 % need B, derivative of force w.r.t input
    B = [-m*RJ'; -m*dexp_'*rhat*R'*RJ'];
end

    function A = D1w(v)
        phiTv = phi'*v;
        A = a(3)*v*phi' + a(4)*mycross(v,phi)*phi' + a(2)*hat(v) ...
          + a(10)*phiTv*(phi*phi') + a(9)*phi*v' + a(9)*phiTv*I3;  
    end

    function A = D1sqw(v,dphi1) % ocmputes matrix D_1^2 w(phi,v)(dphi1,.)
        phiTv = phi'*v;
        phiphiT = phi*phi';
        vTdphi1 = v'*dphi1;
        A = (a(3)*v+a(4)*mycross(v,phi)+a(10)*phiTv*phi)*dphi1' ...
          + (phi'*dphi1)*(a(5)*v*phi' + a(6)*mycross(v,phi)*phi' ...
                             + a(4)*hat(v) + a(11)*phiTv*phiphiT ...
                             + a(10)*phi*v' + a(10)*phiTv*I3) ...
          + a(4)*mycross(v,dphi1)*phi' + a(10)*vTdphi1*phiphiT ...
          + a(9)*vTdphi1*I3 + a(10)*phiTv*dphi1*phi' ...
          + a(9)*dphi1*v';
    end

    function A = D2D1w(dphi1) % computes D_2D_1 w(phi,~)(dphi1,.)
        phiTdphi1 = phi'*dphi1;
        phiphiT = phi*phi';
        A = phiTdphi1*(a(3)*I3 - a(4)*hat(phi) + a(10)*phiphiT) ...
          - a(2)*hat(dphi1) + a(9)*(phi*dphi1' + dphi1*phi');  
    end
end
