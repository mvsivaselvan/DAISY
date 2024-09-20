function [q, J, Q, Qtilde, C] = CableBCTransinCoord(qbar, ...
                                         x0, RJ, RE, r, R0, ...
                                         eta, rho1,rho2)
% Computes the transformation from (d, phi, gamma) -> (p1, p2, vartheta) for
% one end of the cable; NOTE: the main difference from CableBCtrans.m is
% that this function is based on exponential coordinates for the rotation,
% rather than the rotation matrix; also, a joint coordinate system and end
% offset are included, making it a more general.
% Inputs
% qbar = [d; phi; gamm]
%   where d = displacement of joint, 
%         phi = exponential coordinates of rotation of joint
%         gamm = distance between first and second control points
% x0 = reference position of joint (3x1) vector
% RJ = rotation of joint coordinate frame with respect to global
% RE = rotation of cable end with respect to joint coordinate frame
% r = position of cable end relative to joint in joint coordinate system
%     (end offset, 3x1 vector)
% R0 = orientation of cable end in reference configuration
% eta = vector to contract with for second derivative, Q
% rho1 = vector to contract with for second derivative, Qtilde
% rho2 = vector to contract with for third derivative, C
% Outputs:
% q = [p1; p2; vartheta1]
%   where p1, p2 = positions of first and second control points
%         vartheta1 = (twist) rotation of first control point
% J = first derivative of transformation
% Q = second derivative contracted over upper index, 
%     Q_{ib,jb} = eta_i q^i_{ib,jb}
% Qtilde = second derivative contracted over one of the lower indices
%     Qtilde^{i}_{ib} = rho1^{jb} q^i_{ib,jb}
% C = third derivative contracted over two lower indices
%     C^{i}_{ib} = q^i_{ib,jb,kb}rho1^{jb}rho2^{kb}

d = qbar(1:3);
phi = qbar(4:6);
gamm = qbar(7);

if nargout == 1
    a = rothelper(phi, 0);
elseif nargout == 2
    a = rothelper(phi, 1);
elseif (nargout == 3 || nargout == 4)
    a = rothelper(phi, 2);
else
    a = rothelper(phi, 3);
end

I3 = eye(3);
hphi = hat(phi);
hphi2 = hphi^2;
Rphi = I3 + a(1)*hphi + a(2)*hphi2;

Rb = RJ*Rphi*RE;
tau = Rb(:,1);
tau0 = R0(:,1);

tau0_dot_tau = tau0'*tau;
tau0_cross_tau = mycross(tau0,tau);
R = (tau0_dot_tau)*I3 + hat(tau0_cross_tau) + ...
     (tau0_cross_tau/(1+tau0_dot_tau))*tau0_cross_tau'; 
Thet_tilde_e2 = R0'*(R'*Rb(:,2));

rotr = RJ*Rphi*r;

p1 = x0 + RJ*d + rotr;
p2 = p1 + gamm*tau;
vartheta1 = atan2(Thet_tilde_e2(3), Thet_tilde_e2(2));

q = [p1; p2; vartheta1];

if (nargout > 1) % first derivative needed
    dexp_ = I3 - a(2)*hphi + a(9)*hphi2;
    dexpT_ = dexp_';
    
    u = tau + tau0;
    nu_ = norm(u);
    v = 2*u/nu_^2;
    
    RJdexpT_ = RJ*dexpT_;
    
    rotrhat = hat(rotr);
    tauhat = hat(tau);
    
    Drotr = -rotrhat*RJdexpT_;
    Dtau = -tauhat*RJdexpT_;
    Dvarthet1 = v'*RJdexpT_;
    
    J = zeros(7);
    
    J(1:3,1:3) = RJ;
    J(4:6,1:3) = RJ;
    J(1:3,4:6) = Drotr;
    J(4:6,4:6) = Drotr + gamm*Dtau;
    J(7,4:6) = Dvarthet1;
    J(4:6,7) = tau;
end
    
if nargout > 2 % second derivative needed
    dphi1 = rho1(4:6);
    dgamm1 = rho1(7);
    
    % ws1 is short form for w(phi;dphi1)
    ws1 = RJdexpT_*dphi1;
    ws1hat = hat(ws1);
    
    % Dws1 is short form for Dw(phi;dphi1)
    phiTdphi1 = phi'*dphi1;
    dphi1xphi = mycross(dphi1,phi);
    dphi1hat = hat(dphi1);
    Dws1 = RJ*(a(3)*dphi1*phi' - a(4)*dphi1xphi*phi' ...
                - a(2)*dphi1hat + a(10)*phiTdphi1*(phi*phi') ...
                + a(9)*phi*dphi1' + a(9)*phiTdphi1*I3);
    
    % D2rotr is short form for D^2rotr(phi)(dphi1,dot)
    D2rotr = -(rotrhat*Dws1 + ws1hat*rotrhat*RJdexpT_);
    
    % D2tau1 is short form for D^2tau(phi)(dphi1,dot)
    D2tau1 = -(tauhat*Dws1 + ws1hat*tauhat*RJdexpT_);
    
    uu_ = u/nu_; % unit vector in the direction of u
    PP = I3 - uu_*uu_';
    Dv = -(2/nu_^2)*(2*PP-I3)*tauhat*RJdexpT_;
    
    D2varthet1 = ws1'*Dv + v'*Dws1;
    
    Qtilde = zeros(7);
    Qtilde(1:3,4:6) = D2rotr;
    Qtilde(4:6,4:6) = D2rotr + gamm*D2tau1 + dgamm1*Dtau;
    Qtilde(7,4:6) = D2varthet1;
    Qtilde(4:6,7) = Dtau*dphi1;
    
    F1 = eta(1:3);
    F2 = eta(4:6);
    F1F2 = F1 + F2;
    mu = eta(7);
    
    M = mycross(rotr,F1F2) + gamm*mycross(tau,F2) + mu*v;
    m = RJ'*M;
    mTphi = m'*phi;
    MTDws = a(3)*m*phi' - a(4)*mycross(phi,m)*phi' + a(2)*hat(m) ...
          + a(10)*mTphi*(phi*phi') + a(9)*mTphi*I3 + a(9)*phi*m';
    
    Q = zeros(7);
    Q(4:6,4:6) = MTDws ...
               + RJdexpT_'*( (rotr*(F1+F2)'+(gamm*tau)*F2' ...
                         - (rotr'*(F1+F2)+gamm*(tau'*F2))*I3 )*RJdexpT_ ...
                         + mu*Dv );
    Q(7,4:6) = F2'*Dtau;
    Q(4:6,7) = Q(7,4:6)';
end

if nargout > 4 % third derivtive needed
    dphi2 = rho2(4:6);
    dgamm2 = rho2(7);
    
    % ws2 is short form for ws(phi;dphi2)
    ws2 = RJdexpT_*dphi2;
    ws2hat = hat(ws2);
    
    % Dws2 is short form for Dws(phi;dphi2)
    phiTdphi2 = phi'*dphi2;
    dphi2xphi = mycross(dphi2,phi);
    dphi2hat = hat(dphi2);
    Dws2 = RJ*(a(3)*dphi2*phi' - a(4)*dphi2xphi*phi' ...
                - a(2)*dphi2hat + a(10)*phiTdphi2*(phi*phi') ...
                + a(9)*phi*dphi2' + a(9)*phiTdphi2*I3);
    
    % D2ws1 is short form for D^2 ws(phi;dphi1)(dphi2,dot)
    dphi1Tdphi2 = dphi1'*dphi2;
    D2ws1 = RJ*(a(5)*phiTdphi2*dphi1*phi' + a(3)*dphi1*dphi2' ...
               - a(6)*phiTdphi2*dphi1xphi*phi' - a(4)*dphi1xphi*dphi2' ...
                         - a(4)*phiTdphi2*dphi1hat ...
               - a(4)*mycross(dphi1,dphi2)*phi' ...
               + a(11)*phiTdphi1*phiTdphi2*(phi*phi') ...
                    + a(10)*phiTdphi2*phi*dphi1' ...
                    + a(10)*phiTdphi1*phi*dphi2' ...
                    + a(10)*phiTdphi1*phiTdphi2*I3 ...
               + a(10)*dphi1Tdphi2*(phi*phi') + a(9)*dphi1Tdphi2*I3 ...
               + a(10)*phiTdphi1*dphi2*phi' + a(9)*dphi2*dphi1');
           
   % D2v is short form for D^2 v(phi)(dphi2,dot)
   Dtaudphi2 = Dtau*dphi2;
   D2v = -(4/nu_^4)*((2*PP-I3)*Dtaudphi2*u' ...
                     +((u'*Dtaudphi2)*I3+u*Dtaudphi2')*PP)*Dtau ...
         +(2/nu_^2)*(2*PP-I3)*(-tauhat*Dws2+ws2hat*Dtau);
   
   % D3rotr is short form for D^3 rotr(phi)(dphi1,dphi2,dot)
   D3rotr = -(rotrhat*D2ws1 + hat(Dws1*dphi2)*rotrhat*RJdexpT_ ...
              + hat(mycross(ws2,rotr))*Dws1 ...
              + ws1hat*rotrhat*Dws2 ...
              +ws1hat*ws2hat*rotrhat*RJdexpT_);
          
   % D3tau is short form for D^3 tau(phi)(dphi1,dphi2,dot)
   D3tau = -(tauhat*D2ws1 + hat(Dws1*dphi2)*tauhat*RJdexpT_ ...
              + hat(mycross(ws2,tau))*Dws1 ...
              + ws1hat*tauhat*Dws2 ...
              +ws1hat*ws2hat*tauhat*RJdexpT_);
          
   % D3varthet1 is short form for D^3 varthet1(phi)(dphi1,dphi2,dot)
   D3varthet1 = (Dv*dphi2)'*Dws1 + ws1'*D2v ...
              + (Dws1*dphi2)'*Dv + v'*D2ws1;
          
   % D2tau2 is short form for D^2tau(phi)(dphi2,dot)
   D2tau2 = -(tauhat*Dws2 + ws2hat*tauhat*RJdexpT_);
          
   C = zeros(7);
   C(1:3,4:6) = D3rotr;
   C(4:6,4:6) = D3rotr + gamm*D3tau + dgamm1*D2tau2 + dgamm2*D2tau1;
   C(7,4:6) = D3varthet1;
   C(4:6,7) = D2tau1*dphi2;
end
