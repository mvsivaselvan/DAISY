function [F, mu, ...
          F_ij, F_ib, F_ijd, F_ibd, F_ijdd, F_ibdd, ...
          mu_aj, mu_ab, mu_ajd, mu_abd, mu_ajdd, mu_abdd] ...
       = CableInertiaForce(P, P0, Pdot, Pddot, ...
                           varTheta, varThetadot, varThetaddot, ...
                           R0, ...
                           rho, II, ...
                           sg, wg, nel, ...
                           colmat, colmat_brev, ...
                           d, dbrev, alph0)                 %#ok<INUSL>
% INPUTS
% P = spline position DOF (represented as 3*N matrix), \mathscr{P}
% P0 = reference configuration
% Pdot = \dot{P}
% Pddot = \ddot{P}
% varTheta = spline twist DOF (\vartheta)
% varThetadot = \dot{vartheta}
% varThetaddot = \ddot{vartheta}
% R0 = cell array of orientations in the reference configuration at the
%      quadrature points (+ the two end points)
% rho = mass per unit length
% II = body frame mass moment of inertia (3x3 matrix) per unit length
% sg = quadrature points
% wg = quadrature weights
% nel = nel(n) = index of element to which quadrature pt sg(n) belongs to
% colmat = B-spline basis for position and its derivatives (B)
% colmat_brev = B-spline basis for twist and its derivative (Brev)
% d = degree of B-spline basis for position (B)
% dbrev = degree of B-spline basis for twist (Bbrev)
% alph0 = mass-proportional damping coefficient (JB Model)
% OUTPUTS
% F = Vector of inertia forces length 3*N
% mu = Vector of inertia moments length Nbrev
% F_ij, F_ib, F_ijd, F_ibd, F_ijdd, F_ibdd, ...
% mu_aj, mu_ab, mu_ajd, mu_abd, mu_ajdd, mu_abdd = derivatives

% number of basis functions (or control points)
N = size(colmat,2);
Nbrev = size(colmat_brev,2);

% number of quadrature (or collocation) points 
ncolloc = length(wg);

if nargin<18
    alph0 = 0;
end
e1 = [1; 0; 0];
I3 = eye(3);

F = zeros(size(P(:))); % Generalized centrifugal+coriolis force
F_ = zeros(3*(d+1),1); % temp storage when computing F
mu = zeros(size(varTheta)); % Generalized centrifugal+coriolis moment
if (nargout > 4)
    F_ij = zeros(3*N,3*N);
    F_ij_ = zeros(3*(d+1),3*(d+1)); % temp storage
    F_ib = zeros(3*N,Nbrev);
    F_ib_ = zeros(3*(d+1),dbrev+1); % temp storage
    F_ijd = zeros(3*N,3*N);
    F_ijd_ = zeros(3*(d+1),3*(d+1)); % temp storage
    F_ibd  = zeros(3*N,Nbrev);
    F_ibd_ = zeros(3*(d+1),dbrev+1); % temp storage
    F_ijdd = zeros(3*N,3*N);
    F_ijdd_ = zeros(3*(d+1),3*(d+1)); % temp storage
    F_ibdd  = zeros(3*N,Nbrev);
    F_ibdd_ = zeros(3*(d+1),dbrev+1); % temp storage
    mu_aj = zeros(Nbrev,3*N);
    mu_aj_ = zeros(dbrev+1,3*(d+1)); % temp storage
    mu_ab = zeros(Nbrev,Nbrev);
    mu_ab_ = zeros(dbrev+1,dbrev+1); % temp storage
    mu_ajd = zeros(Nbrev,3*N);
    mu_ajd_ = zeros(dbrev+1,3*(d+1)); % temp storage
    mu_abd = zeros(Nbrev,Nbrev);
    mu_abd_ = zeros(dbrev+1,dbrev+1); % temp storage
    mu_ajdd = zeros(Nbrev,3*N);
    mu_ajdd_ = zeros(dbrev+1,3*(d+1)); % temp storage
    mu_abdd = zeros(Nbrev,Nbrev);
    mu_abdd_ = zeros(dbrev+1,dbrev+1); % temp storage
end

for n = 1:ncolloc
    nel_ = nel(n);
    
    B = colmat(3*(n-1)+1,:)';
    Bp = colmat(3*(n-1)+2,:)';
    
    Bbrev = colmat_brev(2*(n-1)+1,:)';
    
    xi0 = P0*B; %#ok<NASGU>
    xi0p = P0*Bp;
    nxi0p = norm(xi0p);
    tau0 = xi0p/nxi0p;
    
    xid = Pdot*B;
    xidd = Pddot*B; 
    xip = P*Bp;
    xipd = Pdot*Bp;
    xipdd = Pddot*Bp;
    nxip = norm(xip);
    tau = xip/nxip;
    
    thet = varTheta'*Bbrev;
    thetd = varThetadot'*Bbrev;
    thetdd = varThetaddot'*Bbrev;
    
    st = sin(thet);
    ct = cos(thet);
    htau0 = hat(tau0);

    R0_ = R0(:,3*n+1:3*n+3);
    Theta = I3 + st*htau0 + (1-ct)*(htau0^2);
    
    PP = I3 - tau*tau';
    tau_dot_tau0 = tau'*tau0;
    tau0_x_tau = mycross(tau0, tau);
    G = hat(tau) - (2*tau0+tau)*tau0_x_tau'/(1+tau_dot_tau0);
    W = (1/nxip)*R0_'*Theta'*G;
    l = II*(W*xipd + e1*thetd); % angular momentum
    
    G1 = (2*tau0 + tau)*tau0'/(1+tau_dot_tau0) - I3;
    M1 = htau0/(1+tau_dot_tau0);
    M2 = (2*tau0+tau);
    M3 = (tau*tau0')/(1+tau_dot_tau0)-I3;
    
    alph = W*xipdd + e1*thetdd;
    Ialph = II*alph;
    
    L_WTp_l = L_WTp(l);
    L_WTthet_l = L_WTthet(l);
    L_Wp_xipd = L_Wp(xipd);
    L_Wthet_xipd = L_Wthet(xipd);
    
    Wd_xipd = L_Wp_xipd*xipd + L_Wthet_xipd*thetd; % Wdot*xipd
    IWd_xipd = II*Wd_xipd;
    F0 = rho*xidd;
    F1 = W'*(Ialph + IWd_xipd) + L_WTp_l*xipd + L_WTthet_l*thetd ...
         - L_Wp_xipd'*l + alph0*rho*xid;
    
    for dd = 0:d
        F_(dd*3+(1:3)) = (B(nel_+dd)*F0 + Bp(nel_+dd)*F1)*nxi0p*wg(n);
    end
    F((nel_-1)*3+1:(nel_+d)*3) = F((nel_-1)*3+1:(nel_+d)*3) + F_;
    
    mu0 = Ialph(1) + IWd_xipd(1) - L_Wthet_xipd'*l;
    mu(nel_+(0:dbrev)) = mu(nel_+(0:dbrev)) + ...
                            Bbrev(nel_+(0:dbrev))*mu0*nxi0p*wg(n);
            
    if nargout > 4
        L_lp = II*L_Wp_xipd;
        L_lthet = II*L_Wthet_xipd;
        L_lpd = II*W;
        L_lthetd = II(:,1); 
        Q2_WTp_xipd = Q2_WTp(xipd);
        Q2_WTthet_thetd = Q2_WTthet(thetd);
        Wd_xipd_p = Q1_Wpp(xipd,xipd) + Q1_Wthetp(xipd,thetd);
        Wd_xipd_thet = Q1_Wpthet(xipd,xipd) + Q1_Wthetthet(xipd,thetd);
        Q2_Wp_xipd = Q2_Wp(xipd);
        Q2_Wthet_thetd = Q2_Wthet(thetd);
        
        %%%% F_ij %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Q1 = L_WTp(Ialph) + W'*II*L_Wp(xipdd) ...
           + Q1_WTpp(l,xipd) + Q2_WTp_xipd*L_lp ...
           + Q1_WTthetp(l,thetd) + Q2_WTthet_thetd*L_lp ...
           + L_WTp(IWd_xipd) ...
           + W'*II*Wd_xipd_p ...
           - Q1T_Wpp(xipd,l) - L_Wp_xipd'*L_lp;
        for ddc = 0:d
            for ddr = 0:d
                F_ij_(ddr*3+(1:3),ddc*3+(1:3)) = ...
                    (Bp(nel_+ddr)*Bp(nel_+ddc))*Q1*nxi0p*wg(n);
            end
        end
        F_ij((nel_-1)*3+1:(nel_+d)*3,(nel_-1)*3+1:(nel_+d)*3) ...
            = F_ij((nel_-1)*3+1:(nel_+d)*3,(nel_-1)*3+1:(nel_+d)*3)...
            + F_ij_;
        
        %%%% F_ib %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Q1 = L_WTthet(Ialph) + W'*II*L_Wthet(xipdd) ...
           + Q1_WTpthet(l,xipd) + Q2_WTp_xipd*L_lthet ...
           + Q1_WTthetthet(l,thetd) + Q2_WTthet_thetd*L_lthet ...
           + L_WTthet(IWd_xipd) ...
           + W'*II*Wd_xipd_thet ...
           - Q1T_Wpthet(xipd,l) - L_Wp_xipd'*L_lthet;
        for ddc = 0:dbrev
            for ddr = 0:d
                F_ib_(ddr*3+(1:3),ddc+1) = ...
                    (Bp(nel_+ddr)*Bbrev(nel_+ddc))*Q1*nxi0p*wg(n);
            end
        end
        F_ib((nel_-1)*3+1:(nel_+d)*3,nel_+(0:dbrev)) = ...
            F_ib((nel_-1)*3+1:(nel_+d)*3,nel_+(0:dbrev)) + F_ib_;
        
        %%%% F_ijd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Q1 = Q2_WTp_xipd*L_lpd + L_WTp_l + Q2_WTthet_thetd*L_lpd ...
           + W'*II*(Q2_Wp_xipd + L_Wp_xipd + Q2_Wthet_thetd) ...
           - Q2T_Wp(l) - L_Wp_xipd'*L_lpd;
        for ddc = 0:d
            for ddr = 0:d
                F_ijd_(ddr*3+(1:3),ddc*3+(1:3)) = ...
                    ( (B(nel_+ddr)*B(nel_+ddc))*alph0*rho*I3 ...
                    +Bp(nel_+ddr)*Bp(nel_+ddc))*Q1*nxi0p*wg(n);
            end
        end
        F_ijd((nel_-1)*3+1:(nel_+d)*3,(nel_-1)*3+1:(nel_+d)*3) ...
            = F_ijd((nel_-1)*3+1:(nel_+d)*3,(nel_-1)*3+1:(nel_+d)*3)...
            + F_ijd_;
        
        %%%% F_ibd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Q1 = Q2_WTp_xipd*L_lthetd + Q2_WTthet_thetd*L_lthetd + L_WTthet_l...
           + W'*L_lthet - L_Wp_xipd'*L_lthetd; 
        for ddc = 0:dbrev
            for ddr = 0:d
                F_ibd_(ddr*3+(1:3),ddc+1) = ...
                    (Bp(nel_+ddr)*Bbrev(nel_+ddc))*Q1*nxi0p*wg(n);
            end
        end
        F_ibd((nel_-1)*3+1:(nel_+d)*3,nel_+(0:dbrev)) = ...
            F_ibd((nel_-1)*3+1:(nel_+d)*3,nel_+(0:dbrev)) + F_ibd_;
        
        %%%% F_ijdd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Q1 = W'*II*W;
        for ddc = 0:d
            for ddr = 0:d
                F_ijdd_(ddr*3+(1:3),ddc*3+(1:3)) = ...
                    ( (B(nel_+ddr)*B(nel_+ddc))*rho*I3 ...
                     +(Bp(nel_+ddr)*Bp(nel_+ddc))*Q1)*nxi0p*wg(n);
            end
        end
        F_ijdd((nel_-1)*3+1:(nel_+d)*3,(nel_-1)*3+1:(nel_+d)*3) ...
            = F_ijdd((nel_-1)*3+1:(nel_+d)*3,(nel_-1)*3+1:(nel_+d)*3)...
            + F_ijdd_;
        
        %%%% F_ibdd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Q1 = W'*II(:,1); 
        for ddc = 0:dbrev
            for ddr = 0:d
                F_ibdd_(ddr*3+(1:3),ddc+1) = ...
                    (Bp(nel_+ddr)*Bbrev(nel_+ddc))*Q1*nxi0p*wg(n);
            end
        end
        F_ibdd((nel_-1)*3+1:(nel_+d)*3,nel_+(0:dbrev)) = ...
            F_ibdd((nel_-1)*3+1:(nel_+d)*3,nel_+(0:dbrev)) + F_ibdd_;
        
        %%%% mu_aj %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Q1 = e1'*II*Wd_xipd_p - Q1T_Wthetp(xipd,l) - L_Wthet_xipd'*L_lp...
            + e1'*II*L_Wp(xipdd);
        for ddc = 0:d
            for ddr = 0:dbrev
                mu_aj_(ddr+1,ddc*3+(1:3)) = ...
                    (Bbrev(nel_+ddr)*Bp(nel_+ddc))*Q1*nxi0p*wg(n);
            end
        end
        mu_aj(nel_+(0:dbrev),(nel_-1)*3+1:(nel_+d)*3) = ...
            mu_aj(nel_+(0:dbrev),(nel_-1)*3+1:(nel_+d)*3) + mu_aj_;
        
        %%%% mu_ab %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Q1 = e1'*II*Wd_xipd_thet - Q1T_Wthetthet(xipd,l) ...
           - L_Wthet_xipd'*L_lthet + e1'*II*L_Wthet(xipdd);
        for ddc = 0:dbrev
            for ddr = 0:dbrev
                mu_ab_(ddr+1,ddc+1) = ...
                    (Bbrev(nel_+ddr)*Bbrev(nel_+ddc))*Q1*nxi0p*wg(n);
            end
        end
        mu_ab(nel_+(0:dbrev),nel_+(0:dbrev)) = ...
            mu_ab(nel_+(0:dbrev),nel_+(0:dbrev)) + mu_ab_;
        
        %%%% mu_ajd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Q1 = e1'*II*(Q2_Wp_xipd + L_Wp_xipd + Q2_Wthet_thetd) ...
           - Q2T_Wthet(l) - L_Wthet_xipd'*L_lpd;
        for ddc = 0:d
            for ddr = 0:dbrev
                mu_ajd_(ddr+1,ddc*3+(1:3)) = ...
                    (Bbrev(nel_+ddr)*Bp(nel_+ddc))*Q1*nxi0p*wg(n);
            end
        end
        mu_ajd(nel_+(0:dbrev),(nel_-1)*3+1:(nel_+d)*3) = ...
            mu_ajd(nel_+(0:dbrev),(nel_-1)*3+1:(nel_+d)*3) + mu_ajd_;
        
        %%%% mu_abd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Q1 = e1'*L_lthet - L_Wthet_xipd'*L_lthetd;
        for ddc = 0:dbrev
            for ddr = 0:dbrev
                mu_abd_(ddr+1,ddc+1) = ...
                    (Bbrev(nel_+ddr)*Bbrev(nel_+ddc))*Q1*nxi0p*wg(n);
            end
        end
        mu_abd(nel_+(0:dbrev),nel_+(0:dbrev)) = ...
            mu_abd(nel_+(0:dbrev),nel_+(0:dbrev)) + mu_abd_;
        
        %%%% mu_ajdd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Q1 = II(1,:)*W;
        for ddc = 0:d
            for ddr = 0:dbrev
                mu_ajdd_(ddr+1,ddc*3+(1:3)) = ...
                    (Bbrev(nel_+ddr)*Bp(nel_+ddc))*Q1*nxi0p*wg(n);
            end
        end
        mu_ajdd(nel_+(0:dbrev),(nel_-1)*3+1:(nel_+d)*3) = ...
            mu_ajdd(nel_+(0:dbrev),(nel_-1)*3+1:(nel_+d)*3) + mu_ajdd_;
        
        %%%% mu_abdd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Q1 = II(1,1);
        for ddc = 0:dbrev
            for ddr = 0:dbrev
                mu_abdd_(ddr+1,ddc+1) = ...
                    (Bbrev(nel_+ddr)*Bbrev(nel_+ddc))*Q1*nxi0p*wg(n);
            end
        end
        mu_abdd(nel_+(0:dbrev),nel_+(0:dbrev)) = ...
            mu_abdd(nel_+(0:dbrev),nel_+(0:dbrev)) + mu_abdd_;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS FOR DERIVATIVE COMPUTATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function LL = L_Wp(upsilon)
        LL = (1/nxip)*(-(W*upsilon)*tau' + R0_'*Theta'*L_Gp(upsilon));
    end

    function LL = L_Gp(upsilon)
        G2 = G2_(upsilon);
        LL = (1/nxip)*(G1*G2*PP);
    end

    function G2 = G2_(upsilon)
        G2 = hat(upsilon) + (tau0_x_tau'*upsilon)/(1+tau_dot_tau0)*I3;
    end

    function LL = L_WTp(upsilon)
        LL = (1/nxip)*(-(W'*upsilon)*tau' + L_GTp(Theta*R0_*upsilon));
    end

    function LL = L_GTp(upsilon)
        G3 = G3_(upsilon);
        LL = (1/nxip)*(G3*PP);
    end

    function G3 = G3_(upsilon)
        G3 = hat(upsilon) + M1*((M2'*upsilon)*M3-tau*upsilon');
    end

    function LL = L_Wthet(upsilon)
        LL = -mycross(e1,W*upsilon);
    end

    function LL = L_WTthet(upsilon)
        LL = W'*mycross(e1,upsilon);
    end

    function QQ = Q1_Wpp(upsilon,pie)
        QQ = -L_Wp(upsilon)*(pie*tau' + (tau'*pie)*I3) ...
           - (1/nxip)*(W*upsilon)*pie'*PP ...
           + R0_'*Theta'*Q1_Gpp(upsilon,pie);
        QQ = QQ/nxip;
    end

    function QQ = Q2_Wp(pie)
        QQ = (-(tau'*pie)*W + R0_'*Theta'*Q2_Gp(pie))/nxip;
    end

    function QQ = Q1_Gpp(upsilon,pie)
        G2 = G2_(upsilon);
        Ppie = PP*pie;
        QQ = -L_Gp(upsilon)*pie*tau' + L_G1p(G2*Ppie) ...
           + G1*L_G2p(Ppie,upsilon) + G1*G2*L_PPp(pie);
        QQ = QQ/nxip;
    end

    function QQ = Q2_Gp(pie)
        Ppie = PP*pie;
        QQ = (G1/nxip)*(-hat(Ppie) + (Ppie*tau0_x_tau')/(1+tau_dot_tau0));
    end

    function LL = L_G1p(upsilon)
        LL = -(tau0'*upsilon)/nxip/(1+tau_dot_tau0)*G1*PP;
    end

    function LL = L_G2p(pie,upsilon)
        G2 = G2_(upsilon);
        LL = -1/nxip/(1+tau_dot_tau0)*pie*tau0'*G2*PP;
    end

    function LL = L_PPp(upsilon)
        LL = -1/nxip*(tau*upsilon' + (tau'*upsilon)*I3)*PP;
    end

    function QQ = Q1_Wpthet(upsilon,pie)
        QQ = -1/nxip*((tau'*pie)*L_Wthet(upsilon) ...
                      + mycross(e1,R0_'*Theta'*L_Gp(upsilon)*pie));
    end

    function QQ = Q1T_Wpp(upsilon,pie)
        mm = tau*pie'*L_Wp(upsilon);
        QQ = -(mm + mm') - (pie'*W*upsilon)/nxip*PP ...
           + Q1T_Gpp(upsilon,Theta*R0_*pie);
        QQ = QQ/nxip;
    end

    function QQ = Q2T_Wp(pie)
        QQ = 1/nxip*(-tau*pie'*W + Q2T_Gp(Theta*R0_*pie));
    end

    function QQ = Q1T_Gpp(upsilon,pie)
        G2 = G2_(upsilon);
        G1Tpie = G1'*pie;
        QQ = -L_Gp(upsilon)'*pie*tau' + L_PPp(G2'*G1Tpie) ...
           + PP*(L_G2p(G1Tpie,upsilon) + G2'*L_G1Tp(pie));
        QQ = QQ/nxip;
    end

    function QQ = Q2T_Gp(pie)
        G1Tpie = G1'*pie;
        QQ = 1/nxip*PP*(hat(G1Tpie)+(G1Tpie*tau0_x_tau')/(1+tau_dot_tau0));
    end

    function LL = L_G1Tp(upsilon)
        LL = -1/nxip/(1+tau_dot_tau0)*tau0*upsilon'*G1*PP;
    end

    function QQ = Q1T_Wpthet(upsilon,pie)
        QQ = -1/nxip*(tau*pie'*L_Wthet(upsilon) ...
            - L_Gp(upsilon)'*Theta*R0_*mycross(e1,pie));
    end

    function QQ = Q1_WTpp(upsilon,pie)
        QQ = -L_WTp(upsilon)*(pie*tau' + (pie'*tau)*I3) ...
           - (W'*upsilon)*pie'*PP/nxip + Q1_GTpp(Theta*R0_*upsilon,pie);
        QQ = QQ/nxip;
    end

    function QQ = Q2_WTp(pie)
        QQ = 1/nxip*(-(tau'*pie)*W' + Q2_GTp(pie)*Theta*R0_);
    end

    function QQ = Q1_GTpp(upsilon,pie)
        G3 = G3_(upsilon);
        QQ = -L_GTp(upsilon)*pie*tau'+L_G3p(PP*pie,upsilon)+G3*L_PPp(pie);
        QQ = QQ/nxip;
    end

    function QQ = Q2_GTp(pie)
        Ppie = PP*pie;
        QQ = -hat(Ppie) ...
           + M1*(((tau0'*Ppie)/(1+tau_dot_tau0)*tau-Ppie)*M2'-tau*pie'*PP);
        QQ = QQ/nxip;
    end
        
    function LL = L_G3p(pie,upsilon)
        LL = (M1/nxip)*M3*((pie*upsilon' + (pie'*upsilon)*I3) ...
                           -(M2'*upsilon)/(1+tau_dot_tau0)*...
                             (pie*tau0'+(pie'*tau0)*I3))*PP;
    end

    function QQ = Q1_WTpthet(upsilon,pie)
        QQ = 1/nxip*(-(tau'*pie)*L_WTthet(upsilon) ...
                     + Q2_GTp(pie)*Theta*R0_*mycross(e1,upsilon));
    end

    function QQ = Q1_Wthetp(upsilon,pie)
        QQ = -hat(e1)*L_Wp(upsilon)*pie;
    end

    function QQ = Q2_Wthet(pie)
        % QQ = -mycross(e1,W*pie);
        QQ = -hat(e1)*W*pie;
    end

    function QQ = Q1_Wthetthet(upsilon,pie)
        QQ = -hat(e1)*L_Wthet(upsilon)*pie;
    end

    function QQ = Q1_WTthetp(upsilon,pie)
        QQ = L_WTp(mycross(e1,upsilon))*pie;
    end

    function QQ = Q2_WTthet(pie)
        QQ = W'*hat(e1)*pie;
    end

    function QQ = Q1_WTthetthet(upsilon,pie)
        QQ = L_WTthet(mycross(e1,upsilon))*pie;
    end

    function QQ = Q1T_Wthetp(upsilon,pie)
        QQ = mycross(e1,pie)'*L_Wp(upsilon);
    end

    function QQ = Q2T_Wthet(pie)
        QQ = mycross(e1,pie)'*W;
    end

    function QQ = Q1T_Wthetthet(upsilon,pie)
        QQ = mycross(e1,pie)'*L_Wthet(upsilon);
    end
end
