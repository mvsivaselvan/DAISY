function [F, mu, ...
          F_ij, F_ib, F_ijd, F_ibd, ...
          mu_aj, mu_ab, mu_ajd, mu_abd] ...
       = CableForce(P, P0, Pdot, varTheta, varThetadot, R0, ...
                    rho, EA, EI, GJ, betAX, betBEND, betTOR, ...
                    sg, wg, nel, ...
                    colmat, colmat_brev, colmat_bar, ...
                    d, dbrev, dbar, ...
                    Mbar, ...
                    u, ...
                    Kbar11, Dbar11) %#ok<INUSL>
% INPUTS
% P = spline position DOF (represented as 3*N matrix), \mathscr{P}
% P0 = reference configuration
% Pd = \dot{P}
% varTheta = spline twist DOF (\vartheta)
% varThetdot = \dot{vartheta}
% R0 = cell array of orientations in the reference configuration at the
%      quadrature points (+ the two end points)
% rho = mass per unit length
% EA = section axial rigidity
% EI = section bending rigidity
% GJ = section torsional rigidity
% betAX = damping coeff associated with rate of axial deformation
% betBEND = damping coeff associated with rate of bending deformation
% betTOR = damping coeff associated with rate of torsional deformation
% sg = quadrature points
% wg = quadrature weights
% nel = nel(n) = index of element to which quadrature pt sg(n) belongs to
% colmat = B-spline basis for position and its derivatives (B)
% colmat_brev = B-spline basis for twist and its derivative (Brev)
% colmat_bar = B-spline basis for strain projection (Bbar)
% d = degree of B-spline basis for position (B)
% dbrev = degree of B-spline basis for twist (Bbrev)
% dbar = degree of B-spline basis for strain projection (Bbar)
% Mbar = "mass matrix" for strain projection
% u = force per unit mass (acceleration due to gravity and ground) 
% Kbar11, Dbar11 = as defined in equation (33) [EQ NUM needs to be updated
%                  to be consistent with document]
% OUTPUTS
% F = Force vector of length 3*N
% mu = Moment vector of length Nbrev
% F_ij, F_ib, F_ijd, F_ibd, mu_aj, mu_ab, mu_ajd, mu_abd = stiffness terms

% number of basis functions (or control points)
N = size(colmat,2);
Nbrev = size(colmat_brev,2);
Nbar = size(colmat_bar,2);

% Constitutive matrices
Em11 = EA;
Em12 = zeros(1,3);
Em22 = diag([GJ, EI, EI]);
Dm11 = betAX*EA;
Dm12 = zeros(1,3);
Dm22 = diag([betTOR*GJ, betBEND*EI, betBEND*EI]);

% number of quadrature (or collocation) points 
ncolloc = length(wg);

I3 = eye(3);

% Compute coefficients for the projected strain
bepsbar = zeros(Nbar,1);
bepsdbar = zeros(Nbar,1);
for n = 1:ncolloc
    Bp = colmat(3*(n-1)+2,:)';
    
    Bbar = colmat_bar(n,:)';
    
    xi0p = P0*Bp;
    nxi0p = norm(xi0p);
    xip = P*Bp;
    xipd = Pdot*Bp;
    nxip = norm(xip);
    tau = xip/nxip;
    
    nel_ = nel(n);
    bepsbar(nel_+(0:dbar)) = bepsbar(nel_+(0:dbar)) + ...
                             Bbar(nel_+(0:dbar))*((nxip - nxi0p)*wg(n));
    bepsdbar(nel_+(0:dbar)) = bepsdbar(nel_+(0:dbar)) + ...
                             Bbar(nel_+(0:dbar))*(tau'*xipd*wg(n));
end
cepsbar = Mbar\bepsbar;
cepsdbar = Mbar\bepsdbar;

% In the first pass, projection Nbar is computed,
% and in the second pass, force vector and stiffness are 
% assembled. The code is organized in this way because 
% there is significant repetition of code between the two
% passes
%----------------------------------------------------------------------
% The following two allocations are to define size for code generation
bNbar = zeros(Nbar, 1); % allocate, computed in pass 1
cNbar = zeros(Nbar, 1); % allocate, computed in pass 1
F_ = zeros(3*(d+1),1); % allocate, computed in pass 2
F = zeros(size(P(:))); % Generalized force
mu = zeros(size(varTheta)); % Generalized moment
if (nargout > 2)
    bepscomma = zeros(Nbar,3*N);
    bepscomma_ = zeros(dbar+1,3*(d+1)); % temp storage
    bepsdcomma = zeros(Nbar,3*N);
    bepsdcomma_ = zeros(dbar+1,3*(d+1)); % temp storage
    F_ij = zeros(3*N,3*N);
    F_ij_ = zeros(3*(d+1),3*(d+1)); % temp storage
    F_ib = zeros(3*N,Nbrev);
    F_ib_ = zeros(3*(d+1),dbrev+1); % temp storage
    F_ijd = zeros(3*N,3*N);
    F_ijd_ = zeros(3*(d+1),3*(d+1)); % temp storage
    F_ibd  = zeros(3*N,Nbrev);
    F_ibd_ = zeros(3*(d+1),dbrev+1); % temp storage
    mu_aj = zeros(Nbrev,3*N);
    mu_aj_ = zeros(dbrev+1,3*(d+1)); % temp storage
    mu_ab = zeros(Nbrev,Nbrev);
    mu_ab_ = zeros(dbrev+1,dbrev+1); % temp storage
    mu_ajd = zeros(Nbrev,3*N);
    mu_abd = zeros(Nbrev,Nbrev);
    mu_abd_ = zeros(dbrev+1,dbrev+1); % temp storage
end
%----------------------------------------------------------------------
for pass = 1:2 
    if pass == 1
        bNbar = zeros(Nbar,1);
    else % pass == 2
        F = zeros(size(P(:))); % Generalized force
        F_ = zeros(3*(d+1),1); % temp storage when computing F
        mu = zeros(size(varTheta)); % Generalized moment
        if (nargout > 2)
            bepscomma = zeros(Nbar,3*N);
            bepscomma_ = zeros(dbar+1,3*(d+1)); % temp storage
            bepsdcomma = zeros(Nbar,3*N);
            bepsdcomma_ = zeros(dbar+1,3*(d+1)); % temp storage
            F_ij = zeros(3*N,3*N);
            F_ij_ = zeros(3*(d+1),3*(d+1)); % temp storage
            F_ib = zeros(3*N,Nbrev); 
            F_ib_ = zeros(3*(d+1),dbrev+1); % temp storage
            F_ijd = zeros(3*N,3*N);
            F_ijd_ = zeros(3*(d+1),3*(d+1)); % temp storage
            F_ibd  = zeros(3*N,Nbrev);
            F_ibd_ = zeros(3*(d+1),dbrev+1); % temp storage
            mu_aj = zeros(Nbrev,3*N);
            mu_aj_ = zeros(dbrev+1,3*(d+1)); % temp storage
            mu_ab = zeros(Nbrev,Nbrev);
            mu_ab_ = zeros(dbrev+1,dbrev+1); % temp storage
            mu_ajd = zeros(Nbrev,3*N);
            mu_abd = zeros(Nbrev,Nbrev);
            mu_abd_ = zeros(dbrev+1,dbrev+1); % temp storage
        end
    end
    for n = 1:ncolloc
        nel_ = nel(n);
        
        B = colmat(3*(n-1)+1,:)';
        Bp = colmat(3*(n-1)+2,:)';
        Bpp = colmat(3*(n-1)+3,:)';

        Bbrev = colmat_brev(2*(n-1)+1,:)';
        Bbrevp = colmat_brev(2*(n-1)+2,:)';

        Bbar = colmat_bar(n,:)';

        xi0 = P0*B; %#ok<NASGU>
        xi0p = P0*Bp;
        xi0pp = P0*Bpp;
        nxi0p = norm(xi0p);
        tau0 = xi0p/nxi0p;
        K0 = mycross(xi0p,xi0pp)/nxi0p^2;

        xi = P*B; %#ok<NASGU>
        xip = P*Bp;
        xipd = Pdot*Bp;
        xipp = P*Bpp;
        xippd = Pdot*Bpp;
        nxip = norm(xip);
        tau = xip/nxip;
        K = mycross(xip,xipp)/nxip^2;
        G = (tau0'*K)*(2*tau0+tau)- (K0'*tau)*tau0;
        rr = K - K0 - G/(1+tau0'*tau);
        PP = I3 - tau*tau';
        L_tauI = PP/nxip;
        L_KI = ((-2/nxip)*tau)*K' + hat(xipp/nxip^2);
        L_GI = (L_KI*tau0)*(2*tau0+tau)' + (tau0'*K)*L_tauI ...
             - (L_tauI*K0)*tau0';
        L_rrI = L_KI + ((L_tauI*tau0)/(1+tau0'*tau)^2)*G' ...
              - L_GI/(1+tau0'*tau);
        L_KII = -hat(xip/nxip^2);
        L_GII = (L_KII*tau0)*(2*tau0+tau)';
        L_rrII = L_KII - L_GII/(1+tau0'*tau);

        thet = varTheta'*Bbrev;
        thetd = varThetadot'*Bbrev;
        thetp = varTheta'*Bbrevp;
        thetpd = varThetadot'*Bbrevp;

        st = sin(thet);
        ct = cos(thet);
        htau0 = hat(tau0);
        tau0p = (xi0pp - (tau0'*xi0pp)*tau0)/nxi0p;

        R0_ = R0(:,3*n+1:3*n+3);
        Theta = I3 + st*htau0 + (1-ct)*(htau0^2);
        ww = thetp*tau0 + st*tau0p - (1-ct)*mycross(tau0,tau0p);

        L_chiO = (mycross(rr,tau0) + tau0p)';
        L_chiI = tau0';

        % calculate epsilon, epsilond, kappa, kappad
        epsbar = cepsbar'*Bbar;
        kappa = (R0_'*(Theta'*rr + ww))/nxi0p;
        epsdbar = cepsdbar'*Bbar;
        kappad = (R0_'*Theta'*(L_rrI'*xipd + L_rrII'*xippd ...
                + L_chiO'*thetd + L_chiI'*thetpd))/nxi0p;
            
        if pass == 1
            NN = Em11*epsbar + Em12*kappa + Dm11*epsdbar + Dm12*kappad;
            bNbar(nel_+(0:dbar)) = bNbar(nel_+(0:dbar)) + ...
                                   Bbar(nel_+(0:dbar))*NN*nxi0p*wg(n);
        else % pass == 2
            Nbar = cNbar'*Bbar;
            MM = Em12'*epsbar + Em22*kappa + Dm12'*epsdbar + Dm22*kappad;
            mm = Theta*R0_*MM;
            F1 = -(rho*nxi0p)*u;
            F2 = tau*Nbar + L_rrI*mm;
            F3 = L_rrII*mm;
            mu1 = L_chiO*mm;
            mu2 = L_chiI*mm;
            for dd = 0:d
                F_(dd*3+(1:3)) = ...
                  (B(nel_+dd)*F1 + Bp(nel_+dd)*F2 + Bpp(nel_+dd)*F3)*wg(n);
            end
            F((nel_-1)*3+1:(nel_+d)*3) = F((nel_-1)*3+1:(nel_+d)*3) + F_;
            mu(nel_+(0:dbrev)) = mu(nel_+(0:dbrev)) + ...
              (Bbrev(nel_+(0:dbrev))*mu1+Bbrevp(nel_+(0:dbrev))*mu2)*wg(n);
            if nargout > 2
                thetad_htau0 = thetd*htau0;
                taud = (PP*xipd)/nxip;
                rrd = L_rrI'*xipd + L_rrII'*xippd;
                XX = mycross(mycross(rr,tau0),tau0) - mycross(tau0,tau0p);
                ThetaR0 = Theta*R0_;
                E22 = (ThetaR0*Em22*ThetaR0')/nxi0p;
                D22 = (ThetaR0*Dm22*ThetaR0')/nxi0p;
                UU = Qtil_rr_I(xipd) + Qtil_rr_IV(xippd) ...
                                     - thetad_htau0*L_rrI';
                VV = Qtil_rr_III(xipd) - thetad_htau0*L_rrII';
                WW = mycross(rrd,tau0) + Qtil_chi_O(thetd);
                YY = mycross(tau0,mm) + E22*L_chiO';
                
                %%%% the b terms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for ddr = 0:dbar
                    for ddc = 0:d
                        bepscomma_(ddr+1,ddc*3+(1:3)) = ...
                            (Bbar(nel_+ddr)*Bp(nel_+ddc)*wg(n))*tau';
                        bepsdcomma_(ddr+1,ddc*3+(1:3)) = ...
                            (Bbar(nel_+ddr)*Bp(nel_+ddc)*wg(n))*taud';
                    end
                end
                bepscomma(nel_+(0:dbar),(nel_-1)*3+1:3*(nel_+d)) ...
                  = bepscomma(nel_+(0:dbar),(nel_-1)*3+1:3*(nel_+d)) + ...
                      bepscomma_;
                bepsdcomma(nel_+(0:dbar),(nel_-1)*3+1:3*(nel_+d)) ...
                  = bepsdcomma(nel_+(0:dbar),(nel_-1)*3+1:3*(nel_+d)) + ...
                      bepsdcomma_;
                  
                %%%% F_ij %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Q1 = (Nbar/nxip)*PP + Q_rr_I(mm) + L_rrI*E22*L_rrI';
                Q2 = Q_rr_III(mm) + L_rrII*E22*L_rrI';
                Q3 = L_rrII*E22*L_rrII';
                D22UU = D22*UU;
                D22VV = D22*VV;
                Q4 = L_rrI*D22UU;
                Q5 = L_rrI*D22VV;
                Q6 = L_rrII*D22UU;
                Q7 = L_rrII*D22VV;
                for ddc = 0:d
                    for ddr = 0:d
                        F_ij_(ddr*3+(1:3),ddc*3+(1:3)) = ...
                           ((Bp(nel_+ddr)*Bp(nel_+ddc))*(Q1 + Q4) + ...
                            (Bp(nel_+ddr)*Bpp(nel_+ddc))*(Q2' + Q5) + ...
                            (Bpp(nel_+ddr)*Bp(nel_+ddc))*(Q2 + Q6) + ...
                            (Bpp(nel_+ddr)*Bpp(nel_+ddc))*(Q3 + Q7))*wg(n);
                    end
                end
                F_ij((nel_-1)*3+1:(nel_+d)*3,(nel_-1)*3+1:(nel_+d)*3) ...
                 = F_ij((nel_-1)*3+1:(nel_+d)*3,(nel_-1)*3+1:(nel_+d)*3)...
                   + F_ij_;
                
                %%%% F_ib %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                E22LchiI = E22*L_chiI';
                D22WW = D22*WW;
                Q8 = L_rrI*YY;
                Q9 = L_rrI*E22LchiI;
                Q10 = L_rrII*YY;
                Q11 = L_rrII*E22LchiI;
                Q12 = L_rrI*D22WW;
                Q13 = L_rrII*D22WW;
                for ddc = 0:dbrev
                    for ddr = 0:d
                        F_ib_(ddr*3+(1:3),ddc+1) = ...
                           ((Bp(nel_+ddr)*Bbrev(nel_+ddc))*(Q8+Q12) + ...
                            (Bp(nel_+ddr)*Bbrevp(nel_+ddc))*Q9 + ...
                            (Bpp(nel_+ddr)*Bbrev(nel_+ddc))*(Q10+Q13) + ...
                            (Bpp(nel_+ddr)*Bbrevp(nel_+ddc))*Q11)*wg(n);
                    end
                end
                F_ib((nel_-1)*3+1:(nel_+d)*3,nel_+(0:dbrev)) = ...
                    F_ib((nel_-1)*3+1:(nel_+d)*3,nel_+(0:dbrev)) + F_ib_;
                
                %%%% F_ijd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Q14 = L_rrI*D22*L_rrI';
                Q15 = L_rrI*D22*L_rrII';
                Q16 = L_rrII*D22*L_rrII';
                for ddc = 0:d
                    for ddr = 0:d
                        F_ijd_(ddr*3+(1:3),ddc*3+(1:3)) = ...
                           ((Bp(nel_+ddr)*Bp(nel_+ddc))*Q14 + ...
                            (Bp(nel_+ddr)*Bpp(nel_+ddc))*Q15 + ...
                            (Bpp(nel_+ddr)*Bp(nel_+ddc))*Q15' + ...
                            (Bpp(nel_+ddr)*Bpp(nel_+ddc))*Q16)*wg(n);
                    end
                end
                F_ijd((nel_-1)*3+1:(nel_+d)*3,(nel_-1)*3+1:(nel_+d)*3) ...
                = F_ijd((nel_-1)*3+1:(nel_+d)*3,(nel_-1)*3+1:(nel_+d)*3)...
                   + F_ijd_;
               
               %%%% F_ibd and mu_ajd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               Q17 = L_rrI*D22*L_chiO';
               Q18 = L_rrI*D22*L_chiI';
               Q19 = L_rrII*D22*L_chiO';
               Q20 = L_rrII*D22*L_chiI';
               for ddc = 0:dbrev
                   for ddr = 0:d
                       F_ibd_(ddr*3+(1:3),ddc+1) = ...
                           ((Bp(nel_+ddr)*Bbrev(nel_+ddc))*Q17 + ...
                            (Bp(nel_+ddr)*Bbrevp(nel_+ddc))*Q18 + ...
                            (Bpp(nel_+ddr)*Bbrev(nel_+ddc))*Q19 + ...
                            (Bpp(nel_+ddr)*Bbrevp(nel_+ddc))*Q20)*wg(n);
                   end
               end
               F_ibd((nel_-1)*3+1:(nel_+d)*3,nel_+(0:dbrev)) = ...
                    F_ibd((nel_-1)*3+1:(nel_+d)*3,nel_+(0:dbrev)) + F_ibd_;
                
               %%%% mu_aj %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               Q21 = L_chiO*D22UU;
               Q22 = L_chiI*D22UU;
               Q23 = L_chiO*D22VV;
               Q24 = L_chiI*D22VV;
               for ddc = 0:d
                   for ddr = 0:dbrev
                       mu_aj_(ddr+1,ddc*3+(1:3)) = ...
                       ((Bbrev(nel_+ddr)*Bp(nel_+ddc))*(Q8'+Q21) + ...
                        (Bbrevp(nel_+ddr)*Bp(nel_+ddc))*(Q9'+Q22) + ...
                        (Bbrev(nel_+ddr)*Bpp(nel_+ddc))*(Q10'+Q23) + ...
                        (Bbrevp(nel_+ddr)*Bpp(nel_+ddc))*(Q11'+Q24))*wg(n);
                   end
               end
               mu_aj(nel_+(0:dbrev),(nel_-1)*3+1:(nel_+d)*3) = ...
                   mu_aj(nel_+(0:dbrev),(nel_-1)*3+1:(nel_+d)*3) + mu_aj_;
                
               %%%% mu_ab %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               Q25 = L_chiO*E22*L_chiO' + Q_chi_O(mm);
               Q26 = L_chiO*E22LchiI;
               Q27 = L_chiI*E22LchiI;
               Q28 = L_chiO*D22WW;
               Q29 = L_chiI*D22WW;
               for ddc = 0:dbrev
                   for ddr = 0:dbrev
                       mu_ab_(ddr+1,ddc+1) = ...
                       ((Bbrev(nel_+ddr)*Bbrev(nel_+ddc))*(Q25+Q28) + ...
                        (Bbrev(nel_+ddr)*Bbrevp(nel_+ddc))*Q26 + ...
                        (Bbrevp(nel_+ddr)*Bbrev(nel_+ddc))*(Q26'+Q29) + ...
                        (Bbrevp(nel_+ddr)*Bbrevp(nel_+ddc))*Q27)*wg(n);
                   end
               end
               mu_ab(nel_+(0:dbrev),nel_+(0:dbrev)) = ...
                   mu_ab(nel_+(0:dbrev),nel_+(0:dbrev)) + mu_ab_;
               
               %%%% mu_ajd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               % computed at the end as simply F_ibd'
               
               %%%% mu_abd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               Q30 = L_chiO*D22*L_chiO';
               Q31 = L_chiO*D22*L_chiI';
               Q32 = L_chiI*D22*L_chiI';
               for ddc = 0:dbrev
                   for ddr = 0:dbrev
                       mu_abd_(ddr+1,ddc+1) = ...
                       ((Bbrev(nel_+ddr)*Bbrev(nel_+ddc))*Q30 + ...
                        (Bbrev(nel_+ddr)*Bbrevp(nel_+ddc))*Q31 + ...
                        (Bbrevp(nel_+ddr)*Bbrev(nel_+ddc))*Q31' + ...
                        (Bbrevp(nel_+ddr)*Bbrevp(nel_+ddc))*Q32)*wg(n);
                   end
               end
               mu_abd(nel_+(0:dbrev),nel_+(0:dbrev)) = ...
                   mu_abd(nel_+(0:dbrev),nel_+(0:dbrev)) + mu_abd_;
            end
            
        end
    end
    if pass == 1
        cNbar = Mbar\bNbar;
    end
end

if nargout > 2
    F_ij = F_ij + bepscomma'*Kbar11*bepscomma ...
                + bepscomma'*Dbar11*bepsdcomma;
    F_ijd = F_ijd + bepscomma'*Dbar11*bepscomma;
    mu_ajd = F_ibd';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UTILITY FUNCTIONS FOR SECOND DERIVATIVE COMPUTATIONS
% It is assumed in these functions that quantities xip, xipp, 
% their norms (nxip, nxipp), tau0, tau, PP, L_tauI,
% K0, K, L_KI, L_KII, G, L_GI, L_GII have been computed and are available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function related to tau
    function Q = Q_tau_I(eta)
        Peta = PP*eta;
        Q = -(tau*Peta' + Peta*tau' + (tau'*eta)*PP)/(nxip^2);
    end

    function Q = Qtil_tau_I(nu)
        Q = Q_tau_I(nu);
    end

% functions related to K
    function Q = Q_K_I(eta)
        LKIeta = L_KI*eta;
        Q = -(((eta'*K)/nxip)*I3 + tau*LKIeta' + LKIeta*tau')*(2/nxip);
    end

    function Q = Q_K_III(eta)
        Q = -(2/nxip)*((L_KII*eta)*tau') + hat(eta)/nxip^2;
    end

    function Q = Qtil_K_I(nu)
        Q = -(2/nxip)*((K*nu')/nxip + L_KI'*(nu*tau'+(tau'*nu)*I3));
    end

    function Q = Qtil_K_III(nu)
        Q = -(2/nxip)*(tau'*nu)*L_KII' + hat(nu)/nxip^2;
    end

    function Q = Qtil_K_IV(nu)
        Q = -(2/nxip)*(L_KII'*nu)*tau' - hat(nu)/nxip^2;
    end

% functions related to G
    function Q = Q_G_I(eta)
        LKItau0 = L_KI*tau0;
        LtauIeta = L_tauI*eta;
        Q = (eta'*(2*tau0 + tau))*Q_K_I(tau0) ...
          + (LKItau0*LtauIeta' + LtauIeta*LKItau0') ...
          + (tau0'*K)*Q_tau_I(eta) - (eta'*tau0)*Q_tau_I(K0);
    end

    function Q = Q_G_III(eta)
        Q = (L_KII*tau0)*(L_tauI*eta)' + (eta'*(2*tau0+tau))*Q_K_III(tau0);
    end

    function Q = Qtil_G_I(nu)
        LKItau0 = L_KI*tau0;
        Q = (2*tau0+tau)*(Q_K_I(tau0)*nu)' ...
          + (L_tauI'*nu)*LKItau0' + (LKItau0'*nu)*L_tauI' ...
          - tau0*(Q_tau_I(K0)*nu)' + (tau0'*K)*Qtil_tau_I(nu);
    end

    function Q = Qtil_G_III(nu)
        Q = (L_tauI'*nu)*(L_KII*tau0)'+(2*tau0+tau)*(Q_K_III(tau0)*nu)';
    end

    function Q = Qtil_G_IV(nu)
        Q = ((L_KII*tau0)'*nu)*L_tauI'+(2*tau0+tau)*(Q_K_III(tau0)'*nu)';
    end

% functions related to rr
    function Q = Q_rr_I(eta)
        tt0 = 1+tau0'*tau;
        LtauItau0 = L_tauI*tau0;
        LGIeta = L_GI*eta;
        etaG = eta'*G;
        Q = Q_K_I(eta) - Q_G_I(eta)/tt0 ...
            - (2*etaG/tt0^3)*(LtauItau0*LtauItau0') ...
            + (etaG/tt0^2)*Q_tau_I(tau0) ...
            + (LtauItau0*LGIeta' + LGIeta*LtauItau0')/tt0^2;
    end

    function Q = Q_rr_III(eta)
        tt0 = 1+tau0'*tau;
        LtauItau0 = L_tauI*tau0;
        LGIIeta = L_GII*eta;
        Q = Q_K_III(eta) - Q_G_III(eta)/tt0 + (LGIIeta*LtauItau0')/tt0^2;
    end

    function Q = Qtil_rr_I(nu)
        tt0 = 1+tau0'*tau;
        LtauItau0 = L_tauI*tau0;
        LtauItau0nu = LtauItau0'*nu;
        Q = Qtil_K_I(nu) - Qtil_G_I(nu)/tt0 ...
            + (-(2/tt0)*LtauItau0nu*(G*LtauItau0') ...
               + G*(Q_tau_I(tau0)*nu)' + (L_GI'*nu)*LtauItau0' ...
               + LtauItau0nu*L_GI')/tt0^2;
    end

    function Q = Qtil_rr_III(nu)
        tt0 = 1+tau0'*tau;
        LtauItau0 = L_tauI*tau0;
        LtauItau0nu = LtauItau0'*nu;
        Q = Qtil_K_III(nu) - Qtil_G_III(nu)/tt0 ...
            + (LtauItau0nu/tt0^2)*L_GII';
    end

    function Q = Qtil_rr_IV(nu)
        tt0 = 1+tau0'*tau;
        LtauItau0 = L_tauI*tau0;
        Q = Qtil_K_IV(nu) - Qtil_G_IV(nu)/tt0 ...
            + ((L_GII'*nu)*LtauItau0')/tt0^2;
    end

% functions related to chi
    function Q = Q_chi_O(eta)
        Q = eta'*XX;
    end

    function Q = Qtil_chi_O(nu)
        Q = nu*XX;
    end
end
