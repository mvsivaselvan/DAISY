function [eps, eps_projected, kappa, Curv0, Curv] ...
             = CableStrains(P0, P, varTheta, ...
                            colmat, colmat_brev, colmat_bar,...
                            nel, wg, Mbar, dbar, R0, dynamic) %#ok<INUSL>
% Calculates axial strain and curvatures for conductor
% INPUTS
% P0 = undeformed (reference) spline position DOF
% P = deformed spline position DOF,\mathscr{P}
%   static:  3xN matrix
%   dynamic: 3NxNsteps (=length of time array (number of deformed configs))
% varTheta = spline twist DOF (\vartheta)
% colmat = B-spline basis for position and its derivatives (B)
% colmat_brev = B-spline basis for twist and its derivative (Brev)
% colmat_bar = B-spline basis for strain projection (Bbar)
% nel = nel(n) = index of element to which quadrature pt sg(n) belongs to
% wg = quadrature weights
% Mbar = "mass matrix" for strain projection
% dbar = degree of B-spline basis for strain projection (Bbar)
% R0 = cell array of orientations in the reference configuration at the
%      quadrature points
% dynamic = binary argument to define static(0)/dynamic(1) state
% OUTPUTS 
% (all computed at quadrature points of conductor)
% eps = axial strain 
% eps_projected = projected axial strain 
% kappa = curvature change vector 
% Curv0 = curvature vector of reference configuration
% Curv = curvature vector of deformed configuration

if dynamic==1
    Nsteps = size(P,2);
    N = size(P,1)/3;
else % static case
    Nsteps = 1; 
    N = size(P,2);
    P_ = P;
    varTheta_ = varTheta;
end
ncolloc = size(colmat_bar,1);

% initialization
eps = zeros(ncolloc,Nsteps);
eps_projected = zeros(ncolloc,Nsteps);
kappa = zeros(3*ncolloc,Nsteps);
Curv0 = zeros(3*ncolloc,Nsteps);
Curv = zeros(3*ncolloc,Nsteps);

for k = 1:Nsteps
    % spline control points at each time step
    if dynamic==1
        P_ = reshape(P(:,k),3,N);
        varTheta_ = varTheta(:,k);
    end
    % coefficients for axial strain projection 
    bepsbar = zeros(size(colmat_bar,2),1);
    for n = 1:ncolloc
        Bp = colmat(3*(n-1)+2,:)';
        
        Bbar = colmat_bar(n,:)';
        
        xi0p = P0*Bp;
        nxi0p = norm(xi0p);
        xip = P_*Bp;
        nxip = norm(xip);
        nel_ = nel(n);
        bepsbar(nel_+(0:dbar)) = bepsbar(nel_+(0:dbar)) + ...
                                Bbar(nel_+(0:dbar))*((nxip - nxi0p)*wg(n));
    end
    cepsbar = Mbar\bepsbar;
    % compute axial strain and curvatures at each time step    
    for n = 1:ncolloc
        B = colmat(3*(n-1)+1,:)';
        Bp = colmat(3*(n-1)+2,:)';
        Bpp = colmat(3*(n-1)+3,:)';
        
        Bbrev = colmat_brev(2*(n-1)+1,:)';
        Bbrevp = colmat_brev(2*(n-1)+2,:)';
        
        Bbar = colmat_bar(n,:)';
        
        % xi0 = P0*B;
        xi0p = P0*Bp; 
        xi0pp = P0*Bpp;
        nxi0p = norm(xi0p);
        tau0 = xi0p/nxi0p;
        K0 = mycross(xi0p,xi0pp)/nxi0p^2;
        
        % xi = P_*B;
        xip = P_*Bp;
        xipp = P_*Bpp;
        nxip = norm(xip);
        tau = xip/nxip;
        K = mycross(xip,xipp)/nxip^2;
       
        G = (tau0'*K)*(2*tau0+tau)- (K0'*tau)*tau0;
        rr = K - K0 - G/(1+tau0'*tau);
            
        thet = varTheta_'*Bbrev;
        thetp = varTheta_'*Bbrevp;
        
        st = sin(thet);
        ct = cos(thet);
        htau0 = hat(tau0);
        tau0p = (xi0pp - (tau0'*xi0pp)*tau0)/nxi0p;
    
        Theta = eye(3) + st*htau0 + (1-ct)*(htau0^2);
        ww = thetp*tau0 + st*tau0p - (1-ct)*mycross(tau0,tau0p);
        
        % axial strain
        eps(n,k) = (nxip/nxi0p)-1;
        eps_projected(n,k) = Bbar'*cepsbar;
        
        % curvatures
        kappa(3*n-2:3*n,k) = (1/nxi0p)*R0(:,3*n+1:3*n+3)'*(Theta'*rr+ww);
        
        PP0 = eye(3)-tau0*tau0';
        % Ref curvature derived from framing undeformed curve
        Curv0(3*n-2:3*n,k) = mycross(xi0p,PP0*xi0pp)/(nxi0p^3); 
        Curv(3*n-2:3*n,k) = kappa(3*n-2:3*n,k) + Curv0(3*n-2:3*n,k);
    end
end