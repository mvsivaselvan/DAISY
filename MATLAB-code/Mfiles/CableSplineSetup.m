function [xg, wg, nel, knots, colmat, colmat_brev, colmat_bar] ...
    = CableSplineSetup(N, d, dbrev, dbar, Ng)

% INPUTS
% N = number of control points for position (P)
% d = degree of B-spline basis functions for position
% dbrev = degree of B-spline basis functions for torsion
% dbar = degree of B-spline basis functions for strain projection
% Ng = number of Gauss integration points per element
% OUTPUTS
% xg = quadrature points
% wg = quadrature weights
% nel = nel(n) = index of element to which quadrature pt sg(n) belongs to
% knots = knot vector
% colmat = B-spline basis for position and its derivatives (B)
% colmat_brev = B-spline basis for twist and its derivative (Brev)
% colmat_bar = B-spline basis for strain projection (Bbar)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup collocation (integration) points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
knots = [zeros(1,d) linspace(0, 1, N-d+1) ones(1,d)]; % uniform
Nel = N-d; % number of elements
L = 1;
Lel = L/Nel; % element length
[xg,wg] = gaqdm(Ng); % Gauss points \in (0,pi) and weights
xg = flipud(cos(xg)); % Gauss points \in (-1,1)
xg = (1+xg)/2*Lel; % Gauss points \in (0,Lel)
wg = wg*(Lel/2); % weights multiplied by appropriate jacobian 
% replicate to nel elements
xg = repmat(xg,1,Nel) + (0:Nel-1)*Lel; % using implicit expansion for add
xg = reshape(xg, Ng*Nel, 1);
wg = repmat(wg, Nel, 1);
nel = kron(1:Nel,ones(1,Ng));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPLINE COLLOCATION MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position splines
xg_ = repmat(xg',3,1); % multiplicity 3 for each collocation point so as to
                       % get function, derivative and second derivative in
                       % collocation matrix
xg_ = xg_(:);
colmat = spcol(knots, d+1, xg_); % order d+1 (degree d) B-splines

% twist splines
knots_brev = [zeros(1,dbrev), 0:1/(N-d):1, ones(1,dbrev)];
xg_ = repmat(xg',2,1); % multiplicity 2 for each collocation point so as to
                       % get function and derivative in collocation matrix
xg_ = xg_(:);
colmat_brev = spcol(knots_brev, dbrev+1, xg_);

% strain projection splines
knots_bar = [zeros(1,dbar), 0:1/(N-d):1, ones(1,dbar)];
colmat_bar = spcol(knots_bar, dbar+1, xg); 
