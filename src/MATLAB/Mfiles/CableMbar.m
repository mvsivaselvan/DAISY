function [Mbar, Kbar11, Dbar11] = CableMbar(P0, ...
                                            EA, betAX, ...
                                            colmat, colmat_bar, ...
                                            wg)
                                  
% INPUTS
% P0 = reference configuration
% EA = section axial rigidity
% betAX = damping coeff associated with rate of axial deformation
% colmat = B-spline basis for position and its derivatives (B)
% colmat_bar = B-spline basis for strain projection (Bbar)
% wg = quadrature weights
% OUTPUTS
% Mbar = "mass matrix" for strain projection
% Kbar11, Dbar11 = as defined in equation (33) [EQ NUM needs to be updated
%                  to be consistent with document]


Em11 = EA;
Dm11 = betAX*EA;

% number of strain projection basis functions
Nbar = size(colmat_bar,2);

% number of quadrature (or collocation) points 
ncolloc = length(wg);

Mbar = zeros(Nbar);
Kbar11 = zeros(Nbar);
Dbar11 = zeros(Nbar);

for n = 1:ncolloc
    Bp = colmat(3*(n-1)+2,:)';
    
    Bbar = colmat_bar(n,:)';
    
    xi0p = P0*Bp;
    nxi0p = norm(xi0p);
    
    Mbar = Mbar + (Bbar*(nxi0p*wg(n)))*Bbar';
    Kbar11 = Kbar11 + (Bbar*(Em11*nxi0p*wg(n)))*Bbar';
    Dbar11 = Dbar11 + (Bbar*(Dm11*nxi0p*wg(n)))*Bbar';
end

Kbar11 = Mbar\(Kbar11/Mbar);
Dbar11 = Mbar\(Dbar11/Mbar);


