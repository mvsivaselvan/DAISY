function dFdu = CableForceInputDerivative(P0, ...
                                          rho, ...
                                          wg, nel, colmat, d)
% Computes the derivative of cable force, F, w.r.t. input, u, i.e. dF/du
% Since dmu/du = 0, it is not computed.
% INPUTS
% P0 = reference configuration
% rho = mass per unit length
% wg = quadrature weights
% nel = nel(n) = index of element to which quadrature pt sg(n) belongs to
% colmat = B-spline basis for position and its derivatives (B)
% d = degree of B-spline basis for position (B)
% OUTPUT
% dFdu = dF/du (a 3Nx3 matrix, since there are 3N force components and
%               3 input components)

% number of basis functions (or control points)
N = size(colmat,2);

% number of quadrature (or collocation) points 
ncolloc = length(wg);

dFdu = zeros(3*N, 3); % 3DOF per control point, 3 inputs
dFdu_ = zeros(3*(d+1),3); % temp storage when computing dFdu
for n = 1:ncolloc
    nel_ = nel(n);
    
    B = colmat(3*(n-1)+1,:)';
    Bp = colmat(3*(n-1)+2,:)';
    
    xi0p = P0*Bp;
    nxi0p = norm(xi0p);
    
    for dd = 0:d
        dFdu_(dd*3+(1:3),:) = -(B(nel_+dd)*rho*nxi0p*wg(n))*eye(3);
    end
    dFdu((nel_-1)*3+1:(nel_+d)*3,:)=dFdu((nel_-1)*3+1:(nel_+d)*3,:)+dFdu_;
end
