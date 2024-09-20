function [p, err] = SplineApproximation(gamm_, J_, N, xg, wg, colmat)
% Approximate a smooth curve by a splne curve
% Inputs
% gamm = function handle for a parametric curve
% N = number of control points in approximating spline curve
% xg = quadrature points to sample desired curve
% wg = quadrature weights
% colmat = B-Spline basis function evaluated at quadrature points
%          The matrix also contains first and second derivatives
% verbose - 0: no plots; 1: plot gamm and spline spproximation
% Outputs
% p = control points (3xN) array
% knots = uniform knots

M = zeros(N,N);
b = zeros(3*N,1);
for n = 1:length(xg)
    B = colmat(3*(n-1)+1,:);
    M = M + B'*B*J_(n+1)*wg(n); % index n+1 for gamm_ and J because the 
                                % first element is for s=0 (start point)
    b = b + kron(B',gamm_(:,n+1))*J_(n+1)*wg(n);
end
M = kron(M, eye(3));
C = zeros(6,3*N);
C(1:3,1:3) = eye(3);
C(4:6,3*N-3+(1:3)) = eye(3);
A = [M C'; C zeros(6)];
c = [b; gamm_(:,1); gamm_(:,end)];

x = A\c;
p = x(1:3*N);
err = norm(M*p - b);
p = reshape(p, 3, N);
