function [L, M1, M2] = getSplineCurveLengthMoments(P,knots,axisPt,axisDir)
% get length of spline curve, and its first and second moments about a
% specified axis
% INPUTS:
% P = 3xN array of control points 
% knots = knot vector
% axisPt, axisDir - specify the axis about which to compute moments; axisPt
%                   is a point the axis goes through, and sxisDir is the
%                   direction the axis points in, so that points on the
%                   axis are of the form axisPt + t*axisDir
% OUTPUTS:
% L = length
% M1 = first moment
% M2 = second moment
% NOTES: 
% * Do addpath .. to use gaqdm
% * This function is useful in computing swinging frequency as 
%   sqrt(g*M1/M2)/2/pi

N = size(P,2); % number of control points
m = length(knots); % number of knots

d = m - N - 1; % degree of spline

Ng = 3; % number of integration points per element

Nel = N-d; % number of elements (knot intervals)
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

xg_ = repmat(xg',2,1); % multiplicity 2 for each collocation point so as to
                       % get function and derivative
xg_ = xg_(:);
colmat = spcol(knots, d+1, xg_); % order d+1 (degree d) B-splines

L = 0;
M1 = 0;
M2 = 0;
for n = 1:length(xg)
    Bp = colmat(2*(n-1)+2,:)';
    xip = P*Bp;
    nxip = norm(xip);
    L = L + nxip*wg(n);
    
    if (nargout > 1) % want moments
        B = colmat(2*(n-1)+1,:)';
        xi = P*B;
        t = ((xi-axisPt)'*axisDir)/(axisDir'*axisDir); % closest pt on axis
        z = norm(xi-axisPt-t*axisDir);
        M1 = M1 + z*nxip*wg(n);
        M2 = M2 + z^2*nxip*wg(n);
    end
end
