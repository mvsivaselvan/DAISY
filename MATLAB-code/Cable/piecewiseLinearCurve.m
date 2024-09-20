function [x, J] = piecewiseLinearCurve(coord, s)
% Piecewise linear curve connecting given points coord
% The parametrization s \in [0,1] \mapsto x is such that 
% x(0) = the first point in coord
% x(s) is such that the length of the curve from the first point 
%      to x(s) is equal to s*(total length)
% This implies x(1) = last point in coord
% J(s) is the jacobina at x(s), which is simply L!

% Length of piecewise linear curve
dL = sqrt(diff(coord(:,1)).^2 + diff(coord(:,2)).^2 + diff(coord(:,3)).^2);
L = [0; cumsum(dL)]; % so that L(n) = distance from the first point 
                     % to the nth point

x = interp1(L, coord, s*L(end))'; % transpose, so that it is columnwise
J = ones(1,size(x,2))*L(end);
