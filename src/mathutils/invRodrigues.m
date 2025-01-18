function phi = invRodrigues(R)

% R is a rotation matrix
% phi is a vector such that expm(hat(phi)) = R

% Computation follows: https://math.stackexchange.com/a/4044019/688462

if (abs(trace(R)-3) < 1e-10) % R is the identity
    phi = zeros(3,1);
    return
end

% find eigenvector of R corresponding to eigenvalue 1
X = (eye(2,2)-R(1:2,1:2))\R(1:2,3);
x = [X; 1];
x = x/norm(x);

thet = atan2(-trace(hat(x)*R),trace(R)-1);

phi = thet*x;
