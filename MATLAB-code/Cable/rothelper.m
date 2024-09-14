function a = rothelper(phi, flag)

% Computes helper functions associated with rotation matrix
% phi = exponential coordinates of rotation
% if flag == 0, compute a(1:2)
% if flag == 1, compute a(1:4,9)
% if flag == 2, compute a(1:4,9:10)
% if flag == 3, compute a(1:11)

a = zeros(11,1);

nphi = norm(phi);
c = cos(nphi);
s = sin(nphi);

if (nphi > 1e-5)
    nphi2 = nphi^2;
    if flag >= 0
        a(1) = s/nphi;
        a(2) = (1-c)/nphi2;
    end
    if flag >= 1
        nphi3 = nphi*nphi2;
        nphi4 = nphi*nphi3;
        a(3) = (nphi*c-s)/nphi3;
        a(4) = (nphi*s-2*(1-c))/nphi4;
        a(9) = (nphi - s)/nphi3;
    end
    if flag >= 2
        nphi5 = nphi*nphi4;
        a(10) = (3*s - nphi*(c+2))/nphi5;
    end
    if flag >= 3
        nphi6 = nphi*nphi5;
        nphi7 = nphi*nphi6;
        a(5) = -(3*nphi*c+(-3+nphi^2)*s)/nphi5;
        a(6) = (8+(-8+nphi^2)*c-5*nphi*s)/nphi6;
        a(11) = (8*nphi+7*nphi*c+(-15+nphi^2)*s)/nphi7;
    end
else
    if flag >= 0
        a(1) = polyval([1/120, 0, -(1/6), 0, 1],nphi);
        a(2) = polyval([1/720, 0, -(1/24), 0, 1/2],nphi);
    end
    if flag >= 1
        a(3) = polyval([-(1/840), 0, 1/30, 0, -(1/3)],nphi);
        a(4) = polyval([-(1/6720), 0, 1/180, 0, -(1/12)],nphi);
        a(9) = polyval([1/5040, 0, -(1/120), 0, 1/6],nphi);
    end
    if flag >= 2
        a(10) = polyval([-(1/60480), 0, 1/1260, 0, -(1/60)],nphi);
    end
    if flag >= 3
        a(5) = polyval([1/7560, 0, -(1/210), 0, 1/15],nphi);
        a(6) = polyval([1/75600, 0, -(1/1680), 0, 1/90],nphi);
        a(11) = polyval([1/831600, 0, -(1/15120), 0, 1/630],nphi);
    end    
end