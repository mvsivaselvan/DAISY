function R = randRot()

% Generate a random rotation matrix

phi = rand(3,1)*(2*pi);

a = rothelper(phi,1);

hphi = hat(phi);

R = eye(3) + a(1)*hphi + a(2)*hphi^2;
