function R = getRandomRotation

v = rand(3,1);
v = v/norm(v); % randomly generate an axis to rotate about
vhat = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
thet_ = (0:3)*pi/2;
thet = thet_(randi(4)); % randomly select an angle to rotate by
R = expm(vhat*thet);