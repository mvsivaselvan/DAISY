function [x, J] = myCircle(s, R, thet1, thet2)
% Parametric circular arc: as s goes from 0 to 1, sweeps a circular arc of
% radius R centered at the origin in the X-Z plane from angle thet1 to
% angle thet2
% Inputs
% s = vector of parameters to evaluate the circular arc
% R = radius of the arc
% thet1, thet2 = start and end angles
% Outputs
% x = (X, Y, Z) coordinates of points on the arc at parameters s
%     (3xlength(s) matrix)
% J = jacobian, ||x'(s)||
% To test, run
% [x, J] = myCircle(linspace(0,1,101),2,-3*pi/4,-pi/4);
% figure, plot3(x(1,:),x(2,:),x(3,:))
% axis equal
% view([0 1 0])

% first generate the arc from 0 to (thet2-thet1), then rotate by this much
c1 = cos(thet1);
s1 = sin(thet1);
rot = [c1 0 -s1; 0 1 0; s1 0 c1]; 

thet = s(:)'*(thet2-thet1); % s(:)' so that we get a row vector
x = rot*[R*cos(thet); zeros(1,length(s)); R*sin(thet)];
% in the following rot doesn't affect norm, so drop
xp = [-R*sin(thet); zeros(1,length(s)); R*cos(thet)]*(thet2-thet1); 
J = sqrt(sum(xp.^2,1));
