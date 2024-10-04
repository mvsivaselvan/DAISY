% Bushing parameters
m = 2.41; % lb-s^2/in
II = 4500; % lb-s^2/in in^2
if 1 % cover plate not stiffened
    kv = 6950; % lb/in
    kr = 3245000; % lb-in/rad
    zetav = 0.006; % damping ratio for vertical motion
    zetar = 0.008; % damping ratio for rocking motion
else % cover plate stiffenned
    kv = 39900; % lb/in
    kr = 7530000; % lb-in/rad
    zetav = 0.004; % damping ratio for vertical motion
    zetar = 0.008; % damping ratio for rocking motion
end
cv = zetav*2*sqrt(kv*m);
cr = zetar*2*sqrt(kr*II);
h = 20; % in
g = 386.4; % in/s^2
thet0 = 20/180*pi; % 1 degree imperfection
if 1
    fprintf('vertical frequency = %gHz\n',sqrt(kv/m)/2/pi);
    fprintf('rocking frequency = %gHz\n',sqrt(kr/II)/2/pi);
    fprintf('frequency ratio: vertical/rocking = %g\n',...
            sqrt(kv/m)/sqrt(kr/II));
end

eqdata = load('eqdata_long.dat');
ux = -eqdata(:,1)*g;
uz = -eqdata(:,3)*g;
dtsample = 1/256;

coupling = @(t,x)(BushingDynamics(t,x,m,II,kv,kr,cv,cr,h,...
                                  g,thet0,ux,uz,dtsample));
nocoupling = @(t,x)(BushingRocking(t,x,m,II,kr,cr,h,g,thet0,...
                                   ux,uz,dtsample));

odeopt = odeset('reltol',1e-6,'abstol',1e-8','maxstep',1e-1);
optimopt = optimoptions('fsolve','display','off');

% initial conditions
delti = -m*g/kv;
theti = fsolve(@(thet_)(kr*thet_-m*g*h*sin(thet_+thet0)),0,optimopt);

tsolve = (0:length(ux)-1)'*dtsample;

fprintf('Solving w/ coupling ...\n')
[T1,X1] = ode45(coupling,tsolve,[delti;theti;0;0],odeopt);

fprintf('Solving w/o coupling ...\n')
[T2,X2] = ode45(nocoupling,tsolve,[theti;0],odeopt);

fprintf('Solving linearized model ...\n')
Mlin = [m -m*h*sin(theti+thet0); -m*h*sin(theti+thet0) II];
Klin = [kv 0; 0 kr - m*g*h*cos(theti+thet0)];
Clin = [cv 0; 0 cr];
Rlin = -[m 0; m*h*cos(theti+thet0) -m*h*sin(theti+thet0)];
A = [zeros(2) eye(2); -Mlin\Klin -Mlin\Clin];
B = [zeros(2); Mlin\Rlin];
C = eye(4);
D = [];
bushingLin = ss(A, B, C, D);
T3 = (0:length(ux)-1)'*dtsample;
X3 = lsim(bushingLin, [ux uz], T3);
X3 = interp1(T3, X3, T2);

fprintf('Max absolute value of theta in the coupled case = %g\n',...
        max(abs(X1(:,2)-theti)));

fprintf('amplification = %g\n',...
        max(abs(X1(:,2)-theti))/max(abs(X2(:,1)-theti)));

figure(100),
    plot(T1, (X1(:,2)-theti)*kr, T2, (X2(:,1)-theti)*kr, T2, X3(:,2)*kr),
    grid on,
    xlabel('Time (s)'),
    ylabel('Base moment (lb-in)'),
    hlegend = legend('w/ coupling','w/o coupling','linearized','location','northwest');
    set(gca, 'fontsize', 12, 'ylim', [-1.2 1.2]*1e5, 'ytick', [-1.2:0.2:1.2]*1e5)
    set(hlegend, 'fontsize', 12)
    
out = [T1 kr*X1(:,2) kr*X2(:,1) kr*(theti+X3(:,2))];
save 'twenty-lin.txt' out -ascii -double