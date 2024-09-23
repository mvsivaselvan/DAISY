%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS SCRIPT IS USED FOR TESTING AND CODE GENERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Template 2 - Catenary --- TR bushing - surge arrester configuration
% Reference configuration = circular arc

% Spline setup
N = 10; % number of control points
d = 3; % cubic spline
dbrev = 2;
dbar = d - 2;
Ng = 3; % number of Gauss points per element
[sg, wg, nel, knots, colmat, colmat_brev, colmat_bar] ...
    = CableSplineSetup(N, d, dbrev, dbar, Ng);

% Reference geometry - cicular arc
span = 125; % in
slack = 2; % in
L = span + slack; % length
R = L/(pi/3); % want a 60-degree arc
curve_gamm = @(s)myCircle(s,R,-2*pi/3,-pi/3);
refGeom = curve_gamm(linspace(0,1,101));
%[gamm_1, J_1] = curve_gamm([0;sg;1]);
[gamm_, J_] = piecewiseLinearCurve(refGeom',[0;sg;1]);
[P0, ~] = SplineApproximation(gamm_, J_, N, sg, wg, colmat);
P0 = P0 - P0(:,1); % translate start end to origin
R0 = getBishopFrame(P0, knots, d, sg);

if (~coder.target('MATLAB')) % generating code, so spcol not called in
                             % getBishopFrame, and loading R0 from a
                             % previous run
    load R0;
end

% Material properties
g = 386.4; % in/s^2
dirg = [0 0 -1];
Acable = 1.57; % in^2
density = 2.54e-4; % (lb-s^2/in)/in^3
Ecable = 1e7; % psi
Icable = 2.18e-3; % in^4
rho = density*Acable;
EI = Ecable*Icable;
EA = Ecable*Acable;
GJ = (2*EI)/(2*(1+0.3)); % I is Imin, i.e., wire-wise, so J = 2*I;
                         % and G = E/2/(1+nu), taking n = 0.3
betBEND = 0.01; % damping coefficient for bending
betAX = 0.01; % damping coefficient for axial
betTOR = 0.01; % damping coefficient for torsion

II = density*(Acable/pi)^2/64*diag([1 1 2]);

% Mbar, Kbar11, Dbar11
[Mbar, Kbar11, Dbar11] ...
     = CableMbar(P0, EA, betAX, colmat, colmat_bar, wg);

% undeformed state
x01 = P0(:,1); r1 = zeros(3,1); 
d10 = zeros(3,1);
Rb10 = R0(:,1:3);
I3 = eye(3);
RJ1 = I3; RE1 = I3;
phi10hat = logm(RJ1'*Rb10*RE1');
phi10 = [phi10hat(3,2); phi10hat(1,3); phi10hat(2,1)];
gamma10 = norm(P0(:,2)-P0(:,1));
x02 = P0(:,end); r2 = zeros(3,1);
d20 = zeros(3,1);
Rb20 = R0(:,(length(sg)+1)*3+(1:3));
RJ2 = I3; RE2 = I3;
phi20hat = logm(RJ2'*Rb20*RE2');
phi20 = [phi20hat(3,2); phi20hat(1,3); phi20hat(2,1)];
gamma20 = -norm(P0(:,end)-P0(:,end-1));
Pmid0 = P0(:,3:end-2);
varTheta0 = zeros(size(colmat_brev,2),1);
varThetamid0 = varTheta0(2:end-1);
d1dot = zeros(3,1);
phi1dot = zeros(3,1);
gamma1dot = 0;
d2dot = zeros(3,1);
phi2dot = zeros(3,1);
gamma2dot = 0;
Pmiddot = zeros(size(Pmid0));
varThetamiddot = zeros(size(varThetamid0));
d1ddot = zeros(3,1);
phi1ddot = zeros(3,1);
gamma1ddot = 0;
d2ddot = zeros(3,1);
phi2ddot = zeros(3,1);
gamma2ddot = 0;
Pmidddot = zeros(size(Pmid0));
varThetamidddot = zeros(size(varThetamid0));
[uXt, uYt, uZt] = deal(0,0,0);
u = g*dirg(:) + [uXt; uYt; uZt]; 
dynamic = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try one step of static analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% line search parameters (see page 33 of Nocedal and Wright)
armijo_gamm = 0.9;
armijo_c1 = 1e-2;
% Rotate the two ends to horizontal (Rb = Identity)
phi1 = [0; 0; 0];
phi2 = [0; 0; 0];
% Horizontal displacement of end 2
d21 = [span-P0(1,end); 0; 0];
x = [gamma10; gamma20; Pmid0(:); varThetamid0];
for k = 1:200 % Newton iterations
%     if k == 35
%         keyboard
%     end
    gamma1k = x(1);
    gamma2k = x(2);
    Pmidk = reshape(x(3:2+3*(N-4)),3,N-4);
    varThetamidk = x(3+3*(N-4):end);
    [Fb, Kb] ...
      = CableForceRotBCinCoord(d10, phi1, gamma1k, ...
                               d21, phi2, gamma2k, ...
                               Pmidk, varThetamidk, ...
                               P0, ...
                               d1dot, phi1dot, gamma1dot, ...
                               d2dot, phi2dot, gamma2dot, ...
                               Pmiddot, varThetamiddot, ...
                               d1ddot, phi1ddot, gamma1ddot, ...
                               d2ddot, phi2ddot, gamma2ddot, ...
                               Pmidddot, varThetamidddot, ...
                               x01, RJ1, RE1, r1, ...
                               x02, RJ2, RE2, r2, ...
                               R0, II, ...
                               rho, EA, EI, GJ, betAX, betBEND, betTOR, ...
                               sg, wg, nel, ...
                               colmat, colmat_brev, colmat_bar, ...
                               d, dbrev, dbar, ...
                               Mbar, ...
                               u, ...
                               Kbar11, Dbar11, dynamic, 0);
    gk = [Fb(7); Fb(14:end)];
    if (norm(gk) <= 1e-8)
        break;
    end
    D1g = Kb([7 14:end], [7 14:end]);
    Dx = -D1g\gk;
    armijo_alph = 1;
    for l = 1:50
        x_ = x + armijo_alph*Dx;
        gamma1_ = x_(1);
        gamma2_ = x_(2);
        Pmid_ = reshape(x_(3:2+3*(N-4)),3,N-4);
        varThetamid_ = x_(3+3*(N-4):end);
        Fb_ = CableForceRotBCinCoord(d10, phi1, gamma1_, ...
                               d21, phi2, gamma2_, ...
                               Pmid_, varThetamid_, ...
                               P0, ...
                               d1dot, phi1dot, gamma1dot, ...
                               d2dot, phi2dot, gamma2dot, ...
                               Pmiddot, varThetamiddot, ...
                               d1ddot, phi1ddot, gamma1ddot, ...
                               d2ddot, phi2ddot, gamma2ddot, ...
                               Pmidddot, varThetamidddot, ...
                               x01, RJ1, RE1, r1, ...
                               x02, RJ2, RE2, r2, ...
                               R0, II, ...
                               rho, EA, EI, GJ, betAX, betBEND, betTOR, ...
                               sg, wg, nel, ...
                               colmat, colmat_brev, colmat_bar, ...
                               d, dbrev, dbar, ...
                               Mbar, ...
                               u, ...
                               Kbar11, Dbar11, dynamic, 0);
        gk_ = [Fb_(7); Fb_(14:end)];
%         if (k == 35)
%             l, norm(gk_)
%         end
        if (0.5*norm(gk_)^2 <= (0.5 - armijo_alph*armijo_c1)*norm(gk)^2)
            break;
        end
        armijo_alph = armijo_alph*armijo_gamm;
    end
    x = x_;
    fprintf('k = %d, ||g|| = %g, ||Dx|| = %g, armijo_alph = %g\n', ...
             k, norm(gk), norm(Dx), armijo_alph);
end

eta = zeros(3,1);
rho1 = zeros(3,1);
rho2 = zeros(3,1);
qbar1 = [d10; phi1; gamma1k];
qbar2 = [d21; phi2; gamma2k];
q1 = CableBCTransinCoord(qbar1, ...
                         x01, RJ1, RE1, r1, Rb10, ...
                         eta, rho1,rho2);
q2 = CableBCTransinCoord(qbar2, ...
                         x02, RJ2, RE2, r2, Rb20, ...
                         eta, rho1,rho2);
P = [q1(1:3) q1(4:6) Pmidk q2(4:6) q2(1:3)];

bb0 = spmak(knots,P0);
bb = spmak(knots,P);
figure(101),
    fnplt(bb0),
    hold on,
    fnplt(bb),
    grid on,
    hold off,
    axis equal,
    view([0,-1,0]),
    xlabel('X position (in)')
    ylabel('Z position (in)')
    legend('Reference','Step 1')

%% Rigid body force
d__ = rand(3,1);
phi__ = rand(3,1)*pi;
dd__ = rand(3,1);
phid__ = rand(3,1)*pi;
ddd__ = rand(3,1);
phidd__ = rand(3,1)*pi;

m__ = rand(1);
II__ = rand(3); II = (II + II');
KT__ = zeros(3);
CT__ = zeros(3);
KR__ = zeros(3);
CR__ = zeros(3);
rr__ = rand(3,1);
RJ__ = randRot();
u__ = rand(3,1);

[F0, K0, C0, M0] = RigidBodyForce(d__,phi__,dd__,phid__,ddd__,phidd__, ...
                                  m__,II__,KT__,CT__,KR__,CR__,rr__, ...
                                  RJ__,u__);
return
%%
% derivative test
[Fb, Kb, Cb, Mb, Kb_FD, Cb_FD, Mb_FD] ...
       = CableRotBCinCoordStiffnessCheckFD(d10, phi1, gamma1k, ...
                               d21, phi2, gamma2k, ...
                               Pmidk, varThetamidk, ...
                               P0, ...
                               d1dot, phi1dot, gamma1dot, ...
                               d2dot, phi2dot, gamma2dot, ...
                               Pmiddot, varThetamiddot, ...
                               d1ddot, phi1ddot, gamma1ddot, ...
                               d2ddot, phi2ddot, gamma2ddot, ...
                               Pmidddot, varThetamidddot, ...
                               x01, RJ1, RE1, r1, ...
                               x02, RJ2, RE2, r2, ...
                               R0, II, ...
                               rho, EA, EI, GJ, betAX, betBEND, betTOR, ...
                               sg, wg, nel, ...
                               colmat, colmat_brev, colmat_bar, ...
                               d, dbrev, dbar, ...
                               Mbar, ...
                               u, ...
                               Kbar11, Dbar11, dynamic);
