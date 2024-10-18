function [qbar1, qbar2, P, varTheta, ...
          qbar1dot, qbar2dot, Pdot, varThetadot, ...
          qbar1ddot, qbar2ddot, Pddot, varThetaddot, ...
          Fc] = ...
            Template2StructureDynamicsOutputs(T, X, Xp, ...
                                                RJ_SA, RJ_B, ...
                                                x01_C, RE1_C, r1_C, ...
                                                x02_C, RE2_C, r2_C, ...
                                                Rb10, Rb20, ...
                                                N, Nbrev, ...
                                                maskSADof, maskBDof, ...
                                                indexSADof, indexBDof, ...
                                                P0, R0, II_C, ...
                                                rho, EA, EI, GJ, ...
                                                betAX, betBEND, betTOR, ...
                                                sg, wg, nel, ...
                                                colmat, colmat_brev, ...
                                                colmat_bar, ...
                                                d, dbrev, dbar, ...
                                                Mbar, ...
                                                u, ...
                                                Kbar11, Dbar11)

% Provides the displacement and velocity outputs from dynamic analysis 
% results (T,X) obtained from ode15i
% INPUTS
% T - time
% X, Xp = solutions obtained from ode15i and deval
% RJ_SA ... R0 = inputs for CableBCTransinCoord
% N = number of translational control points
% Nbrev = number of rotational control points
% maskSADof ... indexBDof = assembly parameters
% OUTPUTS
% qbar1 = joint 1 global disp coordinates (7xNsteps)
%         Nsteps: number of time steps (length of T)
% qbar2 = joint 2 global disp coordinates (7xNsteps)
% P = conductor translational control points coordinates (3NxNsteps)
% varTheta = conductor rotational control points coordinates (NbrevxNsteps)
% qbar1dot ... varThetadot = corresponding velocity components
% qbar1ddot ... varThetaddot = corresponding acceleration components
% Fc = conductor force output

% total number of Dofs
JADof = length(indexSADof) + length(indexBDof); %JADof: Joints Active Dof
NDof = 3*N-12+Nbrev+JADof;

% initialization of displacement and velocity for all time steps
Nsteps = length(T);

qbar1 = zeros(7,Nsteps);
qbar2 = zeros(7,Nsteps);
P = zeros(3*N,Nsteps);
varTheta = zeros(Nbrev,Nsteps);

if nargout > 4 % velocity components wanted
    qbar1dot = zeros(7,Nsteps);
    qbar2dot = zeros(7,Nsteps);
    Pdot = zeros(3*N,Nsteps);
    varThetadot = zeros(Nbrev,Nsteps);
end

if nargout > 8 % acceleration components wanted
    qbar1ddot = zeros(7,Nsteps);
    qbar2ddot = zeros(7,Nsteps);
    Pddot = zeros(3*N,Nsteps);
    varThetaddot = zeros(Nbrev,Nsteps);
end

if nargout > 12 % conducor force wanted
    Fc = zeros(6, Nsteps);
end

for k = 1:Nsteps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISPLACEMENT at each time step
    J_SA = zeros(6,1);
    J_SA(maskSADof)= X(k,indexSADof);
    J_B = zeros(6,1);
    J_B(maskBDof)= X(k,indexBDof);
    gamma1 = X(k,JADof+1);
    gamma2 = X(k,JADof+2);
    Pmid = reshape(X(k,JADof+3:JADof+2+3*(N-4)),3,N-4);
    varThetamid = X(k,JADof+2+3*(N-4)+(1:Nbrev-2))';
    qbar1(:,k) = [J_SA; gamma1];
    qbar2(:,k) = [J_B; gamma2];
    
    % transform displacement into local conductor coordinates (q)
    [q1, J1] = CableBCTransinCoord(qbar1(:,k), ...
                            x01_C, RJ_SA, RE1_C, r1_C, Rb10);
    [q2, J2] = CableBCTransinCoord(qbar2(:,k), ...
                            x02_C, RJ_B, RE2_C, r2_C, Rb20);
    P_ = [q1(1:3) q1(4:6) Pmid q2(4:6) q2(1:3)];
    P(:,k) = P_(:);
    varTheta(:,k) = [q1(end); varThetamid; q2(end)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VELOCITY at each time step
    if nargout > 4 % velocity components wanted
        Jdot_SA = zeros(6,1);
        Jdot_SA(maskSADof) = X(k,NDof+indexSADof);
        Jdot_B = zeros(6,1);
        Jdot_B(maskBDof) = X(k,NDof+indexBDof);
        gamma1dot = X(k,NDof+JADof+1);
        gamma2dot = X(k,NDof+JADof+2);
        Pmiddot = reshape(X(k,NDof+JADof+3:NDof+JADof+2+3*(N-4)),3,N-4);
        varThetamiddot = X(k,NDof+JADof+2+3*(N-4)+(1:Nbrev-2))';
        qbar1dot(:,k) = [Jdot_SA; gamma1dot];
        qbar2dot(:,k) = [Jdot_B; gamma2dot];
        
        % transform velocity into local conductor coordinates (q)
        q1dot = J1*qbar1dot(:,k);
        q2dot = J2*qbar2dot(:,k);
        Pdot_ = [q1dot(1:3) q1dot(4:6) Pmiddot q2dot(4:6) q2dot(1:3)];
        Pdot(:,k) = Pdot_(:);
        varThetadot(:,k) = [q1dot(end); varThetamiddot; q2dot(end)];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ACCELERATION at each time step
    if nargout > 8 % acceleration components wanted
        Jddot_SA = zeros(6,1);
        Jddot_SA(maskSADof) = Xp(k,NDof+indexSADof);
        Jddot_B = zeros(6,1);
        Jddot_B(maskBDof) = Xp(k,NDof+indexBDof);
        gamma1ddot = Xp(k,NDof+JADof+1);
        gamma2ddot = Xp(k,NDof+JADof+2);
        Pmidddot = reshape(Xp(k,NDof+JADof+3:NDof+JADof+2+3*(N-4)),3,N-4);
        varThetamidddot = Xp(k,NDof+JADof+2+3*(N-4)+(1:Nbrev-2))';
        qbar1ddot(:,k) = [Jddot_SA; gamma1ddot];
        qbar2ddot(:,k) = [Jddot_B; gamma2ddot];
        
        [~, ~, ~, Qtilde1] = CableBCTransinCoord(qbar1(:,k), ...
                            x01_C, RJ_SA, RE1_C, r1_C, Rb10, ...
                            qbar1dot(:,k), qbar1dot(:,k));
        [~, ~, ~, Qtilde2] = CableBCTransinCoord(qbar2(:,k), ...
                            x02_C, RJ_B, RE2_C, r2_C, Rb20, ...
                            qbar2dot(:,k), qbar2dot(:,k));
        
        % transform velocity into local conductor coordinates (q)
        q1ddot = J1*qbar1ddot(:,k)+Qtilde1*qbar1dot(:,k);
        q2ddot = J2*qbar2ddot(:,k)+Qtilde2*qbar2dot(:,k);
        Pddot_= [q1ddot(1:3) q1ddot(4:6) Pmidddot q2ddot(4:6) q2ddot(1:3)];
        Pddot(:,k) = Pddot_(:);
        varThetaddot(:,k) = [q1ddot(end); varThetamidddot; q2ddot(end)];       
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Conductor force at each time step
    if nargout > 12 % acceleration components wanted
        Fc_ = ...
        CableForceRotBCinCoord(J_SA(1:3), J_SA(4:6), gamma1, ...
                               J_B(1:3), J_B(4:6), gamma2, ...
                               Pmid, varThetamid, ...
                               P0, ...
                               Jdot_SA(1:3), Jdot_SA(4:6), gamma1dot, ...
                               Jdot_B(1:3), Jdot_B(4:6), gamma2dot, ...
                               Pmiddot, varThetamiddot, ...
                               Jddot_SA(1:3), Jddot_SA(4:6), gamma1ddot,...
                               Jddot_B(1:3), Jddot_B(4:6), gamma2ddot, ...
                               Pmidddot, varThetamidddot, ...
                               x01_C, RJ_SA, RE1_C, r1_C, ...
                               x02_C, RJ_B, RE2_C, r2_C, ...
                               R0, II_C, ...
                               rho, EA, EI, GJ, betAX, betBEND, betTOR, ...
                               sg, wg, nel, ...
                               colmat, colmat_brev, colmat_bar, ...
                               d, dbrev, dbar, ...
                               Mbar, ...
                               u, ...
                               Kbar11, Dbar11, 1, 0);
        Fc(:,k) = Fc_([1:3 8:10]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end