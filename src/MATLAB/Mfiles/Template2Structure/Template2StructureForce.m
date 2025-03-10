function [F, K, C, M] ...
         = Template2StructureForce(d_S1, phi_S1, ...
                        dd_S1, phid_S1, ddd_S1, phidd_S1, ...
                        m_S1, II_S1, KT_S1, CT_S1, KR_S1, CR_S1, r_S1, ...
                        RJ_S1, ...
                        d_S2, phi_S2, ...
                        dd_S2, phid_S2, ddd_S2, phidd_S2, ...
                        m_S2, II_S2, KT_S2, CT_S2, KR_S2, CR_S2, r_S2, ...
                        RJ_S2, ...
                        gamma1, gamma2, ...
                        Pmid, varThetamid, ...
                        P0, ...
                        gamma1dot, gamma2dot, ...
                        Pmiddot, varThetamiddot, ...
                        gamma1ddot, gamma2ddot, ...
                        Pmidddot, varThetamidddot, ...
                        x01_C, RE1_C, r1_C, ...
                        x02_C, RE2_C, r2_C, ...
                        R0, II_C, ...
                        rho, EA, EI, GJ, betAX, betBEND, betTOR, ...
                        sg, wg, nel, ...
                        colmat, colmat_brev, colmat_bar, ...
                        d, dbrev, dbar, ...
                        Mbar, ...
                        u, ...
                        Kbar11, Dbar11, dynamic, ...
                        maskS1Dof, maskS2Dof, maskCDof, ...
                        indexS1Dof, indexS2Dof, indexCDof) 

% Inputs
% d_S1 ... RJ_S1 = inputs for RigidBodyForce (structure 1)
% d_S2 ... RJ_S2 = inputs for RigidBodyForce (structure 2)
% gamma1 ... dynamic = inputs for CableForceRotBCinCoord (conductor)
% maskS1Dof ... indexCDof = assembly parameters
% Outputs
% F = unbalanced force (3*N-12+Nbrev+JADof)x1 - force and moment
% K = stiffness matrix (3*N-12+Nbrev+JADof)x(3*N-12+Nbrev+JADof)
% C = damping matrix - derivative of force w.r.t. velocity components
% M = mass matrix - derivative of force w.r.t. acceleration components

% number of basis functions (or control points)
N = size(colmat,2);
Nbrev = size(colmat_brev,2);

JADof = length(indexS1Dof) + length(indexS2Dof); %JADof: Joints Active Dof
F = zeros(3*N-12+Nbrev+JADof,1);
K = zeros(3*N-12+Nbrev+JADof);

[F_S1, K_S1] = ...
    RigidBodyForce(d_S1, phi_S1, ...
                   dd_S1, phid_S1, ddd_S1, phidd_S1, ...
                   m_S1, II_S1, KT_S1, CT_S1, KR_S1, CR_S1, r_S1, ...
                   RJ_S1, u);
F(indexS1Dof) = F(indexS1Dof) + F_S1(maskS1Dof); 
[F_S2, K_S2] = ...
    RigidBodyForce(d_S2, phi_S2, ...
                   dd_S2, phid_S2, ddd_S2, phidd_S2, ...
                   m_S2, II_S2, KT_S2, CT_S2, KR_S2, CR_S2, r_S2, ...
                   RJ_S2, u);
F(indexS2Dof) = F(indexS2Dof) + F_S2(maskS2Dof); 
[F_C, K_C] ...
    = CableForceRotBCinCoord(d_S1, phi_S1, ...
                             gamma1, ...
                             d_S2, phi_S2, ...
                             gamma2, ...
                             Pmid, varThetamid, ...
                             P0, ...
                             dd_S1, phid_S1, gamma1dot, ...
                             dd_S2, phid_S2, gamma2dot, ...
                             Pmiddot, varThetamiddot, ...
                             ddd_S1, phidd_S1, gamma1ddot, ...
                             ddd_S2, phidd_S2, gamma2ddot, ...
                             Pmidddot, varThetamidddot, ...
                             x01_C, RJ_S1, RE1_C, r1_C, ...
                             x02_C, RJ_S2, RE2_C, r2_C, ...
                             R0, II_C, ...
                             rho, EA, EI, GJ, betAX, betBEND, betTOR, ...
                             sg, wg, nel, ...
                             colmat, colmat_brev, colmat_bar, ...
                             d, dbrev, dbar, ...
                             Mbar, ...
                             u, ...
                             Kbar11, Dbar11, dynamic);
F(indexCDof) = F(indexCDof) + F_C(maskCDof); 

if nargout>1
    K(indexS1Dof,indexS1Dof) = K(indexS1Dof,indexS1Dof) + ...
                            K_S1(maskS1Dof,maskS1Dof);
    K(indexS2Dof,indexS2Dof) = K(indexS2Dof,indexS2Dof) + ...
                            K_S2(maskS2Dof,maskS2Dof);
    K(indexCDof,indexCDof) = K(indexCDof,indexCDof) + ...
                            K_C(maskCDof,maskCDof);
end

if dynamic
    C = zeros(3*N-12+Nbrev+JADof);
    M = zeros(3*N-12+Nbrev+JADof);
    [~, ~, C_S1, M_S1] = ...
        RigidBodyForce(d_S1, phi_S1, ...
                       dd_S1, phid_S1, ddd_S1, phidd_S1, ...
                       m_S1, II_S1, KT_S1, CT_S1, KR_S1, CR_S1, r_S1, ...
                       RJ_S1, u);
    C(indexS1Dof,indexS1Dof) = C(indexS1Dof,indexS1Dof) + ...
                            C_S1(maskS1Dof,maskS1Dof);
    M(indexS1Dof,indexS1Dof) = M(indexS1Dof,indexS1Dof) + ...
                            M_S1(maskS1Dof,maskS1Dof);
    [~, ~, C_S2, M_S2] = ...
        RigidBodyForce(d_S2, phi_S2, ...
                   dd_S2, phid_S2, ddd_S2, phidd_S2, ...
                   m_S2, II_S2, KT_S2, CT_S2, KR_S2, CR_S2, r_S2, ...
                   RJ_S2, u);
    C(indexS2Dof,indexS2Dof) = C(indexS2Dof,indexS2Dof) + ...
                            C_S2(maskS2Dof,maskS2Dof);
    M(indexS2Dof,indexS2Dof) = M(indexS2Dof,indexS2Dof) + ...
                            M_S2(maskS2Dof,maskS2Dof);
    [~, ~, C_C, M_C] ...
        = CableForceRotBCinCoord(d_S1, phi_S1, ...
                             gamma1, ...
                             d_S2, phi_S2, ...
                             gamma2, ...
                             Pmid, varThetamid, ...
                             P0, ...
                             dd_S1, phid_S1, gamma1dot, ...
                             dd_S2, phid_S2, gamma2dot, ...
                             Pmiddot, varThetamiddot, ...
                             ddd_S1, phidd_S1, gamma1ddot, ...
                             ddd_S2, phidd_S2, gamma2ddot, ...
                             Pmidddot, varThetamidddot, ...
                             x01_C, RJ_S1, RE1_C, r1_C, ...
                             x02_C, RJ_S2, RE2_C, r2_C, ...
                             R0, II_C, ...
                             rho, EA, EI, GJ, betAX, betBEND, betTOR, ...
                             sg, wg, nel, ...
                             colmat, colmat_brev, colmat_bar, ...
                             d, dbrev, dbar, ...
                             Mbar, ...
                             u, ...
                             Kbar11, Dbar11, dynamic);
    C(indexCDof,indexCDof) = C(indexCDof,indexCDof) + ...
                            C_C(maskCDof,maskCDof);
    M(indexCDof,indexCDof) = M(indexCDof,indexCDof) + ...
                            M_C(maskCDof,maskCDof);
end