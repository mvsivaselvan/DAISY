function qdot = Bishop(x, q, P, knots, d)

if (coder.target('MATLAB')) % running in MATLAB
    colmat = spcol(knots, d+1, [x;x;x]);
else % generating code
    colmat = spcolC(knots, d+1, x);
end

% B = colmat(1,:)'; % not used
Bp = colmat(2,:)';
Bpp = colmat(3,:)';

xip = P*Bp;
xipp = P*Bpp;

nxip = sqrt(xip'*xip);
tau = xip/nxip;

PP = eye(3) - tau*tau';
taup = (PP*xipp)/nxip;

omega = mycross(tau, taup);
Omega = [0 -omega'; omega hat(omega)];

k = 0.1; % feedback gain to prevent drift of q from unit norm

qTq = q'*q;
qdot = 0.5*Omega*q + k/2/qTq*(1 - qTq)*q;
