function colmat = spcolC(knots, k, colpts, nderiv)

% coder.extrinsic('spcol');
% colmat = spcol(knots, k, repmat(colpts,nderiv,1));
nknots = length(knots);
N = nknots - k;
colmat = zeros(length(colpts)*nderiv,N);
