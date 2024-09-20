function colmat = spcolC(knots, nknots, k, colpts, ncolpts, nderiv)

% This is a dummy function, since MATLAB coder cannot generate for spcol
% This function is called instead by the generator code. 
% It is forced to never inline, so that the this function is explicitly
% called in the nested function Bishop
% Code generated for this function will be deleted. The function spcolC in
% a Spline library will actually compute the necessary
% This is the C prototype:
% spcolC(double* knots, int nknots, int k, double* colpts, int ncolpts, 
%        int nderiv, double* colmat, int ncolmat);

coder.inline("never");

nknots = length(knots);
N = nknots - k;
colmat = rand(length(colpts)*nderiv,N);
