function colmat = spcolC(knots, k, colpts, nderiv)

% This is a dummy function, since MATLAB coder cannot generate for spcol
% This function is called instead by the generator code. 
% It is forced to never inline, so that the this function is explicitly
% called in the nested function Bishop
% Code generated for this function will be deleted. The function spcolC in
% a Spline library will actually compute the necessary

coder.inline("never");

nknots = length(knots);
N = nknots - k;
colmat = rand(length(colpts)*nderiv,N);
