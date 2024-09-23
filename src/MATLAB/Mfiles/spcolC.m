function colmat = spcolC(knots, k, colpt)

% This is a dummy function, since MATLAB coder cannot generate for spcol
% This function is called instead by the generator code. 
% It is forced to never inline, so that the this function is explicitly
% called in the nested function Bishop
% Code generated for this function will be deleted. The function spcolC 
% implemented separately will actually compute the necessary
% For call from getBishopFrame, colpt is always a scalar, and two 
% derivatives are needed, so there are hardcoded
% This is the C prototype:
% void spcolC(const emxArray_real_T *knots, double k, double colpt,
%             emxArray_real_T *colmat)

coder.inline("never");

nknots = length(knots);
N = nknots - k;
colmat = colpt*rand(3,N);
