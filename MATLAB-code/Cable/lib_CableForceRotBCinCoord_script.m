% LIB_CABLEFORCEROTBCINCOORD_SCRIPT   Generate static library
%  CableForceRotBCinCoord from CableForceRotBCinCoord, CableMbar, getBishopFrame, RigidBodyForce,
%  SplineApproximation.
% 
% Script generated from project 'lib_CableForceRotBCinCoord.prj' on 19-Sep-2024.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.EmbeddedCodeConfig'.
cfg = coder.config('lib','ecoder',false);
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = false;
cfg.GenCodeOnly = true;
cfg.SupportNonFinite = false;

%% Define argument types for entry-point 'CableForceRotBCinCoord'.
ARGS = cell(5,1);
ARGS{1} = cell(57,1);
ARGS{1}{1} = coder.typeof(0,[3 1]);
ARGS{1}{2} = coder.typeof(0,[3 1]);
ARGS{1}{3} = coder.typeof(0);
ARGS{1}{4} = coder.typeof(0,[3 1]);
ARGS{1}{5} = coder.typeof(0,[3 1]);
ARGS{1}{6} = coder.typeof(0);
ARGS{1}{7} = coder.typeof(0,[3 Inf],[0 1]);
ARGS{1}{8} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9} = coder.typeof(0,[3 Inf],[0 1]);
ARGS{1}{10} = coder.typeof(0,[3 1]);
ARGS{1}{11} = coder.typeof(0,[3 1]);
ARGS{1}{12} = coder.typeof(0);
ARGS{1}{13} = coder.typeof(0,[3 1]);
ARGS{1}{14} = coder.typeof(0,[3 1]);
ARGS{1}{15} = coder.typeof(0);
ARGS{1}{16} = coder.typeof(0,[3 Inf],[0 1]);
ARGS{1}{17} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{18} = coder.typeof(0,[3 1]);
ARGS{1}{19} = coder.typeof(0,[3 1]);
ARGS{1}{20} = coder.typeof(0);
ARGS{1}{21} = coder.typeof(0,[3 1]);
ARGS{1}{22} = coder.typeof(0,[3 1]);
ARGS{1}{23} = coder.typeof(0);
ARGS{1}{24} = coder.typeof(0,[3 Inf],[0 1]);
ARGS{1}{25} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{26} = coder.typeof(0,[3 1]);
ARGS{1}{27} = coder.typeof(0,[3 3]);
ARGS{1}{28} = coder.typeof(0,[3 3]);
ARGS{1}{29} = coder.typeof(0,[3 1]);
ARGS{1}{30} = coder.typeof(0,[3 1]);
ARGS{1}{31} = coder.typeof(0,[3 3]);
ARGS{1}{32} = coder.typeof(0,[3 3]);
ARGS{1}{33} = coder.typeof(0,[3 1]);
ARGS{1}{34} = coder.typeof(0,[3 Inf],[0 1]);
ARGS{1}{35} = coder.typeof(0,[3 3]);
ARGS{1}{36} = coder.typeof(0);
ARGS{1}{37} = coder.typeof(0);
ARGS{1}{38} = coder.typeof(0);
ARGS{1}{39} = coder.typeof(0);
ARGS{1}{40} = coder.typeof(0);
ARGS{1}{41} = coder.typeof(0);
ARGS{1}{42} = coder.typeof(0);
ARGS{1}{43} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{44} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{45} = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{46} = coder.typeof(0,[Inf Inf],[1 1]);
ARGS{1}{47} = coder.typeof(0,[Inf Inf],[1 1]);
ARGS{1}{48} = coder.typeof(0,[Inf Inf],[1 1]);
ARGS{1}{49} = coder.typeof(0);
ARGS{1}{50} = coder.typeof(0);
ARGS{1}{51} = coder.typeof(0);
ARGS{1}{52} = coder.typeof(0,[Inf Inf],[1 1]);
ARGS{1}{53} = coder.typeof(0,[3 1]);
ARGS{1}{54} = coder.typeof(0,[Inf Inf],[1 1]);
ARGS{1}{55} = coder.typeof(0,[Inf Inf],[1 1]);
ARGS{1}{56} = coder.typeof(0);
ARGS{1}{57} = coder.typeof(0);

%% Define argument types for entry-point 'CableMbar'.
ARGS{2} = cell(6,1);
ARGS{2}{1} = coder.typeof(0,[3 Inf],[0 1]);
ARGS{2}{2} = coder.typeof(0);
ARGS{2}{3} = coder.typeof(0);
ARGS{2}{4} = coder.typeof(0,[Inf Inf],[1 1]);
ARGS{2}{5} = coder.typeof(0,[Inf Inf],[1 1]);
ARGS{2}{6} = coder.typeof(0,[Inf  1],[1 0]);

%% Define argument types for entry-point 'getBishopFrame'.
ARGS{3} = cell(4,1);
ARGS{3}{1} = coder.typeof(0,[3 Inf],[0 1]);
ARGS{3}{2} = coder.typeof(0,[1 Inf],[0 1]);
ARGS{3}{3} = coder.typeof(0);
ARGS{3}{4} = coder.typeof(0,[Inf  1],[1 0]);

%% Define argument types for entry-point 'RigidBodyForce'.
ARGS{4} = cell(15,1);
ARGS{4}{1} = coder.typeof(0,[3 1]);
ARGS{4}{2} = coder.typeof(0,[3 1]);
ARGS{4}{3} = coder.typeof(0,[3 1]);
ARGS{4}{4} = coder.typeof(0,[3 1]);
ARGS{4}{5} = coder.typeof(0,[3 1]);
ARGS{4}{6} = coder.typeof(0,[3 1]);
ARGS{4}{7} = coder.typeof(0);
ARGS{4}{8} = coder.typeof(0,[3 3]);
ARGS{4}{9} = coder.typeof(0,[3 3]);
ARGS{4}{10} = coder.typeof(0,[3 3]);
ARGS{4}{11} = coder.typeof(0,[3 3]);
ARGS{4}{12} = coder.typeof(0,[3 3]);
ARGS{4}{13} = coder.typeof(0,[3 1]);
ARGS{4}{14} = coder.typeof(0,[3 3]);
ARGS{4}{15} = coder.typeof(0,[3 1]);

%% Define argument types for entry-point 'SplineApproximation'.
ARGS{5} = cell(6,1);
ARGS{5}{1} = coder.typeof(0,[Inf  3],[1 0]);
ARGS{5}{2} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{5}{3} = coder.typeof(0);
ARGS{5}{4} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{5}{5} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{5}{6} = coder.typeof(0,[Inf Inf],[1 1]);

%% Invoke MATLAB Coder.
codegen -config cfg -package MATLABelements CableForceRotBCinCoord -args ARGS{1} CableMbar -args ARGS{2} getBishopFrame -args ARGS{3} RigidBodyForce -args ARGS{4} SplineApproximation -args ARGS{5}