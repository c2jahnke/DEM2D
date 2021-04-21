% Simulate particles in two dimensions using the discrete element method

% Particle demonstration code with different collision solvers

% Default: DEM modell with linear spring dashpot interaction (par.algortithm = 'LIN')

% For a first start run DEM2Dmain.m
% To change the number of particles or different contact parameters, etc open DEM2Dparam.m
% To simulate the same particle setting (initial condition), set LoadData = true in DEM2Dmain.m 

% Rotations in beta stage: (some bugs pending), to turn it off set par.considerRotations = false
% Frozen particles for linear penalty based DEM
% Merged particles for linear penalty based DEM in beta stage
% Tool interaction for linear penalty based DEM

% To change the algorithm use par.algorithm = 'LIN' etc.
% Position Based Dynamics (PBD) in beta stage
% nonsmooth Projected Gauss Jacobi (PGJ)in beta stage
% nonsmooth Projected Gauss Seidel (PGS) in beta stage
% penalty based Hertz Mindlin Deresievicz DEM (HMD) in beta state