%% Parameters %%
function par = DEM2Dparam()
par = struct('N',[], 'g',[],'mu', [],'r',[],'bBox',[],'dt',[],'T',[],'step',[],'kN',[],'dN',[]);

%number of particles
par.N = 5;

% gravity
par.g = 0.0981; %[m*kg/s²]
% friction coefficient mu \in [0,1)
par.mu = 0.5;

% mean radius 
par.r = [0.03 0.015]; %[m]
% bounding box, x-length, z-length (height)
par.bBox = [ 2.5 3.5; 
             3 4];

% numerical simulation
par.dt = 5e-5;
par.T = 5e5;
par.step = 0.05/par.dt;

% force parameters

par.kN = 0.5; % [N/m] stiffness
par.kT = 0.2; % adjust accordingly
par.dampN = 0.02; % correct? 2 % of critial damping
par.dampT = 0.02; % tangential damping
end