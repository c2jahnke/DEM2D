%% Parameters %%
function par = DEM2Dparam()
par = struct('N',[], 'g',[],'mu', [],'r',[],'bBox',[],'dt',[],'T',[],'step',[],'kN',[],'dN',[]);

%number of particles
par.N = 2;

% gravity
par.g = 0.1*0.981; %[m/s²]
% friction coefficient mu \in [0,1)
par.mu = 0.5;

% mean radius 
par.r = [0.03 0.025]; %[m]
% bounding box, x-length, z-length (height)
par.bBox = [ 2.5 3.5; 
             3 4];

% numerical simulation
par.dt = 1e-4;
par.T = 2e5;
par.step = 0.05/par.dt;

% force parameters

par.kN = 0.2; % [N/m] stiffness
par.kT = 1e-1; % adjust accordingly
par.dampN = 1e-6; % correct? 2 % of critial damping
par.dampT = 1e-6; % tangential damping
par.wallDistr = 0.1; % coefficient on wall
end