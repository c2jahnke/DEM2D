%% Parameters %%
function par = DEM2Dparam()
par = struct('N',[], 'g',[],'mu', [],'r',[],'bBox',[],'dt',[],'T',[],'step',[],'kN',[],'dN',[]);

%number of particles
par.N = 8;

% gravity
par.g = 9.81; %[m/s²]
% friction coefficient mu \in [0,1)
par.mu = 0.50;

% mean radius 
par.r = [0.045 0.045]; %[m]
% bounding box, x-length, z-length (height)
par.bBox = [ 1.5 2; % x first comp z first comp
             2.0 2.5]; % x second comp, z second comp
% contact detection
par.collisionThreshold = 2;
% numerical simulation
par.dt = 5e-4;%1e-6
par.T = 1e4;%1e6; %2e5
par.step = round(0.05/par.dt);

% force parameters

par.rho = 0.2700; % kg/m³
par.Emodul = 1e4; % should be 1e8
par.kN = par.Emodul*pi/2*par.r(1); % [N/m] stiffness
par.kT = 1/1.2*par.kN; % adjust accordingly
par.dampN = 0.12; % correct? 2 % of critial damping
par.dampT = 0.1; % tangential damping
par.wallDistr = 0.1; % coefficient on wall
% particle wall
par.cohesion = 0;
end