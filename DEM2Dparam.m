%% Parameters %%
function par = DEM2Dparam()
par = struct('N',[], 'g',[],'mu', [],'r',[],'bBox',[],'dt',[],'T',[],'step',[],'kN',[],'dN',[]);

%number of particles
par.N = 3;

% gravity
par.g = 9.81; %[m/s²]
% friction coefficient mu \in [0,1)
par.mu = 0.5;

% mean radius 
par.r = [0.1 0.06]; %[m]
% bounding box, x-length, z-length (height)
par.bBox = [ 1.5 3.5; 
             2.5 4];

% numerical simulation
par.dt = 5e-4;%1e-6
par.T = 1e4;%1e6; %2e5
par.step = round(0.05/par.dt);

% force parameters
par.Emodul = 1e3; % should be 1e8
par.kN = par.Emodul*pi/2*par.r(1); % [N/m] stiffness
par.kT = 0.1*par.kN; % adjust accordingly
par.dampN = 1; % correct? 2 % of critial damping
par.dampT = 1; % tangential damping
par.wallDistr = 0.1; % coefficient on wall
% particle wall
end