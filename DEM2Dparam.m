%% Parameters %%
function par = DEM2Dparam()
    par = struct('N',[], 'g',[],'mu', [],'r',[],'bBox',[],'dt',[],'T',[],'step',[],'kN',[],'dN',[]);

    %number of particles
    par.N = 2;

    % gravity
    par.g = 0.81; %[m/s^2]
    % friction coefficient mu \in [0,1)
    par.mu = 0.3;
    par.muWall = 0.3;

    % mean radius 
    par.r = [0.9 0.9]; %[m]
    % bounding box, x-length, z-length (height)
        par.bBox = [ -2 -2; % x first comp z first comp
                 2 2]; % x second comp, z second comp
    par.spawnBox = [ -2 -2; % x first comp z first comp
                 2 2];
    % contact detection
    par.collisionThreshold = 1.25;
    % numerical simulation
    par.simulationStart = 0;
    par.simulationEnd = 4;
    par.dt = 1e-4;%1e-6
    par.T = round(par.simulationEnd/par.dt); %integrationSteps %1e4; 1e6; %2e5
   
    par.step = round(0.025/par.dt);
    par.VisualizationStep = par.step;
    par.CollisionTime = 5e-4;
    par.CollisionStep = round(par.CollisionTime/par.dt);
    %% force parameters

    par.rho = 27;%2700; % kg/m³
    par.Emodul = 1e8; % should be 1e8
    par.kN = par.Emodul*pi/2*par.r(1); % [N/m] stiffness
   % par.kT = 1/1.2*par.kN; % adjust accordingly
    par.dampN = 0.2; % correct? 2 % of critial damping
    par.dampT = 0.1; % tangential damping
    par.dampTwall = 0.1;
   % par.wallDistr = 0.1; % coefficient on wall
    % particle wall
    par.cohesion = 0;

    % 2 DOF or 3 DOF? Not fully implemented - carefull
    par.considerRotations = true;
    par.Cr = 0.5; % rolling resistance coefficient
    
    %% video parameters
    par.writePdf = false;
    par.writeEps = false;
    par.writePng = false;
    par.writeVid = false;
    par.videoname = 'video-40-rot';%video4-merged';
    par.video_framerate = 20;
    par.videoFontsize = 16;
    
    
    %% merge parameters
    par.merge = false;
    par.mergeThreashold = 10^-3; %Threashold for relative velocity to initialize merg

end