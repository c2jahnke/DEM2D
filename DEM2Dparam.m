%% Parameters %%
function par = DEM2Dparam()
    par = struct('N',[], 'g',[],'mu', [],'r',[],'bBox',[],'dt',[],'T',[],'step',[],'kN',[],'dN',[]);

    %number of particles
    par.N = 2;

    % gravity
    par.g = 9.81; %[m/s²]
    % friction coefficient mu \in [0,1)
    par.mu = 0.3;
    par.muWall = 0.51;

    % mean radius 
    par.r = [0.125 0.128]; %[m]
    % bounding box, x-length, z-length (height)
        par.bBox = [ 2.5 1.0; % x first comp z first comp
                 3.5 3.0]; % x second comp, z second comp
    par.spawnBox = [ 2.5 1.0; % x first comp z first comp
                 3.5 3.0]
%     par.bBox = [ 0.5 0.0; % x first comp z first comp
%                  3.5 4.0]; % x second comp, z second comp
%     par.spawnBox = [ 1.0 0.0; % x first comp z first comp
%                  3.0 3.4]
    % contact detection
    par.collisionThreshold = 2;
    % numerical simulation
    par.simulationStart = 0;
    par.simulationEnd = 4;
    par.dt = 1e-5%1e-6
    par.T = round(par.simulationEnd/par.dt); %integrationSteps %1e4; 1e6; %2e5
   
    par.step = round(0.025/par.dt);
    par.VisualizationStep = par.step;
    par.CollisionTime = 5e-5;
    par.CollisionStep = round(par.CollisionTime/par.dt);
    %% force parameters

    par.rho = 27;%2700; % kg/m³
    par.Emodul = 1e6; % should be 1e8
    par.kN = par.Emodul*pi/2*par.r(1); % [N/m] stiffness
    par.kT = 1/1.2*par.kN; % adjust accordingly
    par.dampN = 0.2; % correct? 2 % of critial damping
    par.dampT = 0.01; % tangential damping
   % par.wallDistr = 0.1; % coefficient on wall
    % particle wall
    par.cohesion = 0;

    % 2 DOF or 3 DOF? Not fully implemented - carefull
    par.considerRotations = true;
    
    %% video parameters
    par.writePdf = false;
    par.writeEps = false;
    par.writePng = false;
    par.writeVid = false;
    par.videoname = 'video';%video4-merged';
    par.video_framerate = 20;
    par.videoFontsize = 16;
    
    
    %% merge parameters
    par.merge = false;
    par.mergeThreashold = 10^-3; %Threashold for relative velocity to initialize merg

end