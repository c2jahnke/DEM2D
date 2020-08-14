%% Parameters %%
function par = DEM2Dparam()
    par = struct;
    par.software = 'MATLAB';%'GNU Octave';%'MATLAB';%'GNU Octave';
    par.PGJ = 0; % use PGJ (non-smooth) scheme or explicit solution
    par.PBD = 0; % use Position Based Dynamics (Müller & Macklin et all)
    %number of particles
    par.N = 1;

    % gravity
    par.g =  0;% -9.81; %[m/s^2]
    par.g_vert = -9.81;
    % friction coefficient mu \in [0,1)
    par.mu = 0.3;
    par.muWall = 0.01;

    % mean radius 
    par.r = [0.5 0.5]; %[m]
    % bounding box, x-length, z-length (height)
    par.bBox = [ -2 -2; % x first comp z first comp
                 2 2]; % x second comp, z second comp
%     par.bBox = [ -0.02 -0.02; % x first comp z first comp
%                      0.02 0.02]; % x second comp, z second comp
    par.spawnBox = [ -1 -2; % x first comp z first comp
                 1 1]; % x second comp, z second comp
             
    par.toolBool = 1;
    par.toolbBox = [ -1.150 -1.0; % x first comp z first comp
                 -1.050 1.2];
    % contact detection
    par.collisionThreshold = 1.1;
    % numerical simulation
    par.simulationStart = 0;
    par.simulationEnd = 2;
    par.dt = 1e-3;%1e-6
    par.T = round(par.simulationEnd/par.dt); %integrationSteps %1e4; 1e6; %2e5
    par.VisualResolution = 0.025;
    par.step = round(par.VisualResolution/par.dt);
    par.VisualizationStep = par.step;
    par.CollisionTime = 1e-3;
    par.CollisionStep = round(par.CollisionTime/par.dt);
    %% force parameters

    par.rho = 2700%;2700; % kg/m³
    par.Emodul = 1e8;%1e8; % should be 1e8
    %par.kN = par.Emodul*pi/2*par.r(1); % [N/m] stiffness
    %     par.kT = 1e8; % adjust accordingly
    par.dampN = 0.2; % correct? 2 % of critial damping
    par.dampT = 0.02; % tangential damping
    par.dampTwall = 0.02;
    
    % particle particle cohesion
    par.cohesion = 0;

    % 2 DOF or 3 DOF? 
    par.considerRotations = true;
    par.Cr = 0.990; % rolling resistance coefficient
    par.CrWall = 0.990;% rolling resistance coefficient for walls
    
    %% video parameters
    par.writePdf = false;
    par.writeEps = false;
    par.writePng = false;
    par.writeVid = false;
    par.videoname = 'Simulation';
    par.video_framerate = 20;
    par.videoFontsize = 16;
    par.videoPartFontsize = 10;
    
    
    %% merge parameters
    par.merge = false;
    par.mergeThreashold = 10^-3; %Threashold for relative velocity to initialize merge

end