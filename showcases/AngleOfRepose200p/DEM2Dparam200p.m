%% Parameters %%
function par = DEM2Dparam()
    par = struct;
    par.software = 'MATLAB';%'GNU Octave';%'MATLAB';%'GNU Octave';
    par.PGJ = 0; % use PGJ (non-smooth) scheme or explicit solution
    par.PBD = 0; % use Position Based Dynamics (Müller & Macklin et all) (prototype)
    par.HMD = 0; % use Hertz-Mindlin and Deresievicz DEM contact model (not implemented prototype)
    par.Frozen = 0; % frozen particles 
   %number of particles
    par.N = 200;

    % gravity
    par.g = -9.81; %[m/s^2]
    par.g_vert = 0;
    % friction coefficient mu \in [0,1)
    par.mu = 0.5;
    par.muWall = 0.5;

    % mean radius 
    par.r = [0.027 0.047]; %[m]
    % bounding box, x-length, z-length (height)
    par.bBox = [ -1.5 -2; % x first comp z first comp
                 1.5 0]; % x second comp, z second comp
%     par.bBox = [ -0.02 -0.02; % x first comp z first comp
%                      0.02 0.02]; % x second comp, z second comp
    %par.spawnBox = [ -2 -2; % x first comp z first comp
     %            2 2];
    par.toolBool = 0;
    par.toolbBox = [ 2.050 -1.18; % x first comp z first comp
                 2.150 -0.1];
    par.toolSpeed = [-0.03;0.005];
    % contact detection
    par.collisionThreshold = 1.25;
    % numerical simulation
    par.simulationStart = 0;
    par.simulationEnd = 10;
    par.dt = 1e-4;%1e-6
    par.T = round(par.simulationEnd/par.dt); %integrationSteps %1e4; 1e6; %2e5
    par.VisualResolution = 0.05;
    par.step = round(par.VisualResolution/par.dt);
    par.VisualizationStep = par.step;
    par.CollisionTime = 5e-4;
    par.CollisionStep = round(par.CollisionTime/par.dt);
    %% force parameters

    par.rho = 2700;%2700; % kg/m³
    par.Emodul = 1e8;%1e8; % should be 1e8
    %par.kN = par.Emodul*pi/2*par.r(1); % [N/m] stiffness
    %     par.kT = 1e8; % adjust accordingly
    par.dampN = 0.2; % correct? 2 % of critial damping
    par.dampT = 0.02; % tangential damping
    par.dampTwall = 0.02;
    
    % particle particle cohesion
    par.cohesion = 0;

    % 2 DOF or 3 DOF? 
    par.considerRotations = false;
    par.Cr = 0.990; % rolling resistance coefficient
    par.CrWall = 0.990;% rolling resistance coefficient for walls
    
    %% video parameters
    par.writePdf = true;
    par.writeEps = false;
    par.writePng = true;
    par.writeVid = true;
    par.videoname = 'video-200p-Schuettwinkel-mu0-5-dt1e-4_10s';%video4-merged';
    par.video_framerate = 20;
    par.videoFontsize = 16;
    par.videoPartFontsize = 3;
    
    
    %% merge parameters
    par.merge = false;
    par.mergeThreashold = 10^-3; %Threashold for relative velocity to initialize merge

end