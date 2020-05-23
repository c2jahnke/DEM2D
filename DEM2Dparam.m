%% Parameters %%
function par = DEM2Dparam()
    par = struct;
    par.software = 'MATLAB';%'GNU Octave';%'MATLAB';%'GNU Octave';
    %number of particles
    par.N = 2;

    % gravity
    par.g = 0;% -9.81; %[m/s^2]
    par.g_vert = 0;
    % friction coefficient mu \in [0,1)
    par.mu = 0.5;
    par.muWall = 0.3;

    % mean radius 
    par.r = [1 1.0]; %[m]
    % bounding box, x-length, z-length (height)
    par.bBox = [ -2 -2; % x first comp z first comp
                 2 2]; % x second comp, z second comp
%     par.bBox = [ -0.02 -0.02; % x first comp z first comp
%                      0.02 0.02]; % x second comp, z second comp
    %par.spawnBox = [ -2 -2; % x first comp z first comp
     %            2 2];
    % contact detection
    par.collisionThreshold = 1.1;
    % numerical simulation
    par.simulationStart = 0;
    par.simulationEnd = 1;
    par.dt = 1e-3;%1e-6
    par.T = round(par.simulationEnd/par.dt); %integrationSteps %1e4; 1e6; %2e5
    par.VisualResolution = 0.025;
    par.step = round(par.VisualResolution/par.dt);
    par.VisualizationStep = par.step;
    par.CollisionTime = 1e-3;
    par.CollisionStep = round(par.CollisionTime/par.dt);
    %% force parameters

    par.rho = 1/pi%;2700; % kg/m³
    par.Emodul = 1e5;%1e8; % should be 1e8
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
    par.videoname = 'video-cohesion';%video4-merged';
    par.video_framerate = 20;
    par.videoFontsize = 16;
    par.videoPartFontsize = 3;
    
    
    %% merge parameters
    par.merge = false;
    par.mergeThreashold = 10^-3; %Threashold for relative velocity to initialize merge

end