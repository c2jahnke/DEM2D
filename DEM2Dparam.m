%% Parameters for the simulation%%
function par = DEM2Dpar()
    par = struct;
    %% Choose contact algorithm)
    par.software = 'MATLAB';%'GNU Octave';%'MATLAB';%'GNU Octave';
    par.PGJ = 0; % use PGJ (non-smooth) scheme (prototype)
    par.PBD = 0; % use Position Based Dynamics (Müller & Macklin et all) (prototype)
    par.HMD = 0; % use Hertz-Mindlin and Deresievicz DEM contact model (not implemented prototype)
    par.Frozen = 0; % frozen particles 
    %number of particles
    par.N = 200;
    % gravity
    par.g = -50;% -9.81;%[m/s^2]
    par.g_vert = 0; 
    % friction coefficient mu \in [0,1)
    par.mu = 0.3;
    par.muWall = 0.3;

    % particle radius radius 
    par.r = [0.630 0.630]; %[m]
    %% bounding Box for particle container
    % x-length, z-length (height)
    par.bBox = [ -2 -2; % x first comp z first comp
                 2 0]; % x second comp, z second comp
    %par.spawnBox = [ -2 -2; % x first comp z first comp
     %            2 2];
    %% Tool in DEM simulation
    par.toolBool = 0; % works only for force-based DEM
    par.toolbBox = [2.050 -0.808; % x first comp z first comp
                 2.10 -0.02];
    par.toolSpeed = [-0.05;0.001];
    % contact detection
    par.collisionThreshold = 1.25;
    %% numerical time stepping
    par.simulationStart = 0;
    par.simulationEnd = 4;
    par.dt = 1e-3;%1e-6
    par.T = round(par.simulationEnd/par.dt); %integrationSteps %1e4; 1e6; %2e5
    par.VisualResolution = 0.1 % visualization
    par.step = round(par.VisualResolution/par.dt);
    par.VisualizationStep = par.step;
    par.CollisionTime = 1e-2; % collision detection, must coincide for prototypes
    par.CollisionStep = round(par.CollisionTime/par.dt);
    %% model parameters
    par.rho = 2700;%2700; % kg/m³
    par.Emodul = 1e8;%1e8; % should be 1e8
    par.dampN = 0.2; % correct? 2 % of critial damping
    par.dampT = 0.02; % tangential damping
    par.dampTwall = 0.02;
    % particle particle cohesion
    par.cohesion =0;
    %% Hertz-Mindlin Deresievicz model parameters
    par.nu = 0.3; % Poisson-ratio, Querkontraktionszahl
    par.e = 0.5; % coefficient of restitution
    par.G = 1; % shear modulus
    
    % 2 DOF or 3 DOF? 
    par.considerRotations = true;
    par.Cr = 0.990; % rolling resistance coefficient
    par.CrWall = 0.990;% rolling resistance coefficient for walls
    
    %% video parameters
    par.writePdf = false;
    par.writeEps = false;
    par.writePng = false;
    par.writeVid = false;
    par.videoname = '2 PBD';%video4-merged';
    par.video_framerate = 20;
    par.videoFontsize = 16;
    par.videoPartFontsize = 8;
    
    %% merge parameters
    par.merge = true;
    par.mergeThreashold = 10^-3; % Threashold for relative velocity to initialize merge

end