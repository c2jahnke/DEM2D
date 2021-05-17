%% Parameters for the simulation%%
function par = DEM2Dparam()
    par = struct;
    par.software = 'MATLAB';%'GNU Octave';
    
    %% Choose contact algorithm
    par.algorithm = 'LIN';%'PGJ';%'PBD';'PGJ';'HMD';'LIN'; 'PBD-LIN'
    %   % use PGJ (non-smooth) Projected Gau� Jacobi scheme (prototype)
    %   % use PGS (non-smooth) Projected Gau� Seidel scheme (prototype)
    %   % use Position Based Dynamics (Mueller & Macklin et all) (prototype)
    %   % use Hertz-Mindlin and Deresievicz DEM contact model (prototype)
    %   % use PDM (PBD and DEM around tool) hybrid approach (so far not properly implemented prototype)
    
    % only for linear DEM 'LIN'
    par.Frozen = 1; % frozen particles 
    % number of particles
    par.N = 50;
    % gravity
    par.g = -9.81;% -9.81;%[m/s^2]
    par.g_vert = 0; 
    % friction coefficient mu \in [0,\infty)
    par.mu = 0.5; %[]
    par.muWall = 0.0; %[]

    % particle radius radius 
    par.r = [0.3 0.3]; %[m]
    %% bounding Box for particle container
    % x-length, z-length (height)
    par.bBox = [ -2 -2; % [m] x first comp z first comp
                 2 2]; % x second comp, z second comp
    %par.spawnBox = [ -2 -2; % x first comp z first comp
     %            2 2];
    %% Tool in DEM simulation
    par.toolBool = 1; % works only for linear penalty-based DEM
    par.toolbBox = [2.050 -0.808; % x first comp z first comp
                 2.10 -0.02];
    par.toolSpeed = [-0.02;0.001];
    % contact detection
    par.collisionThreshold = 1.25;
    %% numerical time stepping
    par.simulationStart = 0;% [s] 
    par.simulationEnd = 5.0;% [s] 
    par.dt = 1e-4;% [s] 1e-6
    par.T = round(par.simulationEnd/par.dt); %integrationSteps %1e4; 1e6; %2e5
    par.CollisionTime = 1e-4; % collision detection, must coincide with par.dt for prototypes, for linear DEM it can be larger
    par.CollisionStep = round(par.CollisionTime/par.dt);
    
    par.VisualResolution = 0.05 %[s]  visualization
    par.step = round(par.VisualResolution/par.dt);
    par.VisualizationStep = par.step;

    %% model parameters
    par.rho = 2700;% % kg/m³
    par.Emodul = 1e8;% N/m^2
    par.dampN = 0.2; % 2 % of critial damping
    par.dampT = 0.02; % tangential damping
    par.dampTwall = 0.02;
    % particle particle cohesion
    par.cohesion =0;
    %% Hertz-Mindlin Deresievicz model parameters
    par.nu = 0.3; % Poisson-ratio, Querkontraktionszahl
    par.e = 0.5; % coefficient of restitution
    par.G = par.Emodul/(2*(1+par.nu)) ; % shear modulus for isotromic materials: 2G(1+nu)=E
    
    %% Rotations for linear DEM 2 DOF or 3 DOF? 
    par.considerRotations = false;
    par.Cr = 0.990; % rolling resistance coefficient
    par.CrWall = 0.990;% rolling resistance coefficient for walls
    
    %% Merge particles for linear DEM parameters
    par.merge = false;
    par.mergeThreashold = 10^-3; % Threashold for relative velocity to initialize merge

    %% video parameters
    par.writePdf = true;
    par.writeEps = false;
    par.writePng = false;
    par.writeVid = true;
    par.videotitle = [num2str(par.N) ' part ' par.algorithm ' E = ' num2str(par.Emodul/(1e6)) ' MPa dt =' num2str(par.dt) ' s'];
    par.videoname = [num2str(par.N) '_part_' par.algorithm '_E_' num2str(par.Emodul/(1e6)) '_MPa'];% frozen DEM';%video4-merged';
    par.video_framerate = 20;
    par.videoFontsize = 16;
    par.videoPartFontsize = 8;
end