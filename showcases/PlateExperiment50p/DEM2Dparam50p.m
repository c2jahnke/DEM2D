%% Parameters %%
function par = DEM2Dparam()
    par = struct;
    par.software = 'MATLAB';%'GNU Octave';
    
    %% Choose contact algorithm
    par.algorithm = 'LIN';%'PGJ';%'PBD';'PGJ';'HMD';'LIN'; 'PBD-LIN'
    %   % use PGJ (non-smooth) Projected Gauﬂ Jacobi scheme (prototype)
    %   % use PGS (non-smooth) Projected Gauﬂ Seidel scheme (prototype)
    %   % use Position Based Dynamics (Mueller & Macklin et all) (prototype)
    %   % use Hertz-Mindlin and Deresievicz DEM contact model (prototype)
    %   % use PDM (PBD and DEM around tool) hybrid approach (so far not properly implemented prototype)
    
    % only for linear DEM 'LIN'
    par.Frozen = 0; % frozen particles 
    % number of particles
    par.N = 10;
    % load particle sample
    par.dataStr = 'data_50p.mat'
    % gravity
    par.g =  -9.81;% -9.81; %[m/s^2]
    par.g_vert = 0;
    % friction coefficient mu \in [0,1)
    par.mu = 0.3;
    par.muWall = 0.01;

    % mean radius 
    par.r = [0.3 0.42]; %[m]
    % bounding box, x-length, z-length (height)
    par.bBox = [ -2 -2; % x first comp z first comp
                 2 2]; % x second comp, z second comp
%     par.bBox = [ -0.02 -0.02; % x first comp z first comp
%                      0.02 0.02]; % x second comp, z second comp
    par.spawnBox = [ -1 -2; % x first comp z first comp
                 1 1]; % x second comp, z second comp
 
    %% Tool in DEM simulation            
    par.toolBool = 1;
    par.toolbBox = [ 2.050 -1.18; % x first comp z first comp
                 2.150 -0.1];
    par.toolSpeed = [-0.03;0.005];
    % contact detection
    par.collisionThreshold = 1.1;
    % numerical simulation
    par.simulationStart = 0;
    par.simulationEnd = 2;
    par.dt = 1e-3;%1e-6
    par.T = round(par.simulationEnd/par.dt); %integrationSteps %1e4; 1e6; %2e5
    par.CollisionTime = 1e-3; % collision detection, must coincide with par.dt for prototypes, for linear DEM it can be larger
    par.CollisionStep = round(par.CollisionTime/par.dt);
    
    par.VisualResolution = 0.05 %[s]  visualization
    par.step = round(par.VisualResolution/par.dt);
    par.VisualizationStep = par.step;

    %% model parameters
    par.rho = 2700%;2700; % kg/m¬≥
    par.Emodul = 1e8;%1e8; % should be 1e8
    par.dampN = 0.2; % correct? 2 % of critial damping
    par.dampT = 0.02; % tangential damping
    par.dampTwall = 0.02;
    
    % particle particle cohesion
    par.cohesion = 0;
    %% PGJ / PGS parameters
    par.w = 0.2; % relaxation factor, normally w = 0.2;
    par.maxIter = 1;
    %% PBD parameters
    par.gamma = 8;
    par.nSteps = 1;
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
    par.writePdf = false;
    par.writeEps = false;
    par.writePng = false;
    par.writeFig = false;
    par.writeVid = true;
    par.videotitle = [num2str(par.N) ' part ' par.algorithm ' E = ' num2str(par.Emodul/(1e6)) ' MPa dt =' num2str(par.dt) ' s'];
    par.videoname = [num2str(par.N) '_part_' par.algorithm '_E_' num2str(par.Emodul/(1e6)) '_MPa_mu=' num2str(par.mu) ];% frozen DEM';%video4-merged';
    if( par.algorithm(1)== 'P')
        if(par.algorithm(2) == 'G')
            par.videoname = [par.videoname '_maxIter_' num2str(par.maxIter)];
        else
            par.videoname = [par.videoname '_gamma_' num2str(par.gamma)];
        end
    end
    par.video_framerate = 20;
    par.videoFontsize = 16;
    par.videoPartFontsize = 8;
end