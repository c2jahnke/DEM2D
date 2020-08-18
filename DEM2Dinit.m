function [data,par,SuccessFlag] = DEM2Dinit(par)
 
    SuccessFlag = false;
    data = struct('position',[],'velocity',[],'angular',[],'acceleration',[],'radius',[],'mass',[],'deltaOld',[]);
    
    % simplified positioning as in Kleinert's Code for larger piles
    [X,Y]=meshgrid(par.bBox(1)+par.r(2):2*par.r(2):par.bBox(2)-par.r(2),...
               par.bBox(3)+par.r(2):2*par.r(2):par.bBox(4)-par.r(2));
    U = X + 10^-1* par.r(1)*rand(size(X));
    V = Y + 10^-1* par.r(1)*rand(size(X));
    N_=size(X,1)*size(X,2);
    
    index=randperm(N_); rmax = par.r(2); rmin = par.r(1);
    n=min(N_,par.N);
    if(n < par.N)
       disp(["Bounding Box too small for number of particles.Reducing numer of particles to " num2str(n)])
       par.N = n;
    end
    for j=1:n
        i=index(j);
        data.radius(j)=rand*(rmax-rmin)+rmin;
        data.mass(j) = pi*data.radius(j)^2*par.rho;
        data.position(:,j) = [U(i);V(i)];
    end
    SuccessFlag = true;

    data.velocity =  0.01*(rand(2,par.N)-0.05);
    data.delta = sparse(par.N,par.N);
    data.angular = sparse(2,par.N); % angular position, velocity and acceleration
    
    data.toolbBox = par.toolbBox;
    if(par.PGJ)
        par.dt = par.CollisionStep;
    end
    data.contactsWall = struct;
    data.contactsWall.isInitialized = sparse(par.N,4); % each particle may collide with all 4 walls
    data.contactsWall.contactAge = zeros(par.N,4); % count contact age for particles
    data.contactsWall.actuationPoint = zeros(par.N,2,4); % actuation point for each particle and each wall
    data.contactsWall.globalContactPoint = zeros(par.N,2,4);
    data.contactsWall.localContactPoint = zeros(par.N,2,4);
    data.contactsWall.globalContactPoint2 = zeros(par.N,2,4);
    data.contactsWall.localContactPoint2 = zeros(par.N,2,4);
    data.contactsWall.rollingDeformation = zeros(par.N,2,4);
    data.contactsWall.accumulatedRollingDeformation = zeros(par.N,2,4);
    data.contactsWall.maxContactAge = 3;
    
    
    data.contactsTool = struct;
    data.contactsTool.isInitialized = sparse(par.N,8); % each particle may collide with all 4 walls
    data.contactsTool.contactAge = zeros(par.N,8); % count contact age for particles
    data.contactsTool.actuationPoint = zeros(par.N,2,8); % actuation point for each particle and each wall
    data.contactsTool.globalContactPoint = zeros(par.N,2,8);
    data.contactsTool.localContactPoint = zeros(par.N,2,8);
    data.contactsTool.globalContactPoint2 = zeros(par.N,2,8);
    data.contactsTool.localContactPoint2 = zeros(par.N,2,8);
    data.contactsTool.rollingDeformation = zeros(par.N,2,8);
    data.contactsTool.accumulatedRollingDeformation = zeros(par.N,2,8);
    data.contactsTool.maxContactAge = 3;
    
    data.contactsParticle = struct;
    data.contactsParticle.isInitialized = sparse(par.N,par.N);
    data.contactsParticle.ActiveContactAge = sparse(par.N,par.N);
    data.contactsParticle.PassiveContactAge = sparse(par.N,par.N);
    data.contactsParticle.actuationPoint = zeros(2,par.N,par.N);
    data.contactsParticle.contactPoint = zeros(2,par.N,par.N);
    data.contactsParticle.maxContactAge = 3;
    data.contactsParticle.deactivated = zeros(par.N,1);
    data.contactsParticle.mergedParticles = false;
    data.contactsParticle.rollingDeformation = zeros(2,par.N,par.N);
    data.contactsParticle.accumulatedRollingDeformation = zeros(2,par.N,par.N);
    
    data.contactsGlued = sparse(par.N,par.N);
    
    data.contactsMerged = struct;
    data.contactsMerged.N = 0;
    data.contactsMerged.index = spalloc(ceil(par.N/2),par.N,par.N); 
    % index 1 (i): aggrete body index, index 2 (j): particles in body (i), third component: maximal size
    data.contactsMerged.aggregateSize = spalloc(ceil(par.N/2),1,ceil(par.N/2));
    data.contactsMerged.position = zeros(4,par.N);
    % data.contactsMerged.relativePosition = zeros(4,par.N);
    data.contactsMerged.relativePosition = spalloc(2*ceil(par.N/2),ceil(par.N/2),par.N);
    data.contactsMerged.positionMerged = zeros(2,par.N);
    data.contactsMerged.velocityMerged = zeros(2,par.N);
    data.contactsMerged.angularMerged = zeros(2,par.N);
    data.contactsMerged.inertiaTensor = spalloc(ceil(par.N/2),1,ceil(par.N/2));
    
    data.contactsMerged.timeFlag = zeros(1,par.N);
    data.contactsMerged.timePoint = sparse(par.N,par.N);
    save('data')
end