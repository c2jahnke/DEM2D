function [data,SuccessFlag] = DEM2Dinit(par)
 
    SuccessFlag = false;
    data = struct('position',[],'velocity',[],'acceleration',[],'radius',[],'mass',[],'deltaOld',[]);
    
    data.contactsWall = struct;
    data.contactsWall.isInitialized = zeros(par.N,4); % each particle may collide with all 4 walls
    data.contactsWall.contactAge = zeros(par.N,4); % count contact age for particles
    data.contactsWall.actuationPoint = zeros(par.N,2,4); % actuation point for each particle and each wall
    data.contactsWall.contactPoint = zeros(par.N,2,4);
    data.contactsWall.localContactPoint = zeros(par.N,2,4);
    data.contactsWall.rollingDeformation = zeros(par.N,2,4);
    data.contactsWall.accumulatedRollingDeformation = zeros(par.N,2,4);
    data.contactsWall.maxContactAge = 3;
    
    data.contactsParticle = struct;
    data.contactsParticle.isInitialized = zeros(par.N,par.N);
    data.contactsParticle.ActiveContactAge = zeros(par.N,par.N);
    data.contactsParticle.PassiveContactAge = zeros(par.N,par.N);
    data.contactsParticle.actuationPoint = zeros(2,par.N,par.N);
    data.contactsParticle.contactPoint = zeros(2,par.N,par.N);
    data.contactsParticle.maxContactAge = 3;
    data.contactsParticle.deactivated = zeros(par.N,1);
    data.contactsParticle.mergedParticles = false;
    data.contactsParticle.rollingDeformation = zeros(2,par.N,par.N);
    data.contactsParticle.accumulatedRollingDeformation = zeros(2,par.N,par.N);
    
    
    
    data.contactsGlued = zeros(par.N,par.N);
    
    data.contactsMerged = struct;
    data.contactsMerged.N = 0;
    data.contactsMerged.index = zeros(2,par.N);
    data.contactsMerged.position = zeros(4,par.N);
    data.contactsMerged.relativePosition = zeros(4,par.N);
    data.contactsMerged.positionMerged = zeros(2,par.N);
    data.contactsMerged.velocityMerged = zeros(2,par.N);
    data.contactsMerged.timeFlag = zeros(1,par.N);
    data.contactsMerged.timePoint = zeros(par.N,par.N);
    
    %     % radii randomly distributed

    c = ceil(par.N/2);%*rand(1,1)); %random number between 0 and 100
    C = zeros(par.N,1); C(1:c) = par.r(1); C(c:par.N) = par.r(2);
    data.radius = C;
    data.mass = par.rho*pi*(data.radius).^2;%4/3*(data.radius).^3*pi; % 3D, how to define 2D mass?
    isArranged = 0;
    count = 0;
    data.inertiaTensor = 0.5*data.mass.*data.radius;
    
    % simplified positioning as in Kleinert's Code for larger piles
%     [X,Y]=meshgrid(par.bBox(1)+par.r(1):2*par.r(1):par.bBox(2)-par.r(1),...
%                par.bBox(3)+par.r(1):2*par.r(1):par.bBox(4)-par.r(1));
%            if(min(size(X)) < par.N)
%                disp("Bounding Box too small for number of particles.")
%            end
%     N_=size(X,1)*size(X,2);
%     index=randperm(N_);
%     n=min(N_,par.N);
%     for j=1:n
%     i=index(j);
%    radius=rand*(rmax-rmin)+rmin;
%     mass = pi*radius^2*params.particle_density;
%     data.position(:,j) = [X(i);Y(i)];
%     end
    %% particle positioning for small piles
    while isArranged == 0
        count = count +1;
        disp("Arranging particles.")
        x = [par.spawnBox(1)+par.r(2)+0.1*max(data.radius):2.2*C:par.spawnBox(2)-par.r(2)-0.1*max(data.radius)];
        z = [par.spawnBox(3)+par.r(2)+0.1*max(data.radius):2.2*C:par.spawnBox(4)-par.r(2)-0.1*max(data.radius)];
        %[X,Z] = meshgrid(x,z);

        if(length(x) > length(z))  
            if(length(x) < sqrt(par.N))
                error('Reduce number of particles or increase bounding box.');
            end
        xInd = randperm(length(x),par.N);
        zInd = randi([1 length(z)],1,par.N);
        else
            if(length(z) < sqrt(par.N))
                error('Reduce number of particles, the particle radius or increase bounding box.');
            end 
        end
        maxSize = length(x)*length(z);
        [X,Z] = meshgrid(x,z); % generate mesh
        Zr = reshape(Z,1,maxSize);
        zInd = randperm(length(Zr),par.N);
        %xInd = randi([1 length(Xr)],1,par.N);
       
        % choose positions from indices
        data.position = [X(zInd)+ 0.1*max(data.radius)*(rand(1,par.N)-0.5) ; Z(zInd) + 0.1*max(data.radius)*(rand(1,par.N)-0.5)];

        d = DEM2Ddist(data.position(1,:),data.position(2,:));

%        Sum up radii:
        R = zeros(par.N,par.N);
        for i=1:par.N
            for j=1:par.N
                if i==j
                    R(i,j)=0;
                else
                    R(i,j) = data.radius(i) + data.radius(j) + 0.01;
                end
            end 
        end

        if max(max(d<abs(R)))==0
            isArranged = 1;
            disp('Successfully arranged.')
            SuccessFlag = true;
        elseif count > 1000
            disp("Too many particles for bounding box. 1000 attempts to arange them failed.")
            return;
            %break;
        end
    end
    data.velocity =  0.01*(rand(2,par.N)-0.05);
    data.delta = zeros(par.N,par.N);
    data.angular = zeros(2,par.N); % angular position, velocity and acceleration

    save('data')
end