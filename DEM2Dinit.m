function [data,SuccessFlag] = DEM2Dinit(par)
 
    SuccessFlag = false;
    data = struct('position',[],'velocity',[],'acceleration',[],'radius',[],'mass',[],'deltaOld',[]);
    data.contactsWall = struct;
    data.contactsWall.isInitialized = zeros(par.N,4); % each particle may collide with all 4 walls
    data.contactsWall.contactAge = zeros(par.N,4); % count contact age for particles
    data.contactsWall.actuationPoint = zeros(par.N,2,4); % actuation point for each particle and each wall
    data.contactsWall.maxContactAge = 3;
%     % radii randomly distributed

    c = ceil(par.N/2);%*rand(1,1)); %random number between 0 and 100
    C = zeros(par.N,1); C(1:c) = par.r(1); C(c:par.N) = par.r(2);
    data.radius = C
    data.mass = par.rho*pi*(data.radius).^2;%4/3*(data.radius).^3*pi; % 3D, how to define 2D mass?
    isArranged = 0;
    count = 0;
    
    
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
        x = [par.bBox(1)+par.r(2):1.2*C:par.bBox(2)-par.r(2)];
        z = [par.bBox(3)+par.r(2):1.2*C:par.bBox(4)-par.r(2)];
        %[X,Z] = meshgrid(x,z);
        if(length(x) < par.N)
            error('Reduce number of particles or increase bounding box.');
        end
        xInd = randperm(length(x),par.N);
        zInd = randperm(length(z),par.N);
        % choose positions from indices
        data.position = [x(xInd)+10^(-5)*rand(1,par.N) ; z(zInd)+10^(-5)*rand(1,par.N)];
        %data.position = [(par.bBox(2) - par.bBox(1) -2*par.r(1))*rand(1,par.N)+par.r(1) + par.bBox(1); (par.bBox(4)-par.bBox(3)-2*par.r(1))*rand(1,par.N)+par.r(1)+par.bBox(3)];

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
        elseif count > 10
            disp("Too many particles for bounding box. 10 attempts to arange them failed.")
            return;
            %break;
        end
    end
    data.velocity =  0.01*(rand(2,par.N)-0.05);
    data.acceleration = zeros(2,par.N);
    data.velocity(2,:) = 0.01*ones(1,par.N);
    data.deltaOld = zeros(par.N,par.N);
    data.delta = zeros(par.N,par.N);
    data.angular = zeros(3,par.N); % angular position, velocity and acceleration
    %data.angular(2,:) = 0.05*ones(1,par.N);% 2D: angle, angular velocity
    data.Xc = zeros(2,par.N,par.N); % contact point in local coordinates
    data.Xinorm = zeros(2,par.N,par.N);
    data.XinormOld = zeros(2,par.N,par.N);

    save('data')
end