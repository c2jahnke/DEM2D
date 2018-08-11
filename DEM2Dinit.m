function data = DEM2Dinit(par)

    data = struct('position',[],'velocity',[],'acceleration',[],'radius',[],'mass',[],'deltaOld',[]);
    m = 3/4*(par.r(1))^3*pi;
%     % radii and weight randomly distributed
%     C = rand(1,par.N);
%     data.radius = par.r(1) + (par.r(2)-par.r(1)).*C; % [m], 20 mm = .02 m 
%     data.mass = m + m/3 .*C; % kg
    % two different radii
    c = ceil(par.N*rand(1,1)); %random number between 0 and 100
    C = zeros(par.N,1); C(1:c) = par.r(1); C(c:par.N) = par.r(2);
    data.radius = C;
    data.mass = 3/4*(data.radius).^3*pi;
    cond = 0;
    count = 0;
    while cond == 0
        count = count +1;
        disp("Arranging particles.")

        data.position = [(par.bBox(2) - par.bBox(1) -2*par.r(1))*rand(1,par.N)+par.r(1) + par.bBox(1); (par.bBox(4)-par.bBox(3)-2*par.r(1))*rand(1,par.N)+par.r(1)+par.bBox(3)];

        d = DEM2Ddist(data.position(1,:),data.position(2,:));

        % Sum up radii:
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
            cond = 1;
            disp('Successfully arranged.')
        elseif count > 10
            disp("Too many particles for bounding box. 10 attempts to arange them failed.")
            break;
        end
    end
    data.velocity =  0.001*(rand(2,par.N)-0.5);
    data.acceleration = zeros(2,par.N);
    data.velocity(2,:) = 0.01*ones(1,par.N);
    data.deltaOld = zeros(par.N,par.N);
    data.delta = zeros(par.N,par.N);
    data.angular = zeros(3,par.N); data.angular(2,:) = -0.1*ones(1,par.N);% 2D: angle, angular velocity
    data.Xc = zeros(2,par.N,par.N); % contact point in local coordinates
    data.Xinorm = zeros(2,par.N,par.N);
    data.XinormOld = zeros(2,par.N,par.N);
    save('data')
end