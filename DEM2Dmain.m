clc, clear all
% -------------------------- Initialization -------------------------- %
par = DEM2Dparam();
data = DEM2Dinit(par);

DEM2Dplot(data,par);
%pause
dt = par.dt;
T = par.T;
step = par.step;
P1  = zeros(T/step+1,2,par.N);
V1  = zeros(T/step+1,2,par.N);
A1 = zeros(T/step+1,3,par.N);
R1  = zeros(T/step+1,par.N);
M1  = zeros(T/step+1,par.N);
P1(1,:,:) = data.position;
V1(1,:,:) = data.velocity;
A1(1,:,:) = data.angular;
R1(1,:) = data.radius;
M1(1,:) = data.mass;
n1(1) = par.N;

% ---------------------------- Iteratrion ---------------------------- %
j = 1;
count = 0;
for k = 1:T
    count = count +1;
    [pk,vk,ak] = DEM2Dsolve_expl(data,par);
    
    data.position = pk;
    data.velocity = vk;
    data.angular = ak;
    
    if count == step
        j = j+1;
        count = 0;
        
        P1(j,:,1:par.N) = pk;
        A1(j,:,1:par.N) = data.angular;
    end
end

pause

% ----------------------------- Drawing ----------------------------- %
theta = linspace(0,2*pi,15);
s = sin(theta);
c = cos(theta); 

J = j;

for k=1:J
    figure (2)
    DEM2DplotDyn(data,par,P1,A1,k);
    hold on
    title('Simulation')
    text(par.bBox(1)-0.1,par.bBox(3)-0.1,['t = ',num2str((k-1)*dt*step,'%10.2f')])
    axis([-0.23+par.bBox(1) par.bBox(2)+0.03 -0.23+par.bBox(3) par.bBox(4)+0.03])
    axis equal
    hold off 
    
    F(j) = getframe();
    
    pause(dt)
    
end
save P1
