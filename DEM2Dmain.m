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
R1  = zeros(T/step+1,par.N);
M1  = zeros(T/step+1,par.N);
P1(1,:,:) = data.position;
R1(1,:) = data.radius;
M1(1,:) = data.mass;
n1(1) = par.N;

% ---------------------------- Iteratrion ---------------------------- %
j = 1;
count = 0;
for k = 1:T
    count = count +1;
    [pk,vk] = DEM2Dsolve_expl(data,par);
    
    data.position = pk;
    data.velocity = vk;
    
    
    if count == step
        % Sistemo i contatori
        j = j+1;
        count = 0;
        
        % Definisco matrici per il disegno:
        P1(j,:,1:par.N) = pk;
        R1(j,1:par.N) = data.radius;
        M1(j,1:par.N) = data.mass;
        n1(j) = par.N;
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
    DEM2DplotDyn(data,par,P1,k);
    hold on
    title('Simulation')
    text(par.bBox(1)-0.3,par.bBox(3)-0.1,['t = ',num2str((k-1)*dt*step,'%10.2f')])
    axis([-0.23+par.bBox(1) par.bBox(2)+0.03 -0.23+par.bBox(3) par.bBox(4)+0.03])
    axis equal
    hold off 
    
    if(k ==250)
        pause;
    end

    F(j) = getframe();
    
    pause(dt)
    
end

