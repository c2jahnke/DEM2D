clc, clear all
% -------------------------- Initialization -------------------------- %
par = DEM2Dparam();
%data = DEM2Dinit(par);
data = DEM2Dload();
DEM2Dplot(data,par);
%pause
dt = par.dt;
T = par.T;
VisualizationStep = par.step;
P1  = zeros(T/VisualizationStep+1,2,par.N); P1(1,:,:) = data.position;
V1  = zeros(T/VisualizationStep+1,2,par.N); V1(1,:,:) = data.velocity;
A1 = zeros(T/VisualizationStep+1,3,par.N); A1(1,:,:) = data.angular;
R1  = zeros(T/VisualizationStep+1,par.N); R1(1,:) = data.radius;
M1  = zeros(T/VisualizationStep+1,par.N); M1(1,:) = data.mass;
n1(1) = par.N;

pause()
% ---------------------------- Iteratrion ---------------------------- %
% TODO: locate forgotten semi-cola
j = 1;
counter = 0;
for k = 1:T
    counter = counter +1;
    [pk,vk,ak] = DEM2Dsolve_expl(data,par);
    
    data.position = pk;
    data.velocity = vk;
    data.angular = ak;
    
    if counter == VisualizationStep
        j = j+1;
        counter = 0; 
        
        P1(j,:,1:par.N) = pk; 
        A1(j,:,1:par.N) = data.angular;
    end
end

pause
figure
% ----------------------------- Drawing ----------------------------- %
DEM2DplotSim(P1,A1,par,data,j)

save P1
