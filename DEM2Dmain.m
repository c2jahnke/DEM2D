clc, clear all, close all
% -------------------------- Initialization -------------------------- %
par = DEM2Dparam();
LoadData = false;
if(LoadData == true)
SuccessFlag = true;
data = DEM2Dload();
else
[data,SuccessFlag] = DEM2Dinit(par);
if(SuccessFlag == 0)
    return
end
end

% -------------------------- Plot initial state -------------------------- %
figure(1)
DEM2Dplot(data,par);

dt = par.dt;
T = par.T;
VisualizationStep = par.step; CollisionStep =  par.step;
par.VisualizationStep = VisualizationStep;
P1  = zeros(T/VisualizationStep+1,2,par.N); P1(1,:,:) = data.position;
V1  = zeros(T/VisualizationStep+1,2,par.N); V1(1,:,:) = data.velocity;
A1 = zeros(T/VisualizationStep+1,3,par.N); A1(1,:,:) = data.angular;
% n1(1) = par.N;

% ---------------------------- Iteration ---------------------------- %

j = 1; VisCounter = 0; ColCounter = 0;
c = DEM2Dcontacts(data,par);
for k = 1:T
    VisCounter = VisCounter +1;
    ColCounter = ColCounter +1;
    if ColCounter == CollisionStep
        ColCounter = 0;
        c = DEM2Dcontacts(data,par);
    end
     [pk,vk,ak,data] = DEM2Dsolve_expl(data,par);
 %   [pk,vk,ak,data] = DEM2Dsolve_pgs(data,par,c.contacts);
    data.position = pk;
    data.velocity = vk;
    data.angular = ak;
    
    if VisCounter == VisualizationStep
        j = j+1; VisCounter = 0;
        P1(j,:,1:par.N) = data.position; 
        A1(j,:,1:par.N) = data.angular; 
        V1(j,:,1:par.N) = data.velocity;
    end
    
end
figure(1)
% ----------------------------- Drawing ----------------------------- %
DEM2DplotSim(P1,A1,par,data,j)

%save('data')