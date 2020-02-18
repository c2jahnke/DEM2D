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

% ------------------------ Plot initial state ------------------------ %
DEM2Dplot(data,par);

dt = par.dt;
T = par.T;
VisualizationStep = par.VisualizationStep; CollisionStep =  par.CollisionStep;



P1 = zeros(T/VisualizationStep+1,2,par.N); P1(1,:,:) = data.position;
V1 = zeros(T/VisualizationStep+1,2,par.N); V1(1,:,:) = data.velocity;
A1 = zeros(T/VisualizationStep+1,3,par.N); A1(1,:,:) = data.angular;
% for merged particles
PM = zeros(T/VisualizationStep+1,2,par.N); 
VM = zeros(T/VisualizationStep+1,2,par.N);

% ---------------------------- Iteration ---------------------------- %

j = 1; VisCounter = 0; ColCounter = 0;
c = DEM2Dcontacts(data,par);
for k = 1:T
    VisCounter = VisCounter +1;
    ColCounter = ColCounter +1;
    if ColCounter == par.CollisionStep
        ColCounter = 0;
        c = DEM2Dcontacts(data,par);
    end
    [pk,vk,ak,Pk,Vk,data] = DEM2Dsolve_expl(data,par,c);
%   [pk,vk,ak,data] = DEM2Dsolve_pgs(data,par,c.contacts);
    data.position = pk;
    data.velocity = vk;
    data.angular = ak;
    
    if VisCounter == par.VisualizationStep
        j = j+1; VisCounter = 0;
        P1(j,:,1:par.N) = data.position; 
        A1(j,:,1:par.N) = data.angular; 
        V1(j,:,1:par.N) = data.velocity;
        if(data.contactsParticle.mergedParticles)
            for kk = 1: data.contactsMerged.N
                if(data.contactsMerged.timeFlag(kk))
                    i = data.contactsMerged.index(1,kk); jj = data.contactsMerged.index(2,kk);
                    data.contactsMerged.timePoint(i,jj) = j;
                    data.contactsMerged.timeFlag(kk) = false;
                end
            end
            PM(j,:,1:par.N) = Pk;
            VM(j,:,1:par.N) = Vk;
        end
    end
end
% -------------------------- Plot time series -------------------------- %
DEM2DplotSim(P1,V1,A1,PM,VM,par,data,j)

time = 1:j;
for i = 1:par.N
plot(time,A1(:,2,i)*180/pi)
hold on
title(['Angular Velocity of particle ' num2str(i)])
end
hold off
figure;
for i = 1:par.N
plot(time,A1(:,1,i)*180/pi)
hold on
title(['Angle of particle ' num2str(i)])
end
%save('data')