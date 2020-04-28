% -------------------------- Initialization -------------------------- %
global par;
par = DEM2Dparam();
global data;
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
%% Test 1 bottom
% data.position(1,1) = -1.0;
% data.position(2,1) = -1.08;
% data.velocity(1,1) = 2.5;
% data.velocity(2,1) = 0;
% % %% Test 2 top
% data.position(1,1) = -1.0;
% data.position(2,1) = 1.1;
% data.velocity(1,1) = 1.0;
% data.velocity(2,1) = 0;
% 
% data.angular(2,1) = 0;

% ------------------------ Plot initial state ------------------------ %
DEM2Dplot(data,par);
drawnow;

T = par.T; VisualizationStep = par.VisualizationStep; CollisionStep =  par.CollisionStep;
A = zeros(T/VisualizationStep+1,2,par.N);
P1 = zeros(T/VisualizationStep+1,2,par.N); P1(1,:,:) = data.position;
V1 = zeros(T/VisualizationStep+1,2,par.N); V1(1,:,:) = data.velocity;
A1 = zeros(T/VisualizationStep+1,2,par.N); A1(1,:,:) = data.angular;
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
    [pk,vk,ak,acceleration,Pk,Vk,data] = DEM2Dsolve_expl(data,par,c);
%   [pk,vk,ak,data] = DEM2Dsolve_pgs(data,par,c.contacts);
    data.position = pk;
    data.velocity = vk;
    data.angular = ak;
    data.acceleration = acceleration;
    
    if VisCounter == par.VisualizationStep
        if(mod(j,10) == 0)
            disp(['################## ' sprintf('% 4d',j) '/' num2str(T/VisualizationStep) ' frames <-> ' sprintf('% 4d',floor(round(j/T*VisualizationStep,2)*100)) ' Prozent ##################']);
            data.position
        end
        j = j+1; VisCounter = 0;
        A(j,:,1:par.N) = acceleration;
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

plot(time,A(:,:,1));

for i = 1:par.N
plot(time,A1(1:j,2,i)*180/pi)
hold on
title(['Angular Velocity of particle ' num2str(i)])
end

hold off
figure;

for i = 1:par.N
plot(time*par.dt/0.025,V1(1:j,2*i-1:2*i))
hold on
title(['Velocity of particle ' num2str(i)])
end

% for i = 1:par.N
% plot(time,A1(1:j,1,i)*180/pi)
% hold on
% title(['Angle of particle ' num2str(i)])
% end
% hold off
% figure;
% 
% for i = 1:par.N
% plot(time,0.5*data.mass(i)*(((V1(1:j,2*i,:)).^2)+((V1(1:j,2*i,:)).^2)))
% hold on
% title(['Ekin of particle ' num2str(i)])
% end