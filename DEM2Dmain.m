% -------------------------- Initialization -------------------------- %
global par;
par = DEM2Dparam();
global data;
LoadData = false;
if(LoadData == true)
    SuccessFlag = false;
    [data,par] = DEM2Dload(par);
else
    [data,par,SuccessFlag] = DEM2Dinit(par);
if(SuccessFlag == 0)
    return
end
end

% ------------------------ Plot initial state ------------------------ %
DEM2Dplot(data,par);
drawnow;

T = par.T; visuStep = par.VisualizationStep;
A = zeros(T/visuStep+1,2,par.N);
P1 = zeros(T/visuStep+1,2,par.N); P1(1,:,:) = data.position;
V1 = zeros(T/visuStep+1,2,par.N); V1(1,:,:) = data.velocity;
A1 = zeros(T/visuStep+1,2,par.N); A1(1,:,:) = data.angular;
% for merged particles
PM = zeros(T/visuStep+1,2,par.N); 
VM = zeros(T/visuStep+1,2,par.N);

% ---------------------------- Iteration ---------------------------- %
visuIndex = 1; visuCounter = 0; collisionCounter = 0;
c = DEM2Dcontacts(data,par);
for k = 1:T
    visuCounter = visuCounter +1;
    collisionCounter = collisionCounter +1;
    if collisionCounter == par.CollisionStep
        collisionCounter = 0;
        c = DEM2Dcontacts(data,par);
    end
    [pk,vk,ak,acc,Pk,Vk,data] = DEM2Dsolve_expl(par,data,c);
%   [pk,vk,ak,data] = DEM2Dsolve_pgs(data,par,c.contacts);
    data.position = pk;
    data.velocity = vk;
    data.angular = ak;
    data.acceleration = acc;

    if visuCounter == par.VisualizationStep
        if(mod(visuIndex,25) == 0)
            disp(['################## ' sprintf('% 4d',visuIndex) '/' num2str(T/visuStep) ' frames <-> ' sprintf('% 4d',floor(round(visuIndex/T*visuStep,2)*100)) ' Prozent ##################']);
        end
        visuIndex = visuIndex+1; visuCounter = 0;
        A(visuIndex,:,1:par.N) = data.acceleration;
        P1(visuIndex,:,1:par.N) = data.position; 
        A1(visuIndex,:,1:par.N) = data.angular; 
        V1(visuIndex,:,1:par.N) = data.velocity;
        if(data.contactsParticle.mergedParticles)
            for kk = 1: data.contactsMerged.N
                if(data.contactsMerged.timeFlag(kk))
                    i = data.contactsMerged.index(1,kk); jj = data.contactsMerged.index(2,kk);
                    data.contactsMerged.timePoint(i,jj) = visuIndex;
                    data.contactsMerged.timeFlag(kk) = false;
                end
            end
            PM(visuIndex,:,1:par.N) = Pk;
            VM(visuIndex,:,1:par.N) = Vk;
        end
    end
end
output = struct;
output.acceleration = A;
output.position = P1;
output.angular = A1;
output.velocity = V1;
output.MergePosition = PM;
output.MergeVelocity = VM;
output.timeInc = 1:visuIndex;
output.finalVisuIndex = visuIndex;
% -------------------------- Plot time series -------------------------- %
DEM2DplotSim(P1,V1,A1,PM,VM,par,data,visuIndex)
% DEM3DplotDyn(P1,A1,data,par,visuIndex)


% DEM2DplotAngularVelocity(data,par,output)
DEM2DplotEnergy(data,par,output)