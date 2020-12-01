% -------------------------- Initialization -------------------------- %
clc, clear all
global par;
%addpath('showcases/AngleOfRepose200p/');
par = DEM2Dparam();
global data;

LoadData = 0; % if true, load previous initial data
if(LoadData == true)
    SuccessFlag = false;
    [data,par] = DEM2Dload(par);
else
    [data,par,SuccessFlag] = DEM2Dinit(par);
if(SuccessFlag == 0)
    return
end
end
data.toolbBox = par.toolbBox;
% data.velocity(:,1) = [10;0];
% data.velocity(:,2) = [-2;0];
% data.position(:,1) = [-1;-1.1];
% data.position(:,2) = [1;-1];
% ------------------------ Plot initial state ------------------------ %
DEM2Dplot(data,par);
drawnow;

T = par.T; visuStep = par.VisualizationStep;
A = zeros(T/visuStep+1,2,par.N);
P1 = zeros(T/visuStep+1,2,par.N); P1(1,:,:) = data.position;
V1 = zeros(T/visuStep+1,2,par.N); V1(1,:,:) = data.velocity;
A1 = zeros(T/visuStep+1,2,par.N); A1(1,:,:) = data.angular;
% tool
PT = zeros(T/visuStep+1,2,2); PT(1,:,:) = data.toolbBox;
% for merged particles
PM = zeros(T/visuStep+1,2,par.N); 
VM = zeros(T/visuStep+1,2,par.N);
% for frozen particles
AP = zeros(T/visuStep+1,par.N); AP(1,:) = ones(1,par.N);
visuIndex = 1; visuCounter = 0; collisionCounter = 0;
c = DEM2Dcontacts(data,par);
% ---------------------------- Iteration ---------------------------- %
tic
for k = 1:T
    visuCounter = visuCounter +1;
    collisionCounter = collisionCounter +1;
   if collisionCounter == par.CollisionStep && ~par.PBD
        if par.Frozen
            [partDist, unfrozenIndex] = DEM2DtoolDistance(par,data);
            index = ((partDist > 10*max(data.radius))');
            data.contactsParticle.deactivated = index';
        end
        collisionCounter = 0;
        c = DEM2Dcontacts(data,par);
        if(par.merge && data.contactsMerged.N > 0)
            DEM2Dsplit(data,par)
        end
    end
    if par.PGJ
        [pk,vk,ak,acc,data] = DEM2Dsolve_pgj(par,data,c);
    elseif par.PBD
        [pk,vk,ak,acc,data] = DEM2Dsolve_pbd(par,data,c);
%         DEM2Dplot(data,par);
    elseif par.HMD
        [pk,vk,ak,acc,data] = DEM2Dsolve_hmd(par,data,c);
    else
        [pk,vk,ak,acc,Pk,Vk,data] = DEM2Dsolve_linearDEM(par,data,c);
    end
    data.position = pk;
    data.velocity = vk;
    data.angular = ak;
    data.acceleration = acc;

    if visuCounter == par.VisualizationStep
        if(mod(visuIndex,25) == 0)
            disp(['################## ' sprintf('% 4d',visuIndex) '/' num2str(T/visuStep) ' frames <-> ' sprintf('% 4d',floor(round(visuIndex/T*visuStep*100))) ' Prozent ##################']);
        end
        visuIndex = visuIndex+1; visuCounter = 0;
        A(visuIndex,:,1:par.N) = data.acceleration;
        P1(visuIndex,:,1:par.N) = data.position; 
        A1(visuIndex,:,1:par.N) = data.angular; 
        V1(visuIndex,:,1:par.N) = data.velocity;
        PT(visuIndex,:,:) = data.toolbBox;
        AP(visuIndex,:,:) = data.contactsParticle.deactivated';
        if(data.contactsParticle.mergedParticles)
            for kk = 1: data.contactsMerged.N
                if(data.contactsMerged.timeFlag(kk))
                    tmp = nonzeros(data.contactsMerged.index(kk,:));
                    for jj = 1:data.contactsMerged.aggregateSize(kk) 
                        for ii = jj+1:data.contactsMerged.aggregateSize(kk) 
                            data.contactsMerged.timePoint(tmp(ii),tmp(jj)) = visuIndex;
                            data.contactsMerged.timePoint(tmp(jj),tmp(ii)) = visuIndex;
                            data.contactsMerged.timeFlag(kk) = false;
                        end
                    end
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
output.deactivePart = AP;
output.MergePosition = PM;
output.MergeVelocity = VM;
output.timeInc = 1:visuIndex;
output.finalVisuIndex = visuIndex;

SimeTime = toc
% -------------------------- Plot time series -------------------------- %
DEM2DplotSim(P1,V1,A1,AP,PT,PM,VM,par,data,visuIndex)
% DEM3DplotDyn(P1,A1,data,par,visuIndex)


% DEM2DplotAngularVelocity(data,par,output)
DEM2DplotEnergy(data,par,output)
%pause(3)
close all