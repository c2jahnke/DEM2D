% -------------------------- Initialization -------------------------- %
clc, clear all
global par;
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
% data.velocity(:,1) = [0.8;0];
% data.velocity(:,2) = [-0.8;0];
% data.position(:,1) = [0.1;-1.38];
% data.position(:,2) = [0.3;-1];

data.position(:,1) = [-0.5;-1];
% data.position(:,2) = [-0.3;-1];
% data.position(:,3) = [-0.1;-1];
% data.position(:,4) = [0.1;-1];
% data.position(:,5) = [0.4;-1];
data.velocity(:,1) = [0;0];
% data.velocity(:,2) = [0;0];
% data.velocity(:,3) = [0;0];
% data.velocity(:,4) = [0;0];
% data.velocity(:,5) = [-1.0;0];
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
algorithm = par.algorithm;
switch algorithm
    case 'PBD'
        for k = 1:T
            visuCounter = visuCounter +1;
            c = DEM2Dcontacts(data,par);
            [pk,vk,ak,acc,data] = DEM2Dsolve_pbd(par,data,c);
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
            end
        end
    case 'PGJ'
        for k = 1:T
            visuCounter = visuCounter +1;
            c = DEM2Dcontacts(data,par);
            [pk,vk,ak,acc,data] = DEM2Dsolve_pgj(par,data,c);
            data.position = pk;
            data.velocity = vk;
            data.angular = ak;
            data.acceleration = acc;
            if visuCounter == par.VisualizationStep
                if(mod(visuIndex,25) == 0)
                    disp(['################## ' sprintf('% 4d',visuIndex) '/' num2str(T/visuStep) ' frames <-> ' sprintf('% 4d',floor(round(visuIndex/T*visuStep*100))) ' Prozent ##################']);
                end
                visuIndex = visuIndex+1; visuCounter = 0;
%                 A(visuIndex,:,1:par.N) = data.acceleration;
                P1(visuIndex,:,1:par.N) = data.position; 
                A1(visuIndex,:,1:par.N) = data.angular; 
                V1(visuIndex,:,1:par.N) = data.velocity;
            end
        end
    case 'PGS'
        for k = 1:T
            visuCounter = visuCounter +1;
            c = DEM2Dcontacts(data,par);
            [pk,vk,ak,acc,data] = DEM2Dsolve_pgs(par,data,c);
            data.position = pk;
            data.velocity = vk;
            data.angular = ak;
            data.acceleration = acc;
            if visuCounter == par.VisualizationStep
                if(mod(visuIndex,25) == 0)
                    disp(['################## ' sprintf('% 4d',visuIndex) '/' num2str(T/visuStep) ' frames <-> ' sprintf('% 4d',floor(round(visuIndex/T*visuStep*100))) ' Prozent ##################']);
                end
                visuIndex = visuIndex+1; visuCounter = 0;
%                 A(visuIndex,:,1:par.N) = data.acceleration;
                P1(visuIndex,:,1:par.N) = data.position; 
                A1(visuIndex,:,1:par.N) = data.angular; 
                V1(visuIndex,:,1:par.N) = data.velocity;
            end
        end
    case 'HMD'
        for k = 1:T
            visuCounter = visuCounter +1;
            collisionCounter = collisionCounter +1;
            if collisionCounter == par.CollisionStep
                collisionCounter = 0;
                c = DEM2Dcontacts(data,par);
            end
            [pk,vk,ak,acc,data] = DEM2Dsolve_hmd(par,data,c);
            data.position = pk;
            data.velocity = vk;
            data.angular = ak;
            data.acceleration = acc;
            if visuCounter == par.VisualizationStep
                if(mod(visuIndex,25) == 0)
                    disp(['################## ' sprintf('% 4d',visuIndex) '/' num2str(T/visuStep) ' frames <-> ' sprintf('% 4d',floor(round(visuIndex/T*visuStep*100))) ' Prozent ##################']);
                end
                visuIndex = visuIndex+1; visuCounter = 0;
                A(visuIndex,:,1:par.N) = acc;
                P1(visuIndex,:,1:par.N) = pk; 
                A1(visuIndex,:,1:par.N) = ak; 
                V1(visuIndex,:,1:par.N) = vk;
            end
        end
    case  'LIN'
        for k = 1:T
            visuCounter = visuCounter +1;
            collisionCounter = collisionCounter +1;
           if collisionCounter == par.CollisionStep 
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
                [pk,vk,ak,acc,Pk,Vk,data] = DEM2Dsolve_linearDEM(par,data,c);
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
                            nonZeroIndex = nonzeros(data.contactsMerged.index(kk,:));
                            for jj = 1:data.contactsMerged.aggregateSize(kk) 
                                for ii = jj+1:data.contactsMerged.aggregateSize(kk) 
                                    data.contactsMerged.timePoint(nonZeroIndex(ii),nonZeroIndex(jj)) = visuIndex;
                                    data.contactsMerged.timePoint(nonZeroIndex(jj),nonZeroIndex(ii)) = visuIndex;
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
    case  'PDM' %Hybrid PBD-DEM
        for k = 1:T
            visuCounter = visuCounter +1;
            collisionCounter = collisionCounter +1;
           if collisionCounter == par.CollisionStep 
                % search particles in proximity of tool
                [partDist, unfrozenIndex] = DEM2DtoolDistance(par,data);
                index = ((partDist > 10*max(data.radius))');
                data.contactsParticle.deactivated = index';
                collisionCounter = 0;
                c = DEM2Dcontacts(data,par);
                
           end
            [pk,vk,ak,acc,Pk,Vk,data] = DEM2Dsolve_PBD_LIN(par,data,c);
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
            end
        end
    otherwise
        warning('Contact algorithm not available, check par.algorithm');
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