function TestValue = TEST_RotMain3()
addpath('../')
TestValue = true;
% -------------------------- Initialization -------------------------- %
global par;
par = TEST_RotParam3();
global data;
LoadData = false;
if(LoadData == true)
    SuccessFlag = true;
    data = DEM2Dload();
else
    [data,par,SuccessFlag] = DEM2Dinit(par);
if(SuccessFlag == 0)
    return
end
end
%% Test 1 bottom
data.position(1,1) = -1.0;
data.position(2,1) = -1.1;
data.velocity(1,1) = 1;
data.velocity(2,1) = 0;


data.angular(2,1) = 0;

% ------------------------ Plot initial state ------------------------ %
% DEM2Dplot(data,par);
% drawnow;

T = par.T; VisualizationStep = par.VisualizationStep; CollisionStep =  par.CollisionStep;
A = zeros(T/VisualizationStep+1,2,par.N);
P1 = zeros(T/VisualizationStep+1,2,par.N); P1(1,:,:) = data.position;
V1 = zeros(T/VisualizationStep+1,2,par.N); V1(1,:,:) = data.velocity;
A1 = zeros(T/VisualizationStep+1,2,par.N); A1(1,:,:) = data.angular;
% tool
PT = zeros(T/VisualizationStep+1,2,2);
% for merged particles
PM = zeros(T/VisualizationStep+1,2,par.N); 
VM = zeros(T/VisualizationStep+1,2,par.N);
% for frozen particles
AP = zeros(T/VisualizationStep+1,par.N); AP(1,:) = ones(1,par.N);
visuIndex = 1; visuCounter = 0; collisionCounter = 0;
% ---------------------------- Iteration ---------------------------- %

% ---------------------------- Iteration ---------------------------- %

visuIndex = 1; VisCounter = 0; ColCounter = 0;
c = DEM2Dcontacts(data,par);
for k = 1:T
    VisCounter = VisCounter +1;
    ColCounter = ColCounter +1;
    if ColCounter == par.CollisionStep
        ColCounter = 0;
        c = DEM2Dcontacts(data,par);
    end
    [pk,vk,ak,acc,Pk,Vk,data] = DEM2Dsolve_expl(par,data,c);
%     [pk,vk,ak,acc,data] = DEM2Dsolve_pgs(par,data,c);
    data.position = pk;
    data.velocity = vk;
    data.angular = ak;
    data.acceleration = acc;
    
    if VisCounter == par.VisualizationStep
        if(mod(visuIndex,10) == 0)
            disp(['################## ' sprintf('% 4d',visuIndex) '/' num2str(T/VisualizationStep) ' frames <-> ' sprintf('% 4d',floor(round(visuIndex/T*VisualizationStep,2)*100)) ' Prozent ##################']);
            data.position
        end
        visuIndex = visuIndex+1; VisCounter = 0;
        A(visuIndex,:,1:par.N) = acc;
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
% -------------------------- Plot time series -------------------------- %
% DEM2DplotSim(P1,V1,A1,AP,PT,PM,VM,par,data,visuIndex)

Test = data.angular
Test2 = data.position
if(Test2 - [0.7344; -1.1000] > 1e-1)
    TestValue = false;
end
rmpath('../')
end