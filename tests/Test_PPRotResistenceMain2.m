function TestValue = Test_PPRotResistenceMain2()
addpath('../')
TestValue = true;

% -------------------------- Initialization -------------------------- %
global par;
par = Test_PPRotResistenceParam2();
global data;
LoadData = false;
if(LoadData == true)
    SuccessFlag = false;
    data = DEM2Dload();
else
    [data,par,SuccessFlag] = DEM2Dinit(par);
if(SuccessFlag == 0)
    return
end
end
%% Rotation resistance p-p contact
data.position(1,1) = 0.879;
data.position(2,1) = -1.1;

data.position(1,2) = - 0.879;
data.position(2,2) = 0;
% data.position(1,3) = - 0.879;
% data.position(2,3) = -1;

data.velocity(1,1) = 0;
data.velocity(2,1) = 3.1;

data.velocity(1,2) = 0;
data.velocity(2,2) = 0;

% data.velocity(1,3) = 0;
% data.velocity(2,3) = 0;
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
    [pk,vk,ak,acceleration,Pk,Vk,data] = DEM2Dsolve_expl(par,data,c);
%   [pk,vk,ak,data] = DEM2Dsolve_pgs(data,par,c.contacts);
    data.position = pk;
    data.velocity = vk;
    data.angular = ak;
    data.acceleration = acceleration;
    

    if VisCounter == par.VisualizationStep
        if(mod(j,10) == 0)
            disp(['################## ' sprintf('% 4d',j) '/' num2str(T/VisualizationStep) ' frames <-> ' sprintf('% 4d',floor(round(j/T*VisualizationStep,2)*100)) ' Prozent ##################']);
            
            data.position
                data.position(1,1) = 0.879;
                data.velocity(1,1) = 0;
        end
        j = j+1; VisCounter = 0;
        A(j,:,1:par.N) = data.acceleration;
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
time = 1:j;
% -------------------------- Plot time series -------------------------- %
DEM2DplotSim(P1,V1,A1,PM,VM,par,data,j)
disp(["Slight rotation of particle 2."])

Test = data.angular
if(Test > 1e-1)
    TestValue = false;
end
rmpath('../')
end