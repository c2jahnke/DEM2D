function [pk,vk,ak,acc,data] = DEM2Dsolve_pbd(par,data,c)
% Position Based Dynamics (Miles, Macklin 2007, 2014, etc)
    nStab = 10; gamma = 0.3; nSteps = 5;
    dts = par.dt/nSteps;
    
%     pk = zeros(2,par.N);
%     vk = zeros(2,par.N);
    ak = zeros(2,par.N);
    acc = zeros(2,par.N);
    
    for k = 1:nSteps 
        xOld = data.position;
        v0 = data.velocity;
        data.position = data.position + dts*data.velocity + dts^2*([par.g_vert; par.g].*ones(2,par.N)-gamma*data.velocity);%
        c = DEM2Dcontacts(data,par);
        contacts = c.contacts;
        numcontacts = length(contacts);
        % postion_based_dynamics.jl 92 - 113
        for j = 1:numcontacts
            bs = contacts(j).b;
            as = contacts(j).a;
            overlap = contacts(j).distance/nSteps;
            if overlap > 0
                Dbs = 1/(contacts(j).distance/nSteps+data.radius(bs))*contacts(j).n;
                denom = Dbs'*1/data.mass(bs)*Dbs;
                if(as>0)
                    Das = 1/(contacts(j).distance/nSteps+data.radius(bs)+data.radius(as))*(data.position(:,as) - data.position(:,bs));
                    Dbs = 1/(contacts(j).distance/nSteps+data.radius(bs)+data.radius(as))*(data.position(:,bs) - data.position(:,as));
                    denom = Das'*1/data.mass(as)*Das + Dbs'*1/data.mass(bs)*Dbs;
                end
                deltaLambda = -overlap./denom;
                deltaXbs = 1/data.mass(bs)*Dbs*deltaLambda;
                data.position(:,bs) = data.position(:,bs) + deltaXbs;
                if(as>0)
                    deltaXas = -1/data.mass(as)*Das*deltaLambda;
                    data.position(:,as) = data.position(:,as) + deltaXas;
                    data.position(:,bs) = data.position(:,bs) - 2*deltaXbs;
                end
            end
            data.velocity = (data.position - xOld)/dts;
        end
    end   
    pk = data.position;
    vk = data.velocity;
end