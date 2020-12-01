function [pk,vk,ak,acc,data] = DEM2Dsolve_pbd(par,data,c)
% Position Based Dynamics (Miles, Macklin 2007, 2014, etc)
    nStab = 10; gamma = 0.3; nSteps = 10;
    dts = par.dt/nSteps;
    
    ak = zeros(2,par.N);
    acc = zeros(2,par.N);
    
    for k = 1:nSteps 
        xOld = data.position;
        v0 = data.velocity;
        x = data.position + dts*data.velocity + dts^2*([par.g_vert; par.g].*ones(2,par.N)-gamma*data.velocity);%
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
                    Das = 1/(contacts(j).distance/nSteps+data.radius(bs)+data.radius(as))*(x(:,as) - x(:,bs));
                    Dbs = 1/(contacts(j).distance/nSteps+data.radius(bs)+data.radius(as))*(x(:,bs) - x(:,as));
                    denom = Das'*1/data.mass(as)*Das + Dbs'*1/data.mass(bs)*Dbs;
                end
                deltaLambda = -overlap./denom;
                deltaXbs = 1/data.mass(bs)*Dbs*deltaLambda;
                x(:,bs) = x(:,bs) + deltaXbs;
                if(as>0)
                    deltaXas = -1/data.mass(as)*Das*deltaLambda;
                    x(:,as) = x(:,as) + deltaXas;
                    x(:,bs) = x(:,bs) - 2*deltaXbs;
                    % friction
%                     cf = [xOld(:,bs)-xOld(:,as) - (x(:,bs)-x(:,as))]'*contacts(j).t;
%                     frictionFactor = -(1/norm(x(:,bs)-x(:,as)))*(2 ...
%                         + (1/norm(x(:,bs)-x(:,as))^2)*(xOld(1,bs)-x(1,as) ...
%                         + xOld(2,bs)-xOld(2,as)- norm(x(:,bs)-x(:,as))^2));
% 
%                     Dfas = frictionFactor.*[x(:,as) - x(:,bs)];
%                     Dfbs = frictionFactor.*[x(:,bs) - x(:,as)];
%                     denom = Dfas' *1/data.mass(as)*Dfas + Dfbs'*1/data.mass(bs)*Dfbs;
%                     deltaLf = cf/denom;
%                     deltaLf = sign(deltaLf)*min(par.mu*deltaLambda,abs(deltaLf));
%                     deltaXfas = -1/data.mass(as)*Dfas*deltaLambda;
%                     deltaXfbs = -1/data.mass(bs)*Dfbs*deltaLambda;
%                     x(:,as) = x(:,as) + deltaXfas;
%                     x(:,bs) = x(:,bs) - deltaXfbs;
                end
            end
            data.velocity = (x - xOld)/dts;
        end
    end   
    data.position = x;
    pk = data.position;
    vk = data.velocity;
end