function [pk,vk,ak,acc,data] = DEM2Dsolve_pgs(par,data,c)
% nonsmooth Projected Gau� Seidel method (from Jan Kleinerts code)
    contacts = c.contacts;
    pk = zeros(2,par.N);
    vk = zeros(2,par.N);
    ak = zeros(2,par.N);
    acc = zeros(2,par.N);
    w=par.w; % relaxation factor to control convergence. w<0 -> underrelaxation
           %                                           w>0 -> overrelaxation
    maxiter=20;
    iter = 0;
    numcontacts = length(contacts);
    
%     v0 = data.velocity;
    
    while iter < maxiter
        iter = iter +1;
        for j = 1:numcontacts
            
            bs = contacts(j).b;
            as = contacts(j).a;
            
            ui = contacts(j).rhs + contacts(j).Gb'*data.velocity(:,bs);
            if( as >0)
                 ui=ui+contacts(j).Ga'*data.velocity(:,as);        
            end
            delta = contacts(j).l - w*contacts(j).eta*ui;
            delta = delta + [par.cohesion; 0];
            delta = contacts(j).ConeProj(par,delta);
            delta = delta - contacts(j).l;
            delta = delta + [par.cohesion; 0];
            contacts(j).l = contacts(j).l + delta;
            
            data.velocity(:,bs) = data.velocity(:,bs)+1/contacts(j).redMass*contacts(j).Gb*delta;
            if(as>0)
                data.velocity(:,as)=data.velocity(:,as)+1/contacts(j).redMass*contacts(j).Ga*delta;
            end
        end
    end
    vk = data.velocity + [par.g_vert;par.g]*par.dt;
    pk = data.position + vk*par.dt;
end