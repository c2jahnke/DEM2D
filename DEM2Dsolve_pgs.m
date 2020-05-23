function [pk,vk,ak,acc,data] = DEM2Dsolve_pgs(par,data,c)
    contacts = c.contacts;
    pk = zeros(2,par.N);
    vk = zeros(2,par.N);
    ak = zeros(2,par.N);
    acc = zeros(2,par.N);
    w=0.2; % relaxation factor to control convergence. w<0 -> underrelaxation
           %                                           w>0 -> overrelaxation
    tol=-1;
    maxiter=20;
    iter = 0;
    numcontacts = length(contacts);
    
    v0 = data.velocity;
    
    while iter < maxiter
        iter = iter +1;
        for j = 1:numcontacts
            a = contacts(j).a;
            b = contacts(j).b;
             
            ui = contacts(j).rhs + contacts(j).Gb'*data.velocity(:,b);
            if( a >0)
                 ui=ui+contacts(j).Ga'*data.velocity(:,a);        
            end
            delta = contacts(j).l - w*contacts(j).eta*ui;
            %delta = delta + [par.cohesion; 0];
            contacts(j).l = contacts(j).ConeProj(par,delta);
            %contacts(j).l = contacts(j).l - [par.cohesion; 0];
        end
        data.velocity =  v0;
        for i=1:numcontacts
            a=contacts(i).a;
            b=contacts(i).b;
            data.velocity(:,b) = data.velocity(:,b)+1/contacts(i).redMass*contacts(i).Gb*contacts(i).l;
            if(a>0)
                data.velocity(:,a)=data.velocity(:,a)+1/contacts(i).redMass*contacts(i).Ga*contacts(i).l;
            end
        end
    end
    vk = data.velocity + [par.g_vert;par.g]*par.dt;
    pk = data.position + vk*par.dt;
end