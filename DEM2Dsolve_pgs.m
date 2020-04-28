function [pk,vk,ak,data] = DEM2Dsolve_pgs(data,par,contacts)
    w=0.2; % relaxation factor to control convergence. w<0 -> underrelaxation
           %                                           w>0 -> overrelaxation
    tol=-1;
    maxiter=1;
    iter = 0;
    numcontacts = length(contacts);
    
    while iter < maxiter
          iter = iter +1;
        for k = 1:numcontacts
            a = contacts(k).a;
            b = contacts(k).b;
           ui = contacts(k).rhs + contacts(k).Ga'*data.velocity(:,a);
           delta = contacts(k).l + w*contacts(k).eta*ui;
           delta = delta + [par.cohesion; 0];
           contacts(k).l = contacts(k).ConeProj(delta)
           contacts(k).l = 0;%contacts(k).l - [par.cohesion; 0];
           
        end
    end
end