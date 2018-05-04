function [fwx,fwz] = DEM2DwallForce(x,z,vx,vz,d,r,par,data)
    box = par.bBox;
    fwx = zeros(par.N,1);
    fwz = zeros(par.N,1);
    p = data.position;
    n = [ 0 -1;
        -1  0;
         0  1;
         1  0];    
    deltaW = zeros(4,par.N);
    for i = 1 : par.N
        % top
        if ((z(i)>box(4)-r(i) && z(i)<box(4)) && ( x(i)< box(3) || x(i)>box(1) ))
            deltaW(1,i) = r(i) - abs(z(i)-box(4));
            
            fwz(i) = -31.4*r(i)*deltaW(1,i); % hard coded for concrete
            data.velocity(2,i) = - par.wallDistr*data.velocity(2,i);
        end
        
        % left
        if ((x(i)< box(1)+r(i) && x(i) > box(1)))
           deltaW(2,i) = r(i) - abs(x(i)-box(1));
           fwx(i) = 31.4*r(i)*deltaW(2,i);
           if(data.velocity(1,i) < 0)
                data.velocity(1,i) = -par.wallDistr*data.velocity(1,i);
           end
        end
         % bottom
        if ((z(i)<box(3)+r(i) && z(i)>box(3)) && ( x(i)< box(3) || x(i)>box(1) ))
            deltaW(4,i) = r(i) - abs(z(i)-box(3));
            
            fwz(i) = 31.4*r(i)*deltaW(4,i);
            data.velocity(2,i) = - par.wallDistr*data.velocity(2,i);
            data.acceleration(2,i) = 0;
        end
        % right
        if ((x(i)> box(2)-r(i) && x(i) < box(2)))
           deltaW(4,i) = r(i) - abs(x(i)-box(2));
           fwx(i) = - 31.4*r(i)*deltaW(4,i); % what is 31.4?
           if(data.velocity(1,i) > 0)
                data.velocity(1,i) = - par.wallDistr*data.velocity(1,i);
           end
        end
        %fwx = fwx;% 0.8
        %fwz = fwz;% 0.8
 
        
        
        
    end
    
    
    
end