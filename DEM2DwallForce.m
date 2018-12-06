function [fwx,fwz] = DEM2DwallForce(x,z,vx,vz,d,radius,par,data)
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
    normal_stiffness = par.Emodul*pi/2*radius(i);
        % top
        if ((z(i)>box(4)-radius(i) && z(i)<box(4)) && ( x(i)< box(3) || x(i)>box(1) ))
            deltaW(1,i) = radius(i) - abs(z(i)-box(4)); 
            fwz(i) = -normal_stiffness*radius(i)*deltaW(1,i); %
            actuationPoint = [x(i) box(4)];
            data.velocity(2,i) = 0;%- par.wallDistr*data.velocity(2,i);
            data.acceleration(2,i) = 0;
        end 
        % left
        if ((x(i)< box(1)+radius(i) && x(i) > box(1)))
           deltaW(2,i) = radius(i) - abs(x(i)-box(1));
           fwx(i) = normal_stiffness*radius(i)*deltaW(2,i);
           actuationPoint = [box(1) z(i)];
           if(data.velocity(1,i) < 0)
                data.velocity(1,i) =  0;%-par.wallDistr*data.velocity(1,i);
           end
        end
         % bottom
        if ((z(i)<box(3)+radius(i) && z(i)>box(3)) && ( x(i)< box(3) || x(i)>box(1) ))
            deltaW(4,i) = radius(i) - abs(z(i)-box(3));
            fwz(i) = normal_stiffness*radius(i)*deltaW(4,i);
            actuationPoint = [x(i) box(3)];
            data.velocity(2,i) = 0;% par.wallDistr*data.velocity(2,i);
            data.acceleration(2,i) = 0;
        end
        % right
        if ((x(i)> box(2)-radius(i) && x(i) < box(2)))
           deltaW(4,i) = radius(i) - abs(x(i)-box(2));
           fwx(i) = - normal_stiffness*radius(i)*deltaW(4,i); 
           actuationPoint = [box(2) z(i)];
           if(data.velocity(1,i) > 0)
                data.velocity(1,i) = 0;% - par.wallDistr*data.velocity(1,i);
           end 
        end
    end
end