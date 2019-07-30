function [fwx,fwz,data] = DEM2DwallForce(x,z,vx,vz,d,radius,par,data)
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
        
        % left
        if ((x(i)< box(1)+radius(i) && x(i) > box(1)))
           deltaW(2,i) = radius(i) - abs(x(i)-box(1));
           fwx(i) = normal_stiffness*deltaW(2,i);
           actuationPoint = [box(1) z(i)];
           if(data.velocity(1,i) < 0)
                data.velocity(1,i) =  0;%-par.wallDistr*data.velocity(1,i);
           end
        end

        % right
        if ((x(i)> box(2)-radius(i) && x(i) < box(2)))
           deltaW(4,i) = radius(i) - abs(x(i)-box(2));
           fwx(i) = - normal_stiffness*deltaW(4,i); 
           actuationPoint = [box(2) z(i)];
           if(data.velocity(1,i) > 0)
                data.velocity(1,i) = 0;% - par.wallDistr*data.velocity(1,i);
           end 
        end
       % bottom
        if ((z(i)<box(3)+radius(i) && z(i)>box(3)) && ( x(i)< box(3) || x(i)>box(1) ))
            if(data.contactsWall.isInitialized(i,3))
                data.contactsWall.contactAge(i,3) = 1;
            else
                data.contactsWall.isInitialized(i,3) = true;
                data.contactsWall.actuationPoint(i,:,3) = [x(i) box(3)];
            end
            % normal contact
            deltaW(4,i) = radius(i) - abs(z(i)-box(3));
            fwz(i) = normal_stiffness*deltaW(4,i)
            % tangential contact
            tangentialSpring = (data.contactsWall.actuationPoint(i,:,3) - [x(i) box(3)])*[1;0]
            if(abs(tangentialSpring) < eps)
                fwx(i) = 0;
            else
                %data.contactsWall.actuationPoint(i,:,3) = [x(i) box(3)];
                fwx(i) = par.kT*tangentialSpring
            end
            % friction
            
            if(abs(fwx(i)) > fwz(i)*par.mu)
                fwx(i) = sign(fwx(i))*fwz(i)*par.mu
            end
            data.velocity(2,i) = 0;
            data.acceleration(2,i) = 0;
        elseif(data.contactsWall.isInitialized(i,3)) % initialized but no contact with bottom
            data.contactsWall.contactAge(i,3) = data.contactsWall.contactAge(i,3) + 1;
            if(data.contactsWall.contactAge(i,3) >= data.contactsWall.maxContactAge)
                data.contactsWall.contactAge(i,3) = 0;
                data.contactsWall.isInitialized(i,3) = false;
                data.contactsWall.actuationPoint(i,:,3) = [0 0];
            end
        end

        % top
        if ((z(i)>box(4)-radius(i) && z(i)<box(4)) && ( x(i)< box(3) || x(i)>box(1) ))
            if(data.contactsWall.isInitialized(i,4))
                data.contactsWall.isInitialized(i,4) = true;
                data.contactsWall.contactAge(i,4) = 1;
            end
            deltaW(1,i) = radius(i) - abs(z(i)-box(4)); 
            fwz(i) = -normal_stiffness*deltaW(1,i); %
            actuationPoint = [x(i) box(4)];
            x12 = actuationPoint - [x(i) z(i)];
            dist = norm(x12)
            data.velocity(2,i) = 0;%- par.wallDistr*data.velocity(2,i);
            data.acceleration(2,i) = 0;
        elseif(data.contactsWall.isInitialized(i,4)) % no contact with top
            data.contactsWall.contactAge(i,4) = data.contactsWall.contactAge(i,4) + 1;
            if(data.contactsWall.contactAge(i,4) >= data.contactsWall.maxContactAge)
                data.contactsWall.contactAge(i,4) = 0;
                data.contactsWall.isInitialized(i,4) = false;
            end
        end
    end
end