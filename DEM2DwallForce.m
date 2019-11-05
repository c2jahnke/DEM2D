function [fwx,fwz,data] = DEM2DwallForce(x,z,vx,vz,d,radius,par,data)
    box = par.bBox;
    fwx = zeros(par.N,1); fwz = zeros(par.N,1);
    fwx_l = zeros(par.N,1);fwx_r = zeros(par.N,1);fwx_b = zeros(par.N,1);fwx_t = zeros(par.N,1);
    fwz_l = zeros(par.N,1);fwz_r = zeros(par.N,1);fwz_b = zeros(par.N,1);fwz_t = zeros(par.N,1);
    p = data.position;
    n = [ 0 -1;
        -1  0;
         0  1;
         1  0];    
    deltaW = zeros(4,par.N);
    ddeltaW = zeros(4,par.N);
    for i = 1 : par.N
    normal_stiffness = par.Emodul*pi/2*radius(i);
        
        % left
        if ((x(i)< box(1)+radius(i) && x(i) > box(1)))
           if(data.contactsWall.isInitialized(i,1))
                data.contactsWall.contactAge(i,1) = 1;
            else
                data.contactsWall.isInitialized(i,1) = true;
                data.contactsWall.actuationPoint(i,:,1) = [box(1) z(i)];
            end
            % normal contact
            deltaW(1,i) = radius(i) - abs(x(i)-box(1));
            ddeltaW(1,i) = -vx(i);
            fwx_l(i) = normal_stiffness*deltaW(1,i) + par.dampN*sqrt(data.mass(i)*normal_stiffness)*ddeltaW(1,i);
            
            % tangential contact
            tangentialSpring = (data.contactsWall.actuationPoint(i,:,1) - [box(1) z(i)])*[1;0];
            if(abs(tangentialSpring) < 10*eps)
                fwz_l(i) = 0;
            else
                fwz_l(i) = par.kT*tangentialSpring;
            end
            % friction
            if(abs(fwz_l(i)) > abs(fwx_l(i))*par.mu)
                fwz_l(i) = sign(fwz_l(i))*fwx_l(i)*par.mu;
                data.contactsWall.actuationPoint(i,:,1) = [x(i) box(1)];
            end
        elseif(data.contactsWall.isInitialized(i,2)) % initialized but no contact with bottom
            data.contactsWall.contactAge(i,1) = data.contactsWall.contactAge(i,1) + 1;
            if(data.contactsWall.contactAge(i,1) >= data.contactsWall.maxContactAge)
                data.contactsWall.contactAge(i,1) = 0;
                data.contactsWall.isInitialized(i,1) = false;
                data.contactsWall.actuationPoint(i,:,1) = [0 0];
            end
        end
%            deltaW(1,i) = radius(i) - abs(x(i)-box(1));
%            ddeltaW(1,i) = -vx(i);
%            fwx_l(i) = normal_stiffness*deltaW(1,i) + par.dampN*sqrt(data.mass(i)*normal_stiffness)*ddeltaW(1,i);
% %            disp(['normal force:' num2str(normal_stiffness*deltaW(1,i))]);
% %            disp(['dissipative force:' num2str(par.dampN*sqrt(data.mass(i)*normal_stiffness)*ddeltaW(1,i))]);
%         
%            actuationPoint = [box(1) z(i)];
%            if(data.velocity(1,i) < 0)
%                 data.velocity(1,i) =  0;%-par.wallDistr*data.velocity(1,i);
%            end
%         end

        % right
        if ((x(i)> box(2)-radius(i) && x(i) < box(2)))
            if(data.contactsWall.isInitialized(i,2))
                data.contactsWall.contactAge(i,2) = 1;
            else
                data.contactsWall.isInitialized(i,2) = true;
                data.contactsWall.actuationPoint(i,:,2) = [box(2) z(i)];
            end
            % normal contact
            deltaW(2,i) = radius(i) - abs(x(i)-box(2));
            ddeltaW(2,i) = vx(i);
            fwx_r(i) = - normal_stiffness*deltaW(2,i) - par.dampN*sqrt(data.mass(i)*normal_stiffness)*ddeltaW(2,i);
            
            % tangential contact
            tangentialSpring = (data.contactsWall.actuationPoint(i,:,2) - [box(2) z(i)])*[1;0];
            if(abs(tangentialSpring) < 10*eps)
                fwz_r(i) = 0;
            else
                fwz_r(i) = par.kT*tangentialSpring;
            end
            % friction
            if(abs(fwz_r(i)) > abs(fwx_r(i))*par.mu)
                fwz_r(i) = sign(fwz_r(i))*fwx_r(i)*par.mu;
                data.contactsWall.actuationPoint(i,:,2) = [x(i) box(2)];
            end
        elseif(data.contactsWall.isInitialized(i,2)) % initialized but no contact with bottom
            data.contactsWall.contactAge(i,2) = data.contactsWall.contactAge(i,2) + 1;
            if(data.contactsWall.contactAge(i,2) >= data.contactsWall.maxContactAge)
                data.contactsWall.contactAge(i,2) = 0;
                data.contactsWall.isInitialized(i,2) = false;
                data.contactsWall.actuationPoint(i,:,2) = [0 0];
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
            deltaW(3,i) = radius(i) - abs(z(i)-box(3));
            ddeltaW(3,i) = -vz(i);
            fwz_b(i) = normal_stiffness*deltaW(3,i) + par.dampN*sqrt(data.mass(i)*normal_stiffness)*ddeltaW(3,i);
            
            % tangential contact
            tangentialSpring = (data.contactsWall.actuationPoint(i,:,3) - [x(i) box(3)])*[1;0];
            if(abs(tangentialSpring) < eps)
                fwx_b(i) = 0;
            else
                fwx_b(i) = par.kT*tangentialSpring;
            end
            % friction
            if(abs(fwx_b(i)) > abs(fwz(i))*par.mu)
                fwx_b(i) = sign(fwx_b(i))*fwz_b(i)*par.mu;
                data.contactsWall.actuationPoint(i,:,3) = [x(i) box(3)];
            end
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
            else
                data.contactsWall.isInitialized(i,4) = true;
                data.contactsWall.actuationPoint(i,:,4) = [x(i) box(4)];
            end
            % normal contact
            deltaW(4,i) = radius(i) - abs(z(i)-box(4)); 
            ddeltaW(4,i) = vz(i);
            fwz_t(i) = -normal_stiffness*deltaW(4,i) - par.dampN*sqrt(data.mass(i)*normal_stiffness)*ddeltaW(4,i);%
             
            % tangential contact
            tangentialSpring = (data.contactsWall.actuationPoint(i,:,4) - [x(i) box(4)])*[1;0];
            if(abs(tangentialSpring) < eps)
                fwx_t(i) = 0;
            else
                fwx_t(i) = -par.kT*tangentialSpring;
            end
            % friction
            
            if(abs(fwx_t(i)) > abs(fwz_t(i))*par.mu)
                fwx_t(i) = sign(fwx_t(i))*fwz_t(i)*par.mu;
                data.contactsWall.actuationPoint(i,:,4) = [x(i) box(4)];
            end

        elseif(data.contactsWall.isInitialized(i,4)) % initialized but no contact with top
            data.contactsWall.contactAge(i,4) = data.contactsWall.contactAge(i,4) + 1;
            if(data.contactsWall.contactAge(i,4) >= data.contactsWall.maxContactAge)
                data.contactsWall.contactAge(i,4) = 0;
                data.contactsWall.isInitialized(i,4) = false;
                data.contactsWall.actuationPoint(i,:,3) = [0 0];
            end
        end
    end
fwx = fwx_l + fwx_r + fwx_b + fwx_t;
fwz = fwz_l + fwz_r + fwz_b + fwz_t;
end
