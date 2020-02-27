function [fwx,fwz,twy,data] = DEM2DwallForce(x,z,vx,vz,par,data)
    box = par.bBox; radius = data.radius;
    fwx_l = zeros(par.N,1);fwx_r = zeros(par.N,1);fwx_b = zeros(par.N,1);fwx_t = zeros(par.N,1);
    fwz_l = zeros(par.N,1);fwz_r = zeros(par.N,1);fwz_b = zeros(par.N,1);fwz_t = zeros(par.N,1);
    twy = zeros(par.N,4); 
    deltaW = zeros(4,par.N);
    ddeltaW = zeros(4,par.N);
    for i = 1 : par.N
        normal_stiffness = par.Emodul*pi/2*radius(i);
        tangential_stiffness = 1/1.2*normal_stiffness;
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
           normalConservative_l = normal_stiffness*deltaW(1,i);
           normalDissipative_l = 2*par.dampN*sqrt(data.mass(i)*normal_stiffness)*ddeltaW(1,i);
           fwx_l(i) =  normalConservative_l +normalDissipative_l;

           % tangential contact
           tangentialSpring = (data.contactsWall.actuationPoint(i,:,1) - [box(1) z(i)])*[0;1];
           if(abs(tangentialSpring) < 10*eps)
               fwz_l(i) = 0;
           else
               fwz_l(i) = (tangential_stiffness*tangentialSpring + 2*par.dampTwall*sqrt(data.mass(i)*tangential_stiffness)*vz(i));
               % friction
               if(abs(fwz_l(i)) > abs(fwx_l(i))*par.muWall)
                   fwz_l(i) = sign(fwz_l(i))*abs(fwx_l(i))*par.muWall;
                   data.contactsWall.actuationPoint(i,:,1) = [box(1) z(i)];%[x(i) box(1)];
               end
               
               if(par.considerRotations)
                   twy(i,1) = fwz_l(i)*data.radius(i);
               end
           end
           elseif(data.contactsWall.isInitialized(i,2)) % initialized but no contact with left wall
           data.contactsWall.contactAge(i,1) = data.contactsWall.contactAge(i,1) + 1;
           if(data.contactsWall.contactAge(i,1) >= data.contactsWall.maxContactAge)
               data.contactsWall.contactAge(i,1) = 0;
               data.contactsWall.isInitialized(i,1) = false;
               data.contactsWall.actuationPoint(i,:,1) = [0 0];
           end
        end
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
            normalConservative_r = -normal_stiffness*deltaW(2,i);
            normalDissipative_r = -2*par.dampN*sqrt(data.mass(i)*normal_stiffness)*ddeltaW(2,i);
            fwx_r(i) = (normalConservative_r + normalDissipative_r);
            
            % tangential contact
            tangentialSpring = (data.contactsWall.actuationPoint(i,:,2) - [box(2) z(i)])*[0;1];
            if(abs(tangentialSpring) < 10*eps)
                fwz_r(i) = 0;
            else
                fwz_r(i) = (tangential_stiffness*tangentialSpring + 2*par.dampTwall*sqrt(data.mass(i)*tangential_stiffness)*vz(i));

                % friction
                if(abs(fwz_r(i)) > abs(fwx_r(i))*par.muWall)
                    fwz_r(i) = sign(fwz_r(i))*abs(fwx_r(i))*par.muWall;
                    data.contactsWall.actuationPoint(i,:,2) = [box(2) z(i)];%[x(i) box(2)];
                end
                if(par.considerRotations)
                    twy(i,2) = fwz_r(i)*data.radius(i);
                end
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
                data.contactsWall.contactPoint(i,:,3) = DEM2Drotation(data.angular(2,i)*2*pi*par.dt)*data.contactsWall.contactPoint(i,:,3)'; % rotate contact point
            else
                data.contactsWall.isInitialized(i,3) = true;
                data.contactsWall.actuationPoint(i,:,3) = [x(i) box(3)];
                data.contactsWall.contactPoint(i,:,3) = [x(i) box(3)]; % rolling resistance
            end
            % normal contact
            deltaW(3,i) = radius(i) - abs(z(i)-box(3));
            ddeltaW(3,i) = -vz(i);
            fwz_b(i) = normal_stiffness*deltaW(3,i) + 2*par.dampN*sqrt(data.mass(i)*normal_stiffness)*ddeltaW(3,i);
            
            % tangential contact
            tangentialSpring = (data.contactsWall.actuationPoint(i,:,3) - [x(i) box(3)])*[1;0];
            if(abs(tangentialSpring) < 10*eps)
                fwx_b(i) = 0;
            else
                fwx_b(i) =tangential_stiffness*tangentialSpring +2*par.dampTwall*sqrt(data.mass(i)*tangential_stiffness)*vx(i);
                % friction
                if(abs(fwx_b(i)) > abs(fwz_b(i))*par.muWall)
                    % disp(["Bottom friction, index" num2str(i) "Ft_old =" num2str(fwx_b(i)) ", Fx_new =" num2str(sign(fwx_b(i))*abs(fwz_b(i))*par.muWall)])
                    fwx_b(i) = sign(fwx_b(i))*abs(fwz_b(i))*par.muWall;
                    % data.contactsWall.actuationPoint(i,:,3) = [x(i) box(3)];

                end
                if(par.considerRotations)
                twy(i,3) = fwx_b(i)*data.radius(i);
                end
                % rolling resistance
                if(par.considerRotations)
                    data.contactsWall.rollingDeformation(i,:,3) = data.contactsWall.contactPoint(i,:,3) - data.contactsWall.actuationPoint(i,:,3);
                    data.contactsWall.accumulatedRollingDeformation(i,:,3) = data.contactsWall.rollingDeformation(i,:,3);
                    % data.contactsWall.accumulatedRollingDeformation(i,:,3) + 
                    if(abs(data.contactsWall.accumulatedRollingDeformation(i,1,3)) > abs(fwz_b(i))*par.muWall/tangential_stiffness*0.1)
                        data.contactsWall.accumulatedRollingDeformation(i,:,3) = data.contactsWall.accumulatedRollingDeformation(i,:,3)*abs(fwz_b(i))*par.muWall/tangential_stiffness*0.1;
                    end
                    disp(["twy(i,3) no rolling resistance",twy(i,3)])
                    twy(i,3) = twy(i,3) + tangential_stiffness*data.contactsWall.accumulatedRollingDeformation(i,1,3)*data.radius(i); % projection into tangential plane necessary
                    disp(["twy(i,3) with rolling resistance",twy(i,3)])
                end
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
            fwz_t(i) = -normal_stiffness*deltaW(4,i) - 2*par.dampN*sqrt(data.mass(i)*normal_stiffness)*ddeltaW(4,i);%
             
            % tangential contact
            tangentialSpring = (data.contactsWall.actuationPoint(i,:,4) - [x(i) box(4)])*[1;0];
            if(abs(tangentialSpring) < 10*eps)
                fwx_t(i) = 0;
            else
                fwx_t(i) = tangential_stiffness*tangentialSpring - 2*par.dampTwall*sqrt(data.mass(i)*tangential_stiffness)*vx(i);
            end
            
            % friction            
            if(abs(fwx_t(i)) > abs(fwz_t(i))*par.muWall)
                fwx_t(i) = sign(fwx_t(i))*abs(fwz_t(i))*par.muWall;
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
        % Particles outside box should vanish
    end
fwx = fwx_l + fwx_r + fwx_b + fwx_t;
fwz = fwz_l + fwz_r + fwz_b + fwz_t;
end