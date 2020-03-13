function [fwx,fwz,twy,data] = DEM2DwallForce(x,z,vx,vz,par,data)
    box = par.bBox; radius = data.radius;
    fwx_l = zeros(par.N,1);fwx_r = zeros(par.N,1);fwx_b = zeros(par.N,1);fwx_t = zeros(par.N,1);
    fwz_l = zeros(par.N,1);fwz_r = zeros(par.N,1);fwz_b = zeros(par.N,1);fwz_t = zeros(par.N,1);
    twy = zeros(par.N,4); 
    deltaW = zeros(4,par.N);
    for i = 1 : par.N
        normalStiffness = par.Emodul*pi/2*radius(i);
        tangentialStiffness = 1/1.2*normalStiffness;
        % left
        if ((x(i)< box(1)+radius(i)))% && x(i) > box(1)))
           if(data.contactsWall.isInitialized(i,1))
                data.contactsWall.contactAge(i,1) = 1;               
                data.contactsWall.contactPoint(i,:,1) = (data.contactsWall.localContactPoint(i,:,1))'+ [x(i) z(i)]';%*data.contactsWall.contactPoint(i,:,1)'; % rotate contact point
                
           else
                data.contactsWall.isInitialized(i,1) = true;
                data.contactsWall.actuationPoint(i,:,1) = [box(1) z(i)];
                data.contactsWall.contactPoint(i,:,1) = [box(1) z(i)]; % rolling resistance
                data.contactsWall.localContactPoint(i,:,1) = data.contactsWall.contactPoint(i,:,1) - [x(i) z(i)]; % rolling resistance
           end
           % normal contact
           deltaW(1,i) = radius(i) - abs(x(i)-box(1));
           %ddeltaW(1,i) = -vx(i);
           normalConservative_l = normalStiffness*deltaW(1,i);
           normalDissipative_l = par.dampN*2*sqrt(0.5*data.mass(i)*normalStiffness)*vx(i);
           fwx_l(i) =  normalConservative_l - normalDissipative_l;

           % tangential contact
           tangentialSpring = (data.contactsWall.actuationPoint(i,:,1) - data.contactsWall.contactPoint(i,:,1))*[0;1];
           if(abs(tangentialSpring) < 100*eps)
               fwz_l(i) = 0;
           else
               fwz_l(i) = (tangentialStiffness*tangentialSpring - par.dampTwall*2*sqrt(0.5*data.mass(i)*tangentialStiffness)*vz(i));
               % friction
               if(abs(fwz_l(i)) > abs(fwx_l(i))*par.muWall)
                   fwz_l(i) = sign(fwz_l(i))*abs(fwx_l(i))*par.muWall;
                   %data.contactsWall.actuationPoint(i,:,1) = [box(1) z(i)];%[x(i) box(1)];
                   data.contactsWall.contactPoint(i,:,1) = data.contactsWall.actuationPoint(i,:,1) + tangentialSpring/2*[0;1]';
                   data.contactsWall.localContactPoint(i,:,1) = (data.contactsWall.contactPoint(i,:,1))'-[x(i) z(i)]';
               end
               
               %DEM2Drotation(data.angular(2,i)*par.dt)*
               if(par.considerRotations)
           %        twy(i,1) = fwz_l(i)*norm(data.contactsWall.actuationPoint(i,:,1)'-data.position(:,i));
                   twy(i,1) = fwz_l(i)*(data.contactsWall.actuationPoint(i,1,1)-data.position(1,i));
               end
               if(par.considerRotations)
                    data.contactsWall.rollingDeformation(i,1,1) = data.contactsWall.contactPoint(i,2,1) - data.contactsWall.actuationPoint(i,2,1); % 4.27
                    data.contactsWall.accumulatedRollingDeformation(i,1,1) = data.contactsWall.accumulatedRollingDeformation(i,1,1) + data.contactsWall.rollingDeformation(i,1,1);%4.28
                    
                    if(abs(data.contactsWall.accumulatedRollingDeformation(i,1,1)) > abs(fwx_l(i))*par.muWall/tangentialStiffness*par.Cr) %4.30
                        data.contactsWall.accumulatedRollingDeformation(i,1,1) = data.contactsWall.accumulatedRollingDeformation(i,1,1)/abs(data.contactsWall.accumulatedRollingDeformation(i,1,1))*abs(fwx_l(i))*(par.muWall/tangentialStiffness)*par.Cr;
                    end
%                       disp(["twy(i,1) no rolling resistance",twy(i,1)])
                      twy(i,1) = twy(i,1) + tangentialStiffness*data.contactsWall.accumulatedRollingDeformation(i,1,1)*(data.contactsWall.actuationPoint(i,1,1)'-data.position(1,i)); % projection into tangential plane necessary
%                       disp(["twy(i,1) with rolling resistance",twy(i,1)])
                end

           end
           elseif(data.contactsWall.isInitialized(i,1)) % initialized but no contact with left wall
           data.contactsWall.contactAge(i,1) = data.contactsWall.contactAge(i,1) + 1;
           if(data.contactsWall.contactAge(i,1) >= data.contactsWall.maxContactAge)
               data.contactsWall.contactAge(i,1) = 0;
               data.contactsWall.isInitialized(i,1) = false;
               data.contactsWall.actuationPoint(i,:,1) = [0 0];
           end
        end
        % right
        if ((x(i)> box(2)-radius(i)))%&& x(i) < box(2)))
            if(data.contactsWall.isInitialized(i,2))
                data.contactsWall.contactAge(i,2) = 1; %DEM2Drotation(data.angular(2,i)*par.dt)*
                data.contactsWall.contactPoint(i,:,2) = (data.contactsWall.localContactPoint(i,:,2))'+ [x(i) z(i)]'; % rotate contact point
            else
                data.contactsWall.isInitialized(i,2) = true;
                data.contactsWall.actuationPoint(i,:,2) = [box(2) z(i)];
                data.contactsWall.contactPoint(i,:,2) = [box(2) z(i)]; % rolling resistance
                data.contactsWall.localContactPoint(i,:,2) = data.contactsWall.contactPoint(i,:,2) - [x(i) z(i)]; % rolling resistance
            end
        
            % normal contact
            deltaW(2,i) = radius(i) - abs(x(i)-box(2));
            normalConservative_r = -normalStiffness*deltaW(2,i);
            normalDissipative_r = par.dampN*2*sqrt(0.5*data.mass(i)*normalStiffness)*vx(i);
            fwx_r(i) = (normalConservative_r - normalDissipative_r);
            
            % tangential contact
            tangentialSpring = (data.contactsWall.actuationPoint(i,:,2) - data.contactsWall.contactPoint(i,:,2))*[0;1];
            if(abs(tangentialSpring) < 10*eps)
                fwz_r(i) = 0;
            else
                fwz_r(i) = (tangentialStiffness*tangentialSpring - par.dampTwall*2*sqrt(0.5*data.mass(i)*tangentialStiffness)*vz(i));

                % friction
                if(abs(fwz_r(i)) > abs(fwx_r(i))*par.muWall)
                    fwz_r(i) = sign(fwz_r(i))*abs(fwx_r(i))*par.muWall;
                    data.contactsWall.contactPoint(i,:,2) = data.contactsWall.actuationPoint(i,:,2) + tangentialSpring/2*[0;1]';
                    data.contactsWall.localContactPoint(i,:,2) = (data.contactsWall.contactPoint(i,:,2))'-[x(i) z(i)]';%(data.contactsWall.contactPoint(i,:,2)')-data.position(:,i);
     
                end
             if(par.considerRotations)
                 twy(i,2) = fwz_r(i)*(data.contactsWall.actuationPoint(i,1,2)-data.position(1,i));
             end
            % rolling resistance
            if(par.considerRotations)
                data.contactsWall.rollingDeformation(i,1,2) = data.contactsWall.contactPoint(i,2,2) - data.contactsWall.actuationPoint(i,2,2); % 4.27
                data.contactsWall.accumulatedRollingDeformation(i,1,2) = data.contactsWall.accumulatedRollingDeformation(i,1,2) + data.contactsWall.rollingDeformation(i,1,2);%4.28

                if(abs(data.contactsWall.accumulatedRollingDeformation(i,1,2)) > abs(fwx_r(i))*par.muWall/tangentialStiffness*par.Cr) %4.30
                    data.contactsWall.accumulatedRollingDeformation(i,1,2) = data.contactsWall.accumulatedRollingDeformation(i,1,2)/abs(data.contactsWall.accumulatedRollingDeformation(i,1,2))*abs(fwx_r(i))*par.muWall/tangentialStiffness*par.Cr; %data.contactsWall.accumulatedRollingDeformation(i,1,2)/abs(data.contactsWall.accumulatedRollingDeformation(i,1,2))*
                    %data.contactsWall.contactPoint(i,:,2) = data.contactsWall.actuationPoint(i,:,2);
                end
%                  disp(["twy(i,2) no rolling resistance",twy(i,2)])
                 twy(i,2) = twy(i,2) + tangentialStiffness*data.contactsWall.accumulatedRollingDeformation(i,1,2)*(data.contactsWall.actuationPoint(i,1,2)-data.position(1,i)); % projection into tangential plane necessary
%                  disp(["twy(i,2) with rolling resistance",twy(i,2)])
            end
            end           
        elseif(data.contactsWall.isInitialized(i,2)) % initialized but no contact with bottom
            data.contactsWall.contactAge(i,2) = data.contactsWall.contactAge(i,2) + 1;
            if(data.contactsWall.contactAge(i,2) >= data.contactsWall.maxContactAge)
                data.contactsWall.contactAge(i,2) = 0;
                data.contactsWall.isInitialized(i,2) = false;
                data.contactsWall.actuationPoint(i,:,2) = [0 0];
                data.contactsWall.localContactPoint(i,:,2) = [0 0];
                data.contactsWall.contactPoint(i,:,2) = [0 0];
            end
        end

        % bottom
        if ((z(i)<box(3)+radius(i)))% && z(i)>box(3)) && ( x(i)< box(3) || x(i)>box(1) ))
            if(data.contactsWall.isInitialized(i,3))
                data.contactsWall.contactAge(i,3) = 1;%DEM2Drotation(data.angular(2,i)*par.dt)*
                data.contactsWall.contactPoint(i,:,3) = data.contactsWall.localContactPoint(i,:,3)+ [x(i) z(i)]; %*data.contactsWall.contactPoint(i,:,3)'; % rotate contact point
            else
                data.contactsWall.isInitialized(i,3) = true;
                data.contactsWall.actuationPoint(i,:,3) = [x(i) box(3)];
                data.contactsWall.contactPoint(i,:,3) = [x(i) box(3)]; % rolling resistance
                data.contactsWall.localContactPoint(i,:,3) = data.contactsWall.contactPoint(i,:,3) - [x(i) z(i)]; % rolling resistance
            end
            
            % normal contact
            deltaW(3,i) = radius(i) - abs(z(i)-box(3));
            normalConservative_b = normalStiffness*deltaW(3,i);
            normalDissipative_b = par.dampN*2*sqrt(0.5*data.mass(i)*normalStiffness)*vz(i);
            fwz_b(i) = normalConservative_b  - normalDissipative_b ;

            % tangential contact
            tangentialSpring = (data.contactsWall.actuationPoint(i,:,3) - data.contactsWall.contactPoint(i,:,3))*[1;0];
            if(abs(tangentialSpring) < 100*eps)
                fwx_b(i) = 0;
            else
                fwx_b(i) = tangentialStiffness*tangentialSpring +par.dampTwall*2*sqrt(0.5*data.mass(i)*tangentialStiffness)*vx(i);
                % friction
                if(abs(fwx_b(i)) > abs(fwz_b(i))*par.muWall)
                    fwx_b(i) = sign(fwx_b(i))*abs(fwz_b(i))*par.muWall;
                    data.contactsWall.contactPoint(i,:,3) = data.contactsWall.actuationPoint(i,:,3) + tangentialSpring/2*[0;1]';
                    data.contactsWall.localContactPoint(i,:,3) = (data.contactsWall.contactPoint(i,:,3))'-[x(i) z(i)]';
                end
            end
%             disp(["fwx_b(i)",fwx_b(i)])
            if(par.considerRotations)
               twy(i,3) = fwx_b(i)*(data.contactsWall.actuationPoint(i,2,3)-data.position(2,i));
            end
            % rolling resistance
            if(par.considerRotations)
                data.contactsWall.rollingDeformation(i,1,3) = data.contactsWall.contactPoint(i,1,3) - data.contactsWall.actuationPoint(i,1,3); % 4.27
                data.contactsWall.accumulatedRollingDeformation(i,1,3) = data.contactsWall.accumulatedRollingDeformation(i,1,3) + data.contactsWall.rollingDeformation(i,1,3);%4.28

                if(abs(data.contactsWall.accumulatedRollingDeformation(i,1,3)) > abs(fwz_b(i))*par.muWall/tangentialStiffness*par.Cr) %4.30
                    data.contactsWall.accumulatedRollingDeformation(i,1,3) =data.contactsWall.accumulatedRollingDeformation(i,1,3)/abs(data.contactsWall.accumulatedRollingDeformation(i,1,3))*(fwz_b(i))*par.muWall/tangentialStiffness*par.Cr; %data.contactsWall.accumulatedRollingDeformation(i,1,3)/abs(data.contactsWall.accumulatedRollingDeformation(i,1,3))*

                end
                disp(["twy(i,3) no rolling resistance",twy(i,3)])
                twy(i,3) = twy(i,3) + tangentialStiffness*data.contactsWall.accumulatedRollingDeformation(i,1,3)*(data.contactsWall.actuationPoint(i,2,3)'-data.position(2,i)); % projection into tangential plane necessary           
                disp(["twy(i,3) with rolling resistance",twy(i,3)])
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
        if ((z(i)>box(4)-radius(i) && z(i)<box(4)))%% && ( x(i)< box(3) || x(i)>box(1) ))
            if(data.contactsWall.isInitialized(i,4))
                data.contactsWall.contactAge(i,4) = 1;
                data.contactsWall.contactPoint(i,:,4) = (data.contactsWall.localContactPoint(i,:,4))'+ [x(i) z(i)]';
            else
                data.contactsWall.isInitialized(i,4) = true;
                data.contactsWall.actuationPoint(i,:,4) = [x(i) box(4)];
                data.contactsWall.contactPoint(i,:,4) = [x(i) box(4)]; % rolling resistance
                data.contactsWall.localContactPoint(i,:,4) = data.contactsWall.contactPoint(i,:,4) - [x(i) z(i)];
            end
            
            % normal contact
            deltaW(4,i) = radius(i) - abs(z(i)-box(4)); 
            normalConservative_t = normalStiffness*deltaW(4,i);
            normalDissipative_t = par.dampN*2*sqrt(0.5*data.mass(i)*normalStiffness)*vz(i);
            fwz_t(i) = -normalConservative_t - normalDissipative_t;%

                
            % tangential contact
            tangentialSpring = (data.contactsWall.actuationPoint(i,:,4) - data.contactsWall.contactPoint(i,:,4))*[1;0];
            if(abs(tangentialSpring) < 100*eps)
                fwx_t(i) = 0;
            else
                fwx_t(i) = tangentialStiffness*tangentialSpring + par.dampTwall*2*sqrt(0.5*data.mass(i)*tangentialStiffness)*vx(i);
            end
            
            % sliding friction            
            if(abs(fwx_t(i)) > abs(fwz_t(i))*par.muWall)
                fwx_t(i) = sign(fwx_t(i))*abs(fwz_t(i))*par.muWall;
                data.contactsWall.contactPoint(i,:,4) = data.contactsWall.actuationPoint(i,:,4) + tangentialSpring/2*[0;1]';
                data.contactsWall.localContactPoint(i,:,4) = (data.contactsWall.contactPoint(i,:,4))'-[x(i) z(i)]';

            end
%             disp(["fwx_t(i)",fwx_t(i)])
            if(par.considerRotations)
                       twy(i,4) = -fwx_t(i)*(data.contactsWall.actuationPoint(i,2,4)-data.position(2,i));
            end
           if(par.considerRotations)
                data.contactsWall.rollingDeformation(i,1,4) = data.contactsWall.contactPoint(i,1,4) - data.contactsWall.actuationPoint(i,1,4); % 4.27
                data.contactsWall.accumulatedRollingDeformation(i,1,4) = data.contactsWall.accumulatedRollingDeformation(i,1,4) + data.contactsWall.rollingDeformation(i,1,4);%4.28

                if(abs(data.contactsWall.accumulatedRollingDeformation(i,1,4)) > abs(fwz_t(i))*par.muWall/tangentialStiffness*par.Cr) %4.30
                    data.contactsWall.accumulatedRollingDeformation(i,1,4) = data.contactsWall.accumulatedRollingDeformation(i,1,4)/abs(data.contactsWall.accumulatedRollingDeformation(i,1,4))*abs(fwz_t(i))*(par.muWall/tangentialStiffness)*par.Cr;
                end
                  disp(["twy(i,4) no rolling resistance",twy(i,4)])
                  twy(i,4) = twy(i,4) - tangentialStiffness*data.contactsWall.accumulatedRollingDeformation(i,1,4)*(data.contactsWall.actuationPoint(i,2,4)'-data.position(2,i)); % projection into tangential plane necessary
                  disp(["twy(i,4) with rolling resistance",twy(i,4)])
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
