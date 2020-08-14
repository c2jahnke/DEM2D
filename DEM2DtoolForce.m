function [ftx,ftz,tty,data] = DEM2DtoolForce(par,data,c)
    ftx = zeros(par.N,1);
    ftz = zeros(par.N,1);
    tty = zeros(par.N,4); box = data.toolbBox;
    ftnormal = zeros(par.N,1);     deltaT = zeros(4,par.N);
    fttangential = zeros(par.N,1);
    frictionForce = zeros(2,par.N);
    boundaryContacts = 0;
    Rsparse = sparse(4,par.N);
    for l = 1:length(c.contacts)
        k = -c.contacts(l).a -4;
        if(k > 0)           
            i = c.contacts(l).b;
            if(k< 3)
            data.contactsTool.actuationPoint(i,:,k) = [box(k) data.position(2,i)]; % 1
            deltaT(k,i) = data.radius(i) - abs(box(k)-data.position(1,i));
            elseif(k <5)   
            data.contactsTool.actuationPoint(i,:,k) = [data.position(1,i) box(k)]; % 1
            deltaT(k,i) = data.radius(i) - abs(box(k)-data.position(2,i));
         else
             continue;
         end
            if(deltaT(k,i) > 0)
            boundaryContacts = boundaryContacts +1;
            normalStiffness = par.Emodul*pi/2*data.radius(i);
            tangentialStiffness = 1/1.2*normalStiffness;
           normal = -c.contacts(l).n;
           tangential = -c.contacts(l).t;
           
           if(data.contactsTool.isInitialized(i,k))
                data.contactsTool.contactAge(i,k) = 1;
           else
                data.contactsTool.isInitialized(i,k) = true;
                data.contactsTool.localContactPoint(i,:,k) = data.contactsTool.actuationPoint(i,:,k) - data.position(:,i)';
                data.contactsTool.localContactPoint2(i,:,k) = data.contactsTool.actuationPoint(i,:,k) - [0,0]; % 2,3
           end
           % normal contact
           %ddeltaW(1,i) = -vx(i);
           normalConservative_l = -normalStiffness*deltaT(k,i);
           normalDissipative_l = par.dampN*2*sqrt(0.5*data.mass(i)*normalStiffness)*data.velocity(:,i)'*normal;
           ftnormal(i) =  normalConservative_l - normalDissipative_l;

           % tangential contact
           %data.contactsTool.globalContactPoint(i,:,1) = data.contactsTool.localContactPoint(i,:,1) + [x(i) z(i)];
           data.contactsTool.globalContactPoint(i,:,k) = data.contactsTool.localContactPoint(i,:,k) + data.position(:,i)';
           data.contactsTool.globalContactPoint2(i,:,k) = data.contactsTool.localContactPoint2(i,:,k) + [0,0];%4,5

           tangentialSpring2 = (data.contactsTool.globalContactPoint(i,:,k) - data.contactsTool.globalContactPoint2(i,:,k))' - ((data.contactsTool.globalContactPoint(i,:,k) - data.contactsTool.globalContactPoint2(i,:,k))*normal)*normal;
           tangentialSpringLength = norm(tangentialSpring2);
           scaleSpringLength = abs(ftnormal(i))*par.mu/(tangentialStiffness*tangentialSpringLength+eps);
           isSticking = true;
           if(scaleSpringLength < 1)
               isSticking = false;
               tangentialSpring2 = tangentialSpring2*scaleSpringLength;
               data.contactsTool.globalContactPoint(i,:,k) = data.contactsTool.actuationPoint(i,:,k) + tangentialSpring2'/2;
                   data.contactsTool.globalContactPoint2(i,:,k) = data.contactsTool.actuationPoint(i,:,k) - tangentialSpring2'/2;
                   data.contactsTool.localContactPoint(i,:,k) = (data.contactsTool.globalContactPoint(i,:,k))'-data.position(:,i);
                   data.contactsTool.localContactPoint2(i,:,k) = data.contactsTool.globalContactPoint2(i,:,k) - [0,0];
           end
           frictionForce(:,i) = -tangentialSpring2*tangentialStiffness;
%            disp(['frictionForce(1,i) no damping: ',num2str(frictionForce(1,i))])
           if(isSticking && par.dampTwall > 0)
               relativeTangentialVelocity = data.velocity(:,i) - data.velocity(:,i)'*normal*normal;
               frictionForce(:,i) = frictionForce(:,i) + par.dampTwall*2*sqrt(0.5*data.mass(i)*tangentialStiffness)*relativeTangentialVelocity;
%                disp(['frictionForce(1,i) with damping: ',num2str(frictionForce(1,i))])
           end
           ft = ftnormal(i)*normal +frictionForce(:,i);
           ftx(i) =ftx(i)+ ft(1); ftz(i) =ftz(i)+ ft(2);
          
           if(par.considerRotations)
               tty(i,k) = tty(i,k) + det([(data.contactsTool.actuationPoint(i,:,k)'-data.position(:,i)) frictionForce(:,i)]);
           end

               if(par.considerRotations)
                    data.contactsTool.rollingDeformation(i,:,k) = (data.contactsTool.globalContactPoint2(i,:,k)+data.contactsTool.globalContactPoint(i,:,k))/2 - data.contactsTool.actuationPoint(i,:,k); % 4.27
                    data.contactsTool.accumulatedRollingDeformation(i,:,k) = data.contactsTool.accumulatedRollingDeformation(i,:,k) + data.contactsTool.rollingDeformation(i,:,k);%4.28
                    data.contactsTool.accumulatedRollingDeformation(i,:,k) =  data.contactsTool.accumulatedRollingDeformation(i,:,k) - ( data.contactsTool.accumulatedRollingDeformation(i,:,k)*normal)*normal';
                    if(norm(data.contactsTool.accumulatedRollingDeformation(i,:,k)) > abs(ftnormal(i))*par.muWall/tangentialStiffness*par.CrWall ) %4.30
                        data.contactsTool.accumulatedRollingDeformation(i,:,k) = data.contactsTool.accumulatedRollingDeformation(i,:,k)/norm(data.contactsTool.accumulatedRollingDeformation(i,:,k))*(ftnormal(i))*(par.muWall/tangentialStiffness)*par.CrWall ;
                    end
%                         disp(["twy(i,k) no rolling resistance",twy(i,k)])
                        tty(i,k) = tty(i,k) - det([(data.contactsTool.actuationPoint(i,:,k)'-data.position(:,i)) tangentialStiffness*data.contactsTool.accumulatedRollingDeformation(i,:,k)']); 
%                         disp(["twy(i,k) with rolling resistance",twy(i,k)])
               end
           elseif(data.contactsTool.isInitialized(i,k)) % initialized but no contact with left wall
           data.contactsTool.contactAge(i,k) = data.contactsTool.contactAge(i,k) + 1;
           %data.contactsTool.globalContactPoint(i,:,1) = (data.contactsTool.localContactPoint(i,:,1))'+ [x(i) z(i)]';
           if(data.contactsTool.contactAge(i,k) >= data.contactsTool.maxContactAge)
               data.contactsTool.contactAge(i,k) = 0;
               data.contactsTool.isInitialized(i,k) = false;
               data.contactsTool.actuationPoint(i,:,k) = [0 0];
           end
        end
        else
            continue
        end
    end
end
