function [fwx,fwz,twy,data] = DEM2DwallForce(par,data,c)
    fwx = zeros(par.N,1);
    fwz = zeros(par.N,1);
    twy = zeros(par.N,4); box = par.bBox;
    fwnormal = zeros(par.N,1);     deltaW = zeros(4,par.N);
    fwtangential = zeros(par.N,1);
    frictionForce = zeros(2,par.N);
    boundaryContacts = 0;
    Rsparse = sparse(4,par.N);
    for l = 1:length(c.contacts)
        k = -c.contacts(l).a;
        if(k > 0)           
            i = c.contacts(l).b;
            if(k< 3)
           data.contactsWall.actuationPoint(i,:,k) = [box(k) data.position(2,i)]; % 1
           deltaW(k,i) = data.radius(i) - abs(box(k)-data.position(1,i));
%            ddeltaW = -data.velocity(1,i);
           elseif(k <5)   
           data.contactsWall.actuationPoint(i,:,k) = [data.position(1,i) box(k)]; % 1
           deltaW(k,i) = data.radius(i) - abs(box(k)-data.position(2,i));
           elseif(k<9)
               continue;
           end
            if(deltaW(k,i) > 0)
            boundaryContacts = boundaryContacts +1;
            normalStiffness = par.Emodul*pi/2*data.radius(i);
            tangentialStiffness = 1/1.2*normalStiffness;
           normal = -c.contacts(l).n;
           tangential = -c.contacts(l).t;
           
           if(data.contactsWall.isInitialized(i,k))
                data.contactsWall.contactAge(i,k) = 1;
           else
                data.contactsWall.isInitialized(i,k) = true;
                data.contactsWall.localContactPoint(i,:,k) = data.contactsWall.actuationPoint(i,:,k) - data.position(:,i)';
                data.contactsWall.localContactPoint2(i,:,k) = data.contactsWall.actuationPoint(i,:,k) - [0,0]; % 2,3
           end
           % normal contact
           %ddeltaW(1,i) = -vx(i);
           normalConservative_l = normalStiffness*deltaW(k,i);
           normalDissipative_l = par.dampN*2*sqrt(0.5*data.mass(i)*normalStiffness)*data.velocity(:,i)'*normal;
           fwnormal(i) =  normalConservative_l - normalDissipative_l;

           % tangential contact
           %data.contactsWall.globalContactPoint(i,:,1) = data.contactsWall.localContactPoint(i,:,1) + [x(i) z(i)];
           data.contactsWall.globalContactPoint(i,:,k) = data.contactsWall.localContactPoint(i,:,k) + data.position(:,i)';
           data.contactsWall.globalContactPoint2(i,:,k) = data.contactsWall.localContactPoint2(i,:,k) + [0,0];%4,5

           tangentialSpring2 = (data.contactsWall.globalContactPoint(i,:,k) - data.contactsWall.globalContactPoint2(i,:,k))' - ((data.contactsWall.globalContactPoint(i,:,k) - data.contactsWall.globalContactPoint2(i,:,k))*normal)*normal;
           tangentialSpringLength = norm(tangentialSpring2);
           scaleSpringLength = abs(fwnormal(i))*par.mu/(tangentialStiffness*tangentialSpringLength+eps);
           isSticking = true;
           if(scaleSpringLength < 1)
               isSticking = false;
               tangentialSpring2 = tangentialSpring2*scaleSpringLength;
               data.contactsWall.globalContactPoint(i,:,k) = data.contactsWall.actuationPoint(i,:,k) + tangentialSpring2'/2;
                   data.contactsWall.globalContactPoint2(i,:,k) = data.contactsWall.actuationPoint(i,:,k) - tangentialSpring2'/2;
                   data.contactsWall.localContactPoint(i,:,k) = (data.contactsWall.globalContactPoint(i,:,k))'-data.position(:,i);
                   data.contactsWall.localContactPoint2(i,:,k) = data.contactsWall.globalContactPoint2(i,:,k) - [0,0];
           end
           frictionForce(:,i) = -tangentialSpring2*tangentialStiffness;
%            disp(['frictionForce(1,i) no damping: ',num2str(frictionForce(1,i))])
           if(isSticking && par.dampTwall > 0)
               relativeTangentialVelocity = data.velocity(:,i) - data.velocity(:,i)'*normal*normal;
               frictionForce(:,i) = frictionForce(:,i) + par.dampTwall*2*sqrt(0.5*data.mass(i)*tangentialStiffness)*relativeTangentialVelocity;
%                disp(['frictionForce(1,i) with damping: ',num2str(frictionForce(1,i))])
           end
           fw = fwnormal(i)*normal +frictionForce(:,i);
           fwx(i) =fwx(i)+ fw(1); fwz(i) =fwz(i)+ fw(2);
          
           if(par.considerRotations)
               twy(i,k) = twy(i,k) + det([(data.contactsWall.actuationPoint(i,:,k)'-data.position(:,i)) frictionForce(:,i)]);
           end

               if(par.considerRotations)
                    data.contactsWall.rollingDeformation(i,:,k) = (data.contactsWall.globalContactPoint2(i,:,k)+data.contactsWall.globalContactPoint(i,:,k))/2 - data.contactsWall.actuationPoint(i,:,k); % 4.27
                    data.contactsWall.accumulatedRollingDeformation(i,:,k) = data.contactsWall.accumulatedRollingDeformation(i,:,k) + data.contactsWall.rollingDeformation(i,:,k);%4.28
                    data.contactsWall.accumulatedRollingDeformation(i,:,k) =  data.contactsWall.accumulatedRollingDeformation(i,:,k) - ( data.contactsWall.accumulatedRollingDeformation(i,:,k)*normal)*normal';
                    if(norm(data.contactsWall.accumulatedRollingDeformation(i,:,k)) > abs(fwnormal(i))*par.muWall/tangentialStiffness*par.CrWall ) %4.30
                        data.contactsWall.accumulatedRollingDeformation(i,:,k) = data.contactsWall.accumulatedRollingDeformation(i,:,k)/norm(data.contactsWall.accumulatedRollingDeformation(i,:,k))*(fwnormal(i))*(par.muWall/tangentialStiffness)*par.CrWall ;
                    end
%                         disp(["twy(i,k) no rolling resistance",twy(i,k)])
                        twy(i,k) = twy(i,k) - det([(data.contactsWall.actuationPoint(i,:,k)'-data.position(:,i)) tangentialStiffness*data.contactsWall.accumulatedRollingDeformation(i,:,k)']); 
%                         disp(["twy(i,k) with rolling resistance",twy(i,k)])
               end
           elseif(data.contactsWall.isInitialized(i,k)) % initialized but no contact with left wall
           data.contactsWall.contactAge(i,k) = data.contactsWall.contactAge(i,k) + 1;
           %data.contactsWall.globalContactPoint(i,:,1) = (data.contactsWall.localContactPoint(i,:,1))'+ [x(i) z(i)]';
           if(data.contactsWall.contactAge(i,k) >= data.contactsWall.maxContactAge)
               data.contactsWall.contactAge(i,k) = 0;
               data.contactsWall.isInitialized(i,k) = false;
               data.contactsWall.actuationPoint(i,:,k) = [0 0];
           end
        end
        else
            continue
        end
    end
end
