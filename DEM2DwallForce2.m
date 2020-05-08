function [fwx,fwz,twy,data] = DEM2DwallForce2(vx,vz,par,data,c)
    fwx = zeros(par.N,1);
    fwz = zeros(par.N,1);
    twy = zeros(par.N,4); box = par.bBox;
    fwnormal = zeros(par.N,1);     deltaW = zeros(4,par.N);
    fwtangential = zeros(par.N,1);
    boundaryContacts = 0;
    Rsparse = sparse(4,par.N);
    for l = 1:length(c.contacts)
        k = -c.contacts(l).a;
        if(k > 0)           
            i = c.contacts(l).b;
            if(k< 3)
%                dirtyRollingHack = 1;
           data.contactsWall.actuationPoint(i,:,k) = [box(k) data.position(2,i)]; % 1
           deltaW(k,i) = data.radius(i) - abs(box(k)-data.position(1,i));
%            ddeltaW = -data.velocity(1,i);
           elseif(k <5)   
%                dirtyRollingHack = -1;
           data.contactsWall.actuationPoint(i,:,k) = [data.position(1,i) box(k)]; % 1
           deltaW(k,i) = data.radius(i) - abs(box(k)-data.position(2,i));
%            ddeltaW = data.velocity(2,i);
           end
            if(deltaW(k,i) > 0)
            boundaryContacts = boundaryContacts +1;
            normalStiffness = par.Emodul*pi/2*data.radius(i);
            tangentialStiffness = 1/1.2*normalStiffness;
           normal = c.contacts(l).n;
           tangential = c.contacts(l).t;
           
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
           tangentialSpring = (data.contactsWall.globalContactPoint(i,:,k) - data.contactsWall.globalContactPoint2(i,:,k))*tangential; %6
           % tangentialSpring = -(data.contactsWall.actuationPoint(i,:,1) - data.contactsWall.globalContactPoint(i,:,1))*[0;1]
           if(abs(tangentialSpring) < 100*eps)
               fwtangential(i) = 0;
           else
               fwtangential(i) = (tangentialStiffness*tangentialSpring - par.dampTwall*2*sqrt(0.5*data.mass(i)*tangentialStiffness)*data.velocity(:,i)'*tangential);
               % friction
               if(abs(fwtangential(i)) > abs(fwnormal(i))*par.muWall)
                   fwtangential(i) = sign(fwtangential(i))*abs(fwnormal(i))*par.muWall;
                   tangentialSpring = fwnormal(i)*par.muWall/tangentialStiffness;
                   data.contactsWall.globalContactPoint(i,:,k) = data.contactsWall.actuationPoint(i,:,k) + tangentialSpring/2*tangential';
                   data.contactsWall.globalContactPoint2(i,:,k) = data.contactsWall.actuationPoint(i,:,k) - tangentialSpring/2*tangential';
                   data.contactsWall.localContactPoint(i,:,k) = (data.contactsWall.globalContactPoint(i,:,k))'-data.position(:,i);
                   data.contactsWall.localContactPoint2(i,:,k) = data.contactsWall.globalContactPoint2(i,:,k) - [0,0];
               end

           end
               fw = fwnormal(i)*normal - fwtangential(i)*tangential;
               fwx(i) =fwx(i)+ fw(1); fwz(i) =fwz(i)+ fw(2);
          
               %DEM2Drotation(data.angular(2,i)*par.dt)*
               if(par.considerRotations)
           %        twy(i,1) = fwz_l(i)*norm(data.contactsWall.actuationPoint(i,:,1)'-data.position(:,i));
                   twy(i,k) = twy(i,k) + det([fwtangential(i)*tangential (data.contactsWall.actuationPoint(i,:,k)'-data.position(:,i))]);
               end
               if(par.considerRotations)
                    data.contactsWall.rollingDeformation(i,1,1) = data.contactsWall.globalContactPoint(i,2,1) - data.contactsWall.actuationPoint(i,2,1); % 4.27
                    data.contactsWall.accumulatedRollingDeformation(i,1,1) = data.contactsWall.accumulatedRollingDeformation(i,1,1) + data.contactsWall.rollingDeformation(i,1,1);%4.28
                    if(abs(data.contactsWall.accumulatedRollingDeformation(i,1,1)) > abs(fwnormal(i))*par.muWall/tangentialStiffness*par.Cr) %4.30
                        data.contactsWall.accumulatedRollingDeformation(i,1,1) = data.contactsWall.accumulatedRollingDeformation(i,1,1)/abs(data.contactsWall.accumulatedRollingDeformation(i,1,1))*abs(fwnormal(i))*(par.muWall/tangentialStiffness)*par.Cr;
                    end
                      %disp(["twy(i,k) no rolling resistance",twy(i,k)])
                      twy(i,k) = twy(i,k) + tangentialStiffness*data.contactsWall.accumulatedRollingDeformation(i,1,1)*(data.contactsWall.actuationPoint(i,1,1)'-data.position(1,i)); % projection into tangential plane necessary
                      %disp(["twy(i,k) with rolling resistance",twy(i,k)])
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
