function [pk,vk,ak,acc,Pk,Vk,data] = DEM2Dsolve_expl(par,data,c)
% linear contact modell (Obermayr 2011)
    N = par.N;
    Pk = zeros(2,N);
    Vk = zeros(2,N);

    dt = par.dt;

    pk = data.position;
    vk = zeros(2,N);
    acc = zeros(2,N);
    ak = zeros(2,N);
    vx = data.velocity(1,:);
    vz = data.velocity(2,:);

    r = data.radius;
    m = data.mass;
    ftx = zeros(par.N,1); ftz = zeros(par.N,1); tty = zeros(par.N,1);
    [fx,fz,ty,data] = DEM2DinteractForce(par,data,c);
    [fwx,fwz,twy,data] = DEM2DwallForce(par,data,c);
    if(par.toolBool)
        [ftx,ftz,tty,data] = DEM2DtoolForce(par,data,c);
    end
   for k=1:N
        if(data.contactsParticle.deactivated(k))
            continue
        else
        ax = (sum(fx(k,:)) + fwx(k,:) + ftx(k,:))/m(k) + par.g_vert;
        az = (sum(fz(k,:)) + fwz(k,:) + ftz(k,:))/m(k) + par.g;
        if(par.considerRotations)
            % 2D inertia tensor for spheres around y-axis I = 0.25mr^2
            I = 0.25*data.mass(k)*(data.radius(k)^2);
            data.angular(2,k) = data.angular(2,k) + 1/(I)*(sum(ty(k,:)) + sum(twy(k,:)) + sum(tty(k,:)))*dt;%data.angular(3,k)*dt;
            data.angular(1,k) = data.angular(1,k) + data.angular(2,k)*dt ;
            data.contactsWall.localContactPoint(k,:,1) = DEM2Drotation(data.angular(2,k)*dt)*data.contactsWall.localContactPoint(k,:,1)'; 
            data.contactsWall.localContactPoint(k,:,2) = DEM2Drotation(data.angular(2,k)*dt)*data.contactsWall.localContactPoint(k,:,2)'; 
            data.contactsWall.localContactPoint(k,:,3) = DEM2Drotation(data.angular(2,k)*dt)*data.contactsWall.localContactPoint(k,:,3)'; 
            data.contactsWall.localContactPoint(k,:,4) = DEM2Drotation(data.angular(2,k)*dt)*data.contactsWall.localContactPoint(k,:,4)';
            
            data.contactsTool.localContactPoint(k,:,1) = DEM2Drotation(data.angular(2,k)*dt)*data.contactsTool.localContactPoint(k,:,1)'; 
            data.contactsTool.localContactPoint(k,:,2) = DEM2Drotation(data.angular(2,k)*dt)*data.contactsTool.localContactPoint(k,:,2)'; 
            data.contactsTool.localContactPoint(k,:,3) = DEM2Drotation(data.angular(2,k)*dt)*data.contactsTool.localContactPoint(k,:,3)'; 
            data.contactsTool.localContactPoint(k,:,4) = DEM2Drotation(data.angular(2,k)*dt)*data.contactsTool.localContactPoint(k,:,4)';
            
            if(data.contactsWall.isInitialized) % check
                pk(1,k) = pk(1,k);% + data.angular(2,k)*2*pi*dt;
            end
            
        end
        % standard Euler
%         pk(1,k) = data.position(1,k) + data.velocity(1,k)*dt; 
%         pk(2,k) = data.position(2,k) + data.velocity(2,k)*dt;
        % symplectic Euler
        acc(1,k) = ax; acc(2,k) = az;
        vk(1,k) = vx(k) + ax*dt;
        vk(2,k) = vz(k) + az*dt;
        
        pk(1,k) = data.position(1,k) + vk(1,k)*dt; 
        pk(2,k) = data.position(2,k) + vk(2,k)*dt;
        
        data.toolbBox(:,1) = data.toolbBox(:,1) + par.toolSpeed(1)*par.dt;
        data.toolbBox(:,2) = data.toolbBox(:,2) + par.toolSpeed(2)*par.dt;
        ak = data.angular;
        end
    end
    if(data.contactsParticle.mergedParticles)
        for k = 1:data.contactsMerged.N
            i = nonzeros(data.contactsMerged.index(k,:));
            ax = 0; az = 0; tyMerged = 0;
            for jj = 1:data.contactsMerged.aggregateSize(k)
                
            %% angular momentum
              data.contactsMerged.angularMerged(2,k) = data.contactsMerged.angularMerged(2,k) ...
                  + 1/(data.contactsMerged.inertiaTensor(k))*(sum(ty(i(jj),:)) + sum(twy(i(jj),:)))*dt; ...
%              + 0.01*det([data.contactsMerged.relativePosition(k,:)' [sum(fx(i(jj),:))+fwx(i(jj),:);sum(fz(i(jj),:))+fwz(i(jj))]]); 
              ax = ax + sum(fx(i(jj),:)+fwx(i(jj),:))/m(i(jj));
              az = az + sum(fz(i(jj),:)+fwz(i(jj),:))/m(i(jj));
            % angular momentum
            end
      
            data.contactsMerged.angularMerged(1,k) = data.contactsMerged.angularMerged(1,k) + data.contactsMerged.angularMerged(2,k)*dt;
            data.contactsMerged.velocityMerged(:,k) = data.contactsMerged.velocityMerged(:,k) + [ax + par.g_vert; az + par.g]*dt;
            data.contactsMerged.positionMerged(:,k) = data.contactsMerged.positionMerged(:,k) + data.contactsMerged.velocityMerged(:,k)*dt;
            
            for jj = 1:data.contactsMerged.aggregateSize(k)
                data.contactsMerged.relativePosition(2*jj-1:2*jj,k) = DEM2Drotation(data.contactsMerged.angularMerged(2,k)*dt)*data.contactsMerged.relativePosition(2*jj-1:2*jj,k);
                vk(:,i(jj)) = data.contactsMerged.velocityMerged(:,k);% + vtx(k);
                pk(:,i(jj)) = data.contactsMerged.positionMerged(:,k) + data.contactsMerged.relativePosition(2*jj-1:2*jj,k);
            end
        end
        
        Pk = data.contactsMerged.positionMerged;
        Vk = data.contactsMerged.velocityMerged;
    end
    
end