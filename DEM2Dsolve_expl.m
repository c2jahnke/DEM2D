function [pk,vk,ak,acceleration,Pk,Vk,data] = DEM2Dsolve_expl(data,par,c)
    N = par.N;
    Pk = zeros(2,N);
    Vk = zeros(2,N);

    dt = par.dt;

    pk = zeros(2,N);
    vk = zeros(2,N);
    
    %x = data.position(1,:);
    %z = data.position(2,:);

    vx = data.velocity(1,:);
    vz = data.velocity(2,:);

    r = data.radius;
    m = data.mass;
    d = DEM2Ddist(data.position(1,:),data.position(2,:));
    
    [fx,fz,ty,data] = DEM2DinteractForce(d,r,par,data,c);
    [fwx,fwz,twy,data] = DEM2DwallForce2(vx,vz,par,data,c);

    
   for k=1:N
        if(data.contactsParticle.deactivated(k))
            continue
        else
        ax = (sum(fx(k,:)) + fwx(k,:))/m(k) - par.g_vert;% + par.g;% - par.g;
        az = (sum(fz(k,:)) + fwz(k,:))/m(k) - par.g;
        if(par.considerRotations)
            % 2D inertia tensor for spheres around y-axis I = 0.25mr^2
            I = 0.25*data.mass(k)*(data.radius(k)^2);
            data.angular(2,k) = data.angular(2,k) + 1/(I)*(sum(ty(k,:)) + sum(twy(k,:)))*dt;%data.angular(3,k)*dt;
            data.angular(1,k) = data.angular(1,k) + data.angular(2,k)*dt ;
            data.contactsWall.localContactPoint(k,:,1) = DEM2Drotation(data.angular(2,k)*dt)*data.contactsWall.localContactPoint(k,:,1)'; 
            data.contactsWall.localContactPoint(k,:,2) = DEM2Drotation(data.angular(2,k)*dt)*data.contactsWall.localContactPoint(k,:,2)'; 
            data.contactsWall.localContactPoint(k,:,3) = DEM2Drotation(data.angular(2,k)*dt)*data.contactsWall.localContactPoint(k,:,3)'; 
            data.contactsWall.localContactPoint(k,:,4) = DEM2Drotation(data.angular(2,k)*dt)*data.contactsWall.localContactPoint(k,:,4)';
%             data.contactsWall.globalContactPoint(k,:,1) = (data.contactsWall.localContactPoint(k,:,1))'+ data.position(:,k);
            if(data.contactsWall.isInitialized) % check
                pk(1,k) = pk(1,k);% + data.angular(2,k)*2*pi*dt;
            end
            
        end

%         pk(1,k) = data.position(1,k) + data.velocity(1,k)*dt; 
%         pk(2,k) = data.position(2,k) + data.velocity(2,k)*dt;
        acceleration(1,k) = ax; acceleration(2,k) = az;
        vk(1,k) = vx(k) + ax*dt;% + vtx(k);
        vk(2,k) = vz(k) + az*dt;% + vtz(k);
        
        pk(1,k) = data.position(1,k) + vk(1,k)*dt; 
        pk(2,k) = data.position(2,k) + vk(2,k)*dt;
        ak = data.angular;
        end
    end
    if(data.contactsParticle.mergedParticles)
        for k = 1:data.contactsMerged.N
            
            i = data.contactsMerged.index(1,k); j = data.contactsMerged.index(2,k);

            G = DEM2DvectField(data.contactsMerged.velocityMerged(1,k),data.contactsMerged.velocityMerged(2,k),fx(i,:)+fx(j,:),fz(i,:)+fz(j,:),fwx(i,:)+fwx(j,:),fwz(i,:)+fwz(j,:),data.contactsMerged.mass(k),par);
            
            data.contactsMerged.velocityMerged(1,k) = data.contactsMerged.velocityMerged(1,k) + G(3)*dt;% + vtx(k);
            data.contactsMerged.velocityMerged(2,k) = data.contactsMerged.velocityMerged(2,k) + G(4)*dt;% + vtz(k);

            data.contactsMerged.positionMerged(:,k) = data.contactsMerged.positionMerged(:,k) + data.contactsMerged.velocityMerged(:,k).*dt;  
            
            % angular momentum
            
            vk(:,i) = data.contactsMerged.velocityMerged(:,k);% + vtx(k);
            vk(:,j) = data.contactsMerged.velocityMerged(:,k);% + vtz(k);

            pk(:,i) = data.contactsMerged.positionMerged(:,k) + data.contactsMerged.relativePosition(1:2,k);
            pk(:,j) = data.contactsMerged.positionMerged(:,k) + data.contactsMerged.relativePosition(3:4,k);
        end
        Pk = data.contactsMerged.positionMerged;
        Vk = data.contactsMerged.velocityMerged;

    end
    
end