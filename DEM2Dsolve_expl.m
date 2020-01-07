function [pk,vk,ak,Pk,Vk,data] = DEM2Dsolve_expl(data,par,c)
    
    Pk = [];
    Vk = [];
    N = par.N;
    dt = par.dt;

    pk = zeros(2,N);
    vk = zeros(2,N);

    x = data.position(1,:);
    z = data.position(2,:);

    vx = data.velocity(1,:);
    vz = data.velocity(2,:);

    r = data.radius;
    m = data.mass;
    d = DEM2Ddist(x,z);
    
    % fx = zeros(N,N); fz = zeros(N,N);%
    [fx,fz,ty,data] = DEM2DinteractForce(x,z,vx,vz,d,r,par,data,c);
    [fwx,fwz,twy,data] = DEM2DwallForce(x,z,vx,vz,par,data);

    
    for k=1:N
        if(data.contactsParticle.deactivated(k))
            continue
        else
        % Define F:
        F = DEM2DvectField(vx(k),vz(k),fx(k,:),fz(k,:),fwx(k),fwz(k),m(k),par);
        
        if(par.considerRotations)
            % 2D inertia tensor for spheres around y-axis I = 0.25mr²
            I = 0.25*data.mass(k)*(data.radius(k)^2);
            % data.angular
            data.angular(3,k) = data.angular(3,k) +1/(I)*(sum(ty(k,:)) + sum(twy(k,:)))*dt;
            data.angular(2,k) = data.angular(2,k) + data.angular(3,k)*dt;
            data.angular(1,k) = data.angular(1,k) + data.angular(2,k)*dt ;
        end
        vk(1,k) = vx(k) + F(3)*dt;% +data.angular(3,k)*data.radius(k)*dt;% + vtx(k);
        vk(2,k) = vz(k) + F(4)*dt;% + vtz(k);

        pk(1,k) = x(k) + vk(1,k)*dt;
        pk(2,k) = z(k) + vk(2,k)*dt;

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