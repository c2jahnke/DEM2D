function [pk,vk,ak,acc,Pk,Vk,data] = DEM2Dsolve_expl(par,data,c)
    N = par.N;
    Pk = zeros(2,N);
    Vk = zeros(2,N);

    dt = par.dt;

    pk = zeros(2,N);
    vk = zeros(2,N);
    acc = zeros(2,N);
    ak = zeros(2,N);
    vx = data.velocity(1,:);
    vz = data.velocity(2,:);

    r = data.radius;
    m = data.mass;

    [fx,fz,ty,data] = DEM2DinteractForce(par,data,c);
    [fwx,fwz,twy,data] = DEM2DwallForce(par,data,c);

    
   for k=1:N
        if(data.contactsParticle.deactivated(k))
            continue
        else
        ax = (sum(fx(k,:)) + fwx(k,:))/m(k) + par.g_vert;
        az = (sum(fz(k,:)) + fwz(k,:))/m(k) + par.g;
        if(par.considerRotations)
            % 2D inertia tensor for spheres around y-axis I = 0.25mr^2
            I = 0.25*data.mass(k)*(data.radius(k)^2);
            data.angular(2,k) = data.angular(2,k) + 1/(I)*(sum(ty(k,:)) + sum(twy(k,:)))*dt;%data.angular(3,k)*dt;
            data.angular(1,k) = data.angular(1,k) + data.angular(2,k)*dt ;
            data.contactsWall.localContactPoint(k,:,1) = DEM2Drotation(data.angular(2,k)*dt)*data.contactsWall.localContactPoint(k,:,1)'; 
            data.contactsWall.localContactPoint(k,:,2) = DEM2Drotation(data.angular(2,k)*dt)*data.contactsWall.localContactPoint(k,:,2)'; 
            data.contactsWall.localContactPoint(k,:,3) = DEM2Drotation(data.angular(2,k)*dt)*data.contactsWall.localContactPoint(k,:,3)'; 
            data.contactsWall.localContactPoint(k,:,4) = DEM2Drotation(data.angular(2,k)*dt)*data.contactsWall.localContactPoint(k,:,4)';
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
        ak = data.angular;
        end
    end
    if(data.contactsParticle.mergedParticles)
        for k = 1:data.contactsMerged.N
            i = nonzeros(data.contactsMerged.index(k,:));
            ax = 0; az = 0;
            for jj = 1:data.contactsMerged.aggregateSize(k)
                
            %i = data.contactsMerged.index(k,1); j = data.contactsMerged.index(k,2);
%             data.contactsMerged.inertiaTensor(1,i,j) = 0.5*data.mass(i)*(data.radius(i)^2)...
%                 + data.mass(i)*norm(data.position(:,i) - data.contactsParticle.actuationPoint(:,i,j))^2 ...
%                 + 0.5*data.mass(j)*(data.radius(j)^2)  ...
%                 + data.mass(i)*norm(data.position(:,i) - data.contactsParticle.actuationPoint(:,i,j))^2; % J_i + m_i*d_ij^2 + J_j + m_j*d_ji^2
         data.contactsMerged.inertiaTensor(k) = 0.5*data.mass(i(jj))*(data.radius(i(jj))^2)...
                 + data.mass(i(jj))*norm(data.position(:,i(jj)) - data.contactsMerged.positionMerged(:,k))^2
                
%             data.contactsMerged.angularMerged(2,k) = data.contactsMerged.angularMerged(2,k) + 1/(data.contactsMerged.inertiaTensor(1,i,j))*(sum(ty(i,:)) + sum(twy(j,:)) + sum(ty(i,:)) + sum(twy(j,:)))*dt;
%             data.contactsMerged.angularMerged(1,k) = data.contactsMerged.angularMerged(1,k) + data.contactsMerged.angularMerged(2,k)*dt;
%             ax = (sum(fx(i,:)+fx(j,:)) +fwx(i,:)+fwx(j,:))/m(k) + par.g_vert;
%             az = (sum(fz(i,:)+fz(j,:)) + fwz(i,:)+fwz(j,:))/m(k) + par.g;
%                         % angular momentum
%             
%             vk(:,i) = data.contactsMerged.velocityMerged(:,k);% + vtx(k);
%             vk(:,j) = data.contactsMerged.velocityMerged(:,k);% + vtz(k);
% 
%             pk(:,i) = data.contactsMerged.positionMerged(:,k) + data.contactsMerged.relativePosition(1:2,k);
%             pk(:,j) = data.contactsMerged.positionMerged(:,k) + data.contactsMerged.relativePosition(3:4,k);
            data.contactsMerged.angularMerged(2,k) = data.contactsMerged.angularMerged(2,k) + 1/(data.contactsMerged.inertiaTensor(k))*(sum(ty(i(jj),:)) + sum(twy(i(jj),:)))*dt;
            ax = ax + sum(fx(i(jj),:)+fwx(i(jj),:))/m(i(jj));
            az = az + sum(fz(i(jj),:)+fwz(i(jj),:))/m(i(jj));
                        % angular momentum
            
            end
            data.contactsMerged.angularMerged(1,k) = data.contactsMerged.angularMerged(1,k) + data.contactsMerged.angularMerged(2,k)*dt;
            data.contactsMerged.velocityMerged(:,k) = data.contactsMerged.velocityMerged(:,k) + [ax + par.g_vert; az + par.g]*dt;
            data.contactsMerged.positionMerged(:,k) = data.contactsMerged.positionMerged(:,k) + data.contactsMerged.velocityMerged(:,k)*dt;
            for jj = 1:data.contactsMerged.aggregateSize(k)
                vk(:,i(jj)) = data.contactsMerged.velocityMerged(:,k);% + vtx(k);
                pk(:,i(jj)) = data.contactsMerged.positionMerged(:,k) + data.contactsMerged.relativePosition(2*jj-1:2*jj,k);
            end
        end
        
        Pk = data.contactsMerged.positionMerged;
        Vk = data.contactsMerged.velocityMerged;
    end
    
end