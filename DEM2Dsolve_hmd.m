function [pk,vk,ak,acc,data] = DEM2Dsolve_hmd(par,data,c)
    % Hertz Mindlin Deresievicz model (1882-1953), similar to EDEM-version
    % (beta stage)
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
    [fx,fz,ty,data] = DEM2DinteractHMD(par,data,c);
    [fwx,fwz,twy,data] = DEM2DwallForce(par,data,c);
    for k=1:N
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
        pk(1,k) = data.position(1,k) + data.velocity(1,k)*dt; 
        pk(2,k) = data.position(2,k) + data.velocity(2,k)*dt;
        % symplectic Euler
        acc(1,k) = ax; acc(2,k) = az;
        vk(1,k) = vx(k) + ax*dt;
        vk(2,k) = vz(k) + az*dt;
        
%         pk(1,k) = data.position(1,k) + vk(1,k)*dt; 
%         pk(2,k) = data.position(2,k) + vk(2,k)*dt;
        ak = data.angular;
        end
end