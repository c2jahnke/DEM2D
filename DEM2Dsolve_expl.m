function [pk,vk,ak] = DEM2Dsolve_expl(data,par)

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
    [fx,fz,tx,tz,vtx,vtz] = DEM2DinteractForce(x,z,vx,vz,d,r,par,data);
    [fwx,fwz] = DEM2DwallForce(x,z,vx,vz,d,r,par,data);
    
    
    for k=1:N
    % Define F:
    F = DEM2DvectField(vx(k),vz(k),fx(k,:),fz(k,:),fwx(k),fwz(k),m(k),par);
    
    vk(1,k) = vx(k) + F(3)*dt + vtx(k);
    vk(2,k) = vz(k) + F(4)*dt + vtz(k);
    
    pk(1,k) = x(k) + vk(1,k)*dt;
    pk(2,k) = z(k) + vk(2,k)*dt;
    
    data.angular(1,k) = data.angular(1,k) + data.angular(2,k)*dt;
    data.angular(2,k) = data.angular(2,k) + data.angular(3,k)*dt;
    % data.angular(3,i) = [tx(i); tz(i)]; How to do this?
     
end
    ak = data.angular;
end