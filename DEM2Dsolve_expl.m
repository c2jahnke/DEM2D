function [pk,vk] = DEM2Dsolve_expl(data,par)

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
    [fx,fz] = DEM2DinteractForce(x,z,vx,vz,d,r,par,data);
    [fwx,fwz] = DEM2DwallForce(x,z,vx,vz,d,r,par,data);
    
    
    for i=1:N
    % Define F:
    F = DEM2DvectField(vx(i),vz(i),fx(i,:),fz(i,:),fwx(i),fwz(i),m(i),par);
    
    vk(1,i) = vx(i) + F(3)*dt;
    vk(2,i) = vz(i) + F(4)*dt;
    
    pk(1,i) = x(i) + vk(1,i)*dt;
    pk(2,i) = z(i) + vk(2,i)*dt;
        
end
    

end