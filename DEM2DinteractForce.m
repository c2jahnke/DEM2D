function [fx,fz] = DEM2DinteractForce(x,z,vx,vz,d,r,par,data)

N = par.N;

R = zeros(N,N);

for i=1:N
    for j=i:N
        R(i,j) = r(i) + r(j);
        R(j,i) = R(i,j); % not needed, drop for optimization
        data.delta(i,j) =  R(i,j) - d(i,j);
        data.delta(j,i) = data.delta(i,j); % not needed, drop for optimization
    end 
end
delta = data.delta;

fx = zeros(N,N);
fz = zeros(N,N);
ddeltadt = zeros(N,N); % calculate ddelta/dt'
for i = 1:N
    for j = i:N
        ddeltadt(i,j) = (delta(i,j)-data.deltaOld(i,j))/par.dt; % du = (u^n+1 - u^n)/dt
        ddeltadt(j,i) = ddeltadt(i,j);
    end
end
torqX = zeros(N,N); torqZ = zeros(N,N);
for i=1:N-1
    I = find( delta(i,(i+1):N) > 0 );
    for j=1:length(I)
        % normal direction
        nx = (x(i)-x(i+I(j)))/d(i,i+I(j));
        nz = (z(i)-z(i+I(j)))/d(i,i+I(j));
        %%%%%%%%%%%%%%%%%%% Coulomb Friction %%%%%%%%%%%%%%%%%%%%%%%
        % tangential direction
        tx = nz;
        tz = -nx;
        % tangential interaction in local coordinates
        if(isequal(data.Xc(:,i,i+I(j)),[0;0]) )
        data.Xc(:,i,i+I(j)) = [tx;tz]*r(i);
        data.Xc(:,i+I(j),i) = [tx;tz]*r(i+I(j));
        end
        xc(:,i,i+I(j)) = DEM2Drotation(data.angular(1,i))*data.Xc(:,i,i+I(j)) + data.position(:,i);
        xc(:,i+I(j),i) = DEM2Drotation(data.angular(1,i))*data.Xc(:,i+I(j),i) + data.position(:,i);
        Xi(:,i,i+I(j)) = xc(:,i,i+I(j)) -xc(:,i+I(j),i);
        data.Xinorm(:,i,i+I(j)) = Xi(:,i,i+I(j)) - (Xi(:,i,i+I(j))'*[nx;nz])*[nx;nz];
        % tangential forces
        Ft(:,i,i+I(j)) = - par.kT*data.Xinorm(:,i,i+I(j)) -par.dampT*(data.Xinorm(:,i,i+I(j)) - data.XinormOld(:,i,i+I(j)))/par.dt;
        data.Xinorm = data.XinormOld;
        % calculate actuation point
        xa(:,i,i+I(j)) = data.position(:,i) + r(i)/(r(i) + r(i+I(j)))*(data.position(:,i+I(j)) - data.position(:,i));
        torqX(i,i+I(j)) = (xa(1,i,i+I(j)) -data.position(1,i)).*Ft(1,i,i+I(j)); % wrong? In 3D it's a cross product
        torqZ(i,i+I(j)) = (xa(2,i,i+I(j)) -data.position(2,i)).*Ft(2,i,i+I(j)); % wrong? In 3D it's a cross product
        %TODO implement 2d euler equation
        %%%%%%%%%%%%%%%%%%% END: Coulomb Friction %%%%%%%%%%%%%%%%%%%%%%%
        % normal interaction
        F(i,i+I(j)) = 31.4*(r(i)+r(i+I(j)))/2*delta(i,i+I(j)) + par.dampN * ddeltadt(i,i+I(j)); % consider mass? par.kN = (r(i)+r(i+I(j)))/2
        F(i+I(j),i) = - F(i,i+I(j));
        %% Warning: check tangential forces
        fx(i,i+I(j)) = F(i,i+I(j))*nx + Ft(1,i,i+I(j))*tx; fx(i+I(j),i) = F(i+I(j),i)*nx - Ft(1,i,i+I(j))*tx;
        fz(i,i+I(j)) = F(i,i+I(j))*nz + Ft(2,i,i+I(j))*tz; fz(i+I(j),i) = F(i+I(j),i)*nz - Ft(2,i,i+I(j))*tz;
        % F_Nij = k_N * delta_ij + d_N ddelta_ij/dt
    end 
end

end