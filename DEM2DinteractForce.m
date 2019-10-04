function [fx,fz, torqX, torqZ, vtx,vtz] = DEM2DinteractForce(x,z,vx,vz,d,r,par,data)

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
vtx = zeros(1,N);
vtz = zeros(1,N);
ddeltadt = zeros(N,N); % calculate ddelta/dt'
dist = DEM2Ddist(x,z);
for i = 1:N
    for j = i+1:N
        % numerical delta_dot
        ddeltadt(i,j) = (delta(i,j)-data.deltaOld(i,j))/par.dt; % du = (u^n+1 - u^n)/dt
        % analytical delta_dot
        Analytical_ddeltadt(i,j) = (x(i) - x(j))*(vx(i) - vx(j)) + (z(i) - z(j))*(vz(i) - vz(j))/dist(i,j);
        if(Analytical_ddeltadt(i,j) ~= ddeltadt(i,j))
        disp(['Different delta_dot: Analytical:' num2str(Analytical_ddeltadt(i,j)) 'Numerical:'  num2str(ddeltadt(i,j))]);
        pause(0.1)
        end
        ddeltadt(j,i) = ddeltadt(i,j);
    end
end
torqX = zeros(N,N); torqZ = zeros(N,N);
Ft = zeros(2,par.N,par.N);
xc = zeros(2,par.N,par.N);
Xi = zeros(2,par.N,par.N);
xa = zeros(2,par.N,par.N);
F = zeros(par.N, par.N);
for i=1:N-1
    I = find( delta(i,(i+1):N) > 0 );
    for j=1:length(I)
        % normal direction
        nx = (x(i)-x(i+I(j)))/d(i,i+I(j));
        nz = (z(i)-z(i+I(j)))/d(i,i+I(j));
        % relative mass
        effMass = data.mass(i)*data.mass(j)/(data.mass(i)+data.mass(j));
        % normal stiffness
        normalStiffness = par.Emodul*(r(i)+r(i+I(j)))/2*pi/2;
        %%%%%%%%%%%%%%%%%%% Coulomb Friction %%%%%%%%%%%%%%%%%%%%%%%
        % tangential direction
        tx = nz;
        tz = -nx;
        % tangential interaction in local coordinates
        if(isequal(data.Xc(:,i,i+I(j)),[0;0]) )
        data.Xc(:,i,i+I(j)) = [tx;tz]*r(i);
        data.Xc(:,i+I(j),i) = [tx;tz]*r(i+I(j));
        end
        % local coordinates
        xc(:,i,i+I(j)) = DEM2Drotation(data.angular(1,i))*data.Xc(:,i,i+I(j)) + data.position(:,i);
        xc(:,i+I(j),i) = DEM2Drotation(data.angular(1,i))*data.Xc(:,i+I(j),i) + data.position(:,i);
        Xi(:,i,i+I(j)) = xc(:,i,i+I(j)) -xc(:,i+I(j),i);
        data.Xinorm(:,i,i+I(j)) = Xi(:,i,i+I(j)) - (Xi(:,i,i+I(j))'*[nx;nz])*[nx;nz];
        % tangential forces
        if(data.XinormOld(:,i,i+I(j)) == 0) % fix to avoid strange behaviour by damping term 03.08.2018
            Ft(:,i,i+I(j)) = par.kT*data.Xinorm(:,i,i+I(j))
        else
        Ft(:,i,i+I(j)) = par.kT*data.Xinorm(:,i,i+I(j)) +par.dampT*2*sqrt(par.kT*effMass)*(data.Xinorm(:,i,i+I(j)) - data.XinormOld(:,i,i+I(j)))/par.dt;
        end % end fix
        vtx(i) = -2*Ft(1,i,i+I(j))*par.dt;
        vtz(i) = -2*Ft(2,i,i+I(j))*par.dt;
        % normal interaction % do not hard-code! 
        F(i,i+I(j)) = normalStiffness*delta(i,i+I(j)) + par.dampN*2*sqrt(normalStiffness*effMass)* ddeltadt(i,i+I(j)); % kN = pi/2*E(r(i)+r(i+I(j)))/2 consider mass? par.kN = (r(i)+r(i+I(j)))/2
        F(i+I(j),i) = - F(i,i+I(j));
        if (Ft(1,i,i+I(j)) + Ft(2,i,i+I(j)))^(1/2) > par.mu*F(i,i+I(j))
            % calculate actuation point
            xa(:,i,i+I(j)) = data.position(:,i) + r(i)/(r(i) + r(i+I(j)))*(data.position(:,i+I(j)) - data.position(:,i));
            Xi(:,i,i+I(j)) = par.mu*F(i,i+I(j))/par.kT*data.Xinorm(:,i,i+I(j)); % \Xi_T/norm(Xi_T)
            data.Xinorm(:,i,i+I(j)) = Xi(:,i,i+I(j)) - (Xi(:,i,i+I(j))'*[nx;nz])*[nx;nz];
            Ft(:,i,i+I(j)) = - par.kT*data.Xinorm(:,i,i+I(j)) -par.dampT*(data.Xinorm(:,i,i+I(j)) - data.XinormOld(:,i,i+I(j)))/par.dt;
            torqX(i,i+I(j)) = (xa(1,i,i+I(j)) -data.position(1,i)).*Ft(1,i,i+I(j)); % wrong? In 3D it's a cross product
            torqZ(i,i+I(j)) = (xa(2,i,i+I(j)) -data.position(2,i)).*Ft(2,i,i+I(j)); % wrong? In 3D it's a cross product
        end
        %% Warning: check tangential forces
        %%%%%%%%%%%%%%%%%%% END: Coulomb Friction %%%%%%%%%%%%%%%%%%%%%%%
        fx(i,i+I(j)) = F(i,i+I(j))*nx + Ft(1,i,i+I(j))*tx; fx(i+I(j),i) = F(i+I(j),i)*nx - Ft(1,i,i+I(j))*tx;
        fz(i,i+I(j)) = F(i,i+I(j))*nz + Ft(2,i,i+I(j))*tz; fz(i+I(j),i) = F(i+I(j),i)*nz - Ft(2,i,i+I(j))*tz;
        % F_Nij = k_N * delta_ij + d_N ddelta_ij/dt
    end 
end
data.deltaOld = data.delta;

end