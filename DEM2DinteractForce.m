function [fx,fz,torqY,data] = DEM2DinteractForce(x,z,vx,vz,d,r,par,data,c)

N = par.N;
Rsparse = sparse(N,N);
for i = 1:length(c.contacts)
    k = c.contacts(i).a;
    l = c.contacts(i).b;
    Rsparse(k,l) = r(k) + r(l);
    data.delta(k,l) =  Rsparse(k,l) - d(k,l);
end

delta = data.delta;
fx = sparse(N,N); fz = sparse(N,N); % particle interaction forces

torqY = sparse(N,N);
Ft = zeros(2,par.N,par.N);
F = sparse(par.N, par.N);

for i=1:N-1
    for k = i:N
        if(data.contactsParticle.isInitialized(i,k))
            if(delta(i,k) < 0)
                data.contactsParticle.PassiveContactAge(i,k) = data.contactsParticle.PassiveContactAge(i,k) +1;

                if(data.contactsParticle.PassiveContactAge(i,k)> data.contactsParticle.maxContactAge)
                data.contactsParticle.ActiveContactAge(i,k) = 0;
                data.contactsParticle.isInitialized(i,k) = 0;
                data.contactsParticle.actuationPoint(:,i,k) = [0; 0];
                end
            end
        end
    end
                
    I = find( delta(i,(i+1):N) > 0 );
    for j=1:length(I)
        if(data.contactsGlued(i,i+I(j)))
            continue;
        end
        % normal direction
        nx = (x(i)-x(i+I(j)))/d(i,i+I(j));
        nz = (z(i)-z(i+I(j)))/d(i,i+I(j)); 
        if(data.contactsParticle.isInitialized(i,i+I(j)))
            data.contactsParticle.ActiveContactAge(i,i+I(j)) = data.contactsParticle.ActiveContactAge(i,i+I(j)) + 1;
            data.contactsParticle.contactPoint(:,i,i+I(j)) = DEM2Drotation(data.angular(2,i)*2*pi*par.dt)*data.contactsParticle.contactPoint(:,i,i+I(j));
            data.contactsParticle.contactPoint(:,i+I(j),i) = DEM2Drotation(data.angular(2,i+I(j))*2*pi*par.dt)*data.contactsParticle.contactPoint(:,i+I(j),i);
        else
            data.contactsParticle.isInitialized(i,i+I(j)) = true;
            data.contactsParticle.PassiveContactAge(i,i+I(j)) = 0;
            data.contactsParticle.actuationPoint(:,i,i+I(j)) = (r(i)*[x(i+I(j));z(i+I(j))] + r(i+I(j))*[x(i);z(i)])/Rsparse(i,i+I(j));
            data.contactsParticle.contactPoint(:,i,i+I(j)) = data.contactsParticle.actuationPoint(:,i,i+I(j)) - [x(i+I(j));z(i+I(j))];
            data.contactsParticle.contactPoint(:,i+I(j),i) = data.contactsParticle.actuationPoint(:,i,i+I(j)) - [x(i);z(i)];            
        end
        % effective mass
        effMass = data.mass(i)*data.mass(j)/(data.mass(i)+data.mass(j));
        % normal stiffness
        normalStiffness = par.Emodul*(r(i)+r(i+I(j)))/2*pi/2;
        
        relVel_ij = data.velocity(:,i) - data.velocity(:,i+I(j));
        %%%%%%%%%%%%%%%%%%% Coulomb Friction %%%%%%%%%%%%%%%%%%%%%%%
        % tangential direction
        tx = nz;
        tz = -nx;
        % normal interaction 
        % cohesive force
        Fcohesion = 0;
        if(par.cohesion>0)
            Fcohesion = -par.cohesion*((data.radius(i)+data.radius(i+I(j))));
        end
%         FnormalCons = normalStiffness*delta(i,i+I(j));
%         FnormalDiss = -par.dampN*2*sqrt(normalStiffness*effMass)*  norm(relVel_ij.*[nx; nz]);
        F(i,i+I(j)) = normalStiffness*delta(i,i+I(j)) - par.dampN*4*sqrt(normalStiffness*effMass)*  norm(relVel_ij.*[nx; nz]) + Fcohesion ;%ddeltadt(i,i+I(j)); 
        F(i+I(j),i) = - F(i,i+I(j));
        % tangential interaction
        globalContactPoint_1 = [x(i+I(j));z(i+I(j))] + data.contactsParticle.contactPoint(:,i,i+I(j));
        globalContactPoint_2 = [x(i);z(i)] + data.contactsParticle.contactPoint(:,i+I(j),i);% check sign
        tangentialSpring = (globalContactPoint_2 - globalContactPoint_1)'*[tx;tz];
        if(abs(tangentialSpring) > eps)
            tangentialSpringLength = norm(tangentialSpring);
            tangentialStiffness = 1/1.2*normalStiffness;
            scaleSpringLength =  F(i,i+I(j))*par.mu/(tangentialStiffness*tangentialSpringLength);
            if( scaleSpringLength < 1 && tangentialSpringLength > 10*eps)
                %% slipping, sliding occurs
                %disp('slipping occurs')
                tangentialSpring = scaleSpringLength*tangentialSpring;
                globalContactPoint_1 = data.contactsParticle.actuationPoint(:,i,i+I(j)) - 0.5*tangentialSpring;
                globalContactPoint_2 = data.contactsParticle.actuationPoint(:,i,i+I(j)) + 0.5*tangentialSpring;
                data.contactsParticle.contactPoint(:,i,i+I(j)) = globalContactPoint_1 - [x(i+I(j));z(i+I(j))];
                data.contactsParticle.contactPoint(:,i+I(j),i) = globalContactPoint_2 - [x(i);z(i)];    
            end
            % add tangential dissipation force;
            Ft(i,i+I(j)) = -tangentialSpring*tangentialStiffness;

            if(par.considerRotations)
            torque1 = Ft(i,i+I(j))*norm(data.contactsParticle.actuationPoint(:,i,i+I(j))- [x(i); z(i)]);%Rsparse(i,i+I(j))/2; % check Obermayr S 29
            torque2= -Ft(i,i+I(j))*norm(data.contactsParticle.actuationPoint(:,i,i+I(j))- [x(i+I(j)); z(i+I(j))]);
            torqY(i) = torque1;
            torqY(i+I(j)) = torque2;
            
            % rolling resistence CHECK!!! 19.02.2020
             data.contactsParticle.rollingDeformation(:,i,i+I(j)) = data.contactsParticle.contactPoint(:,i,i+I(j)) - data.contactsParticle.actuationPoint(:,i,i+I(j));
             data.contactsParticle.accumulatedRollingDeformation(:,i,i+I(j)) = data.contactsParticle.accumulatedRollingDeformation(:,i,i+I(j)) + data.contactsParticle.rollingDeformation(:,i,i+I(j));

            if(abs(data.contactsParticle.accumulatedRollingDeformation(:,i,i+I(j))) > abs(F(i,i+I(j)))*par.mu/tangentialStiffness*0.1)
                data.contactsParticle.accumulatedRollingDeformation(:,i,i+I(j)) = data.contactsParticle.accumulatedRollingDeformation(:,i,i+I(j))*abs(F(i,i+I(j)))*par.muWall/tangentialStiffness*0.1;
            end

%              disp(["torqY(i) no rolling resistance",torqY(i)])
             torqY(i) = torqY(i) + tangentialStiffness*norm(data.contactsParticle.accumulatedRollingDeformation(:,i,i+I(j)))*data.radius(i); % projection into tangential plane necessary
             torqY(i+I(j))= torqY(i+I(j)) - tangentialStiffness*norm(data.contactsParticle.accumulatedRollingDeformation(:,i,i+I(j)))*data.radius(i+I(j));
%              disp(["torqY(i) with rolling resistance",torqY(i)])
            end
        end

        % Euler Equations in 3D
        % omegadot1 = 1/I_1*(I_3 - I_2)omega2*omega3 = M_1
        % omegadot2 = 1/I_2*(I_1 - I_3)omega3*omega1 = M_2
        % omegadot3 = 1/I_3*(I_2 - I_1)omega1*omega2 = M_3
        %%%%%%%%%%%%%%%%%%% END: Coulomb Friction %%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%% Transformation to global coordinate system %%%%%%%%%%%%%%%%%%%%%%%
        fx(i,i+I(j)) = F(i,i+I(j))*nx + Ft(1,i,i+I(j))*tx; fx(i+I(j),i) = F(i+I(j),i)*nx - Ft(1,i,i+I(j))*tx;
        fz(i,i+I(j)) = F(i,i+I(j))*nz + Ft(2,i,i+I(j))*tz; fz(i+I(j),i) = F(i+I(j),i)*nz - Ft(2,i,i+I(j))*tz;
        
        %%%%%%%%%%%%%%%%%%%%%%% Merge Particles %%%%%%%%%%%%%%%%%%%%%%%%%
        % check if ActiveContactAge is large, relative velocities small
        % --> glue particles together
        if(data.contactsParticle.ActiveContactAge(i,i+I(j)) > 10 &&  norm(data.velocity(:,i)-data.velocity(:,i+I(j))) < par.mergeThreashold && (data.contactsGlued(i,i+I(j)) == false) && par.merge == true ...
                && sum(data.contactsGlued(i,:)) == false && sum(data.contactsGlued(:,i+I(j))) == false )
           data.contactsGlued(i,i+I(j)) = true;
           data.contactsGlued(i+I(j),i) = true;
           
            fx(i,i+I(j)) = 0;             fx(i+I(j),i) = 0;
            fz(i,i+I(j)) = 0;             fz(i+I(j),i) = 0;
            data.contactsParticle.mergedParticles = true;
            % deactivate glued particles
            data.contactsParticle.deactivated(i) = 1;
            data.contactsParticle.deactivated(i+I(j)) = 1;
            % append merged particles
            data.contactsMerged.N = [data.contactsMerged.N + 1]; nMerged = data.contactsMerged.N ;
            data.contactsMerged.timeFlag(nMerged) = true;
            data.contactsMerged.index(:,nMerged) = [i ;i+I(j)];
            data.contactsMerged.position(:,nMerged) = [data.position(:,i); data.position(:,i+I(j))];
            data.contactsMerged.mass(nMerged) = (data.mass(i)+data.mass(j)); mergedMass = data.contactsMerged.mass(nMerged);
            data.contactsMerged.positionMerged(:,nMerged) = 1/mergedMass*(data.position(:,i)*data.mass(i)+data.position(:,i+I(j))*data.mass(i+I(j)));
            data.contactsMerged.velocityMerged(:,nMerged) =1/mergedMass*(data.velocity(:,i)*data.mass(i)+data.velocity(:,i+I(j))*data.mass(i+I(j)));
            data.contactsMerged.relativePosition(:,nMerged) = [data.position(:,i)-data.contactsMerged.positionMerged(:,nMerged); data.position(:,i+I(j))-data.contactsMerged.positionMerged(:,nMerged)];
            
        end       
    end 
end

end