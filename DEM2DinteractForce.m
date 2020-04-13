function [fx,fz,torqY,data] = DEM2DinteractForce(d,r,par,data,c)
    x = data.position(1,:);
    z = data.position(2,:);
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

torqY = spalloc(par.N,par.N,2*par.N);
%Ft = zeros(2,par.N,par.N);
Ft = sparse(N,N);
F = sparse(N, N);

%for i=1:N-1
for j = 1:length(c.contacts)
    i = c.contacts(j).a;
    k = c.contacts(j).b;

    %for k = i:N
        if(data.contactsParticle.isInitialized(i,k))
            if(data.delta(i,k) < 0)
                data.contactsParticle.PassiveContactAge(i,k) = data.contactsParticle.PassiveContactAge(i,k) +1;

                if(data.contactsParticle.PassiveContactAge(i,k)> data.contactsParticle.maxContactAge)
                data.contactsParticle.ActiveContactAge(i,k) = 0;
                data.contactsParticle.isInitialized(i,k) = 0;
                data.contactsParticle.actuationPoint(:,i,k) = [0; 0];
                end
            end
        %end
    end
                
    %I = find( data.delta(i,(i+1):N) > 0 );
    %for j=1:length(I)
    if(data.delta(i,k) <0)
        continue;
    end
        if(data.contactsGlued(i,k))
            continue;
        end
        % normal direction
        nx = (x(i)-x(k))/d(i,k);
        nz = (z(i)-z(k))/d(i,k); 
        if(data.contactsParticle.isInitialized(i,k))
            data.contactsParticle.ActiveContactAge(i,k) = data.contactsParticle.ActiveContactAge(i,k) + 1;
            data.contactsParticle.localContactPoint(:,i,k) = DEM2Drotation(data.angular(2,i)*par.dt)*data.contactsParticle.localContactPoint(:,i,k);
            data.contactsParticle.localContactPoint(:,k,i) = DEM2Drotation(data.angular(2,k)*par.dt)*data.contactsParticle.localContactPoint(:,k,i);
        else
            data.contactsParticle.isInitialized(i,k) = true;
            data.contactsParticle.PassiveContactAge(i,k) = 0;
            data.contactsParticle.actuationPoint(:,i,k) = (r(i)*[x(k);z(k)] + r(k)*[x(i);z(i)])/Rsparse(i,k);
            data.contactsParticle.actuationPoint(:,k,i) = data.contactsParticle.actuationPoint(:,i,k); % redundant
            data.contactsParticle.localContactPoint(:,i,k) = data.contactsParticle.actuationPoint(:,i,k) - [x(k);z(k)];
            data.contactsParticle.localContactPoint(:,k,i) = data.contactsParticle.actuationPoint(:,i,k) - [x(i);z(i)];            
        end
        % effective mass
        effMass = data.mass(i)*data.mass(k)/(data.mass(i)+data.mass(k));
        % normal stiffness
        normalStiffness = par.Emodul*(r(i)+r(k))/2*pi/2;
        relVel_ik = data.velocity(:,i) - data.velocity(:,k);
        %%%%%%%%%%%%%%%%%%% Coulomb Friction %%%%%%%%%%%%%%%%%%%%%%%
        % tangential direction
        tx = nz;
        tz = -nx;
        % normal interaction 
        % cohesive force
        Fcohesion = 0;
        if(par.cohesion>0)
            Fcohesion = -par.cohesion*pi*((data.radius(i)+data.radius(k))/2)^2;
        end
%         FnormalCons = normalStiffness*delta(i,k);
%         FnormalDiss = -par.dampN*2*sqrt(normalStiffness*effMass)*  norm(relVel_ij.*[nx; nz]);
        F(i,k) = normalStiffness*delta(i,k) - par.dampN*2*sqrt(normalStiffness*effMass)*(relVel_ik'*[nx; nz]) + Fcohesion ;%ddeltadt(i,k); 
        F(k,i) = - F(i,k);
%         disp(["F(k,i)",num2str(F(k,i))])
        % tangential interaction
        globalContactPoint_1 = [x(k);z(k)] + data.contactsParticle.localContactPoint(:,i,k);
        globalContactPoint_2 = [x(i);z(i)] + data.contactsParticle.localContactPoint(:,k,i);% check sign
        tangentialSpring = (globalContactPoint_2 - globalContactPoint_1)'*[tx;tz];
        if(abs(tangentialSpring) > eps)
            tangentialSpringLength = norm(tangentialSpring);
            tangentialStiffness = 1/1.2*normalStiffness;
            scaleSpringLength =  F(i,k)*par.mu/(tangentialStiffness*tangentialSpringLength);
            if( scaleSpringLength < 1 && tangentialSpringLength > 10*eps)
                %% slipping, sliding occurs
                %disp('slipping occurs')
                tangentialSpring = scaleSpringLength*tangentialSpring;
                globalContactPoint_1 = data.contactsParticle.actuationPoint(:,i,k) - 0.5*tangentialSpring;
                globalContactPoint_2 = data.contactsParticle.actuationPoint(:,k,i) + 0.5*tangentialSpring;
                data.contactsParticle.localContactPoint(:,i,k) = globalContactPoint_1 - [x(k);z(k)];
                data.contactsParticle.localContactPoint(:,k,i) = globalContactPoint_2 - [x(i);z(i)];    
            end
            % add tangential dissipation force;
%             Ft(i,k) = -tangentialSpring*tangentialStiffness;
            Ft(i,k) =  -tangentialSpring*tangentialStiffness;
            Ft(k,i) =  tangentialSpring*tangentialStiffness;
            if(par.considerRotations)
                torqY(i,k) = Ft(i,k)*(data.contactsParticle.actuationPoint(:,i,k)- [x(i); z(i)])'*[nx; nz];%Rsparse(i,k)/2; % check Obermayr S 29
                torqY(k,i) = -Ft(k,i)*(data.contactsParticle.actuationPoint(:,k,i)- [x(k); z(k)])'*[nx; nz];


                % rolling resistence CHECK!!! 11.03.2020, compare to DEM2DwallForce
                % consider 1-dimensional rollingDeformation in tangential plane
%                 data.contactsParticle.rollingDeformation(:,i,k) = (globalContactPoint_1 - globalContactPoint_2);% - ((globalContactPoint_1 - globalContactPoint_2)'*[nx;nz])*[nx;nz];
%                 data.contactsParticle.accumulatedRollingDeformation(:,i,k) = data.contactsParticle.accumulatedRollingDeformation(:,i,k) + data.contactsParticle.rollingDeformation(:,i,k); %rolling deformation
%                 data.contactsParticle.accumulatedRollingDeformation(:,i,k) = data.contactsParticle.accumulatedRollingDeformation(:,i,k) - (data.contactsParticle.accumulatedRollingDeformation(:,i,k)'*[nx;nz])*[nx;nz]; %projection
%                 if(abs(data.contactsParticle.accumulatedRollingDeformation(:,i,k)) > abs(F(i,k))*par.mu/tangentialStiffness*par.Cr)
%                     data.contactsParticle.accumulatedRollingDeformation(:,i,k) = sign(data.contactsParticle.accumulatedRollingDeformation(:,i,k))*abs(F(i,k))*par.muWall/tangentialStiffness*0.1;
%                 end
% 
%                 disp(['torqY(i,k) no rolling resistance ',num2str(torqY(i,k))])
%                 rollingResistanceiiIj = - tangentialStiffness*det([data.contactsParticle.accumulatedRollingDeformation(:,i,k) (data.contactsParticle.actuationPoint(:,i,k)-[x(i); z(i)])]);
%                 torqY(i,k) = torqY(i,k) - rollingResistanceiiIj; %2 d cross product
%                 disp(['torqY(i,k) with rolling resistance ',num2str(torqY(i,k))])
%                 torqY(k,i) = torqY(k,i) + tangentialStiffness*det([data.contactsParticle.accumulatedRollingDeformation(:,i,k) (data.contactsParticle.actuationPoint(:,k,i)-[x(k); z(k)])]);
                
            end
        end

        % Euler Equations in 3D
        % omegadot1 = 1/I_1*(I_3 - I_2)omega2*omega3 = M_1
        % omegadot2 = 1/I_2*(I_1 - I_3)omega3*omega1 = M_2
        % omegadot3 = 1/I_3*(I_2 - I_1)omega1*omega2 = M_3
        %%%%%%%%%%%%%%%%%%% END: Coulomb Friction %%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%% Transformation to global coordinate system %%%%%%%%%%%%%%%%%%%%%%%
%         fxOld(i,k) = F(i,k)*nx + Ft(1,i,k)*tx; fxOld(k,i) = F(k,i)*nx - Ft(1,i,k)*tx;
%         fzOld(i,k) = F(i,k)*nz + Ft(2,i,k)*tz; fzOld(k,i) = F(k,i)*nz - Ft(2,i,k)*tz;
        fx(i,k) = F(i,k)*nx + Ft(i,k)*tx; fx(k,i) = F(k,i)*nx - Ft(i,k)*tx;
        fz(i,k) = F(i,k)*nz + Ft(i,k)*tz; fz(k,i) = F(k,i)*nz - Ft(i,k)*tz;
        
%         if ~((fx - fxOld) == sparse(N,N) )
%             disp(['warning'])
%         end%&& (fz - fzNew) == sparse(N,N)
        %%%%%%%%%%%%%%%%%%%%%%% Merge Particles %%%%%%%%%%%%%%%%%%%%%%%%%
        % check if ActiveContactAge is large, relative velocities small
        % --> glue particles together
        if(data.contactsParticle.ActiveContactAge(i,k) > 10 &&  norm(data.velocity(:,i)-data.velocity(:,k)) < par.mergeThreashold && (data.contactsGlued(i,k) == false) && par.merge == true ...
                && sum(data.contactsGlued(i,:)) == false && sum(data.contactsGlued(:,k)) == false )
            data.contactsGlued(i,k) = true;
            data.contactsGlued(k,i) = true;
            fx(i,k) = 0;             fx(k,i) = 0;
            fz(i,k) = 0;             fz(k,i) = 0;
            data.contactsParticle.mergedParticles = true;
            % deactivate glued particles
            data.contactsParticle.deactivated(i) = 1;
            data.contactsParticle.deactivated(k) = 1;
            % append merged particles
            data.contactsMerged.N = [data.contactsMerged.N + 1]; nMerged = data.contactsMerged.N ;
            data.contactsMerged.timeFlag(nMerged) = true;
            data.contactsMerged.index(:,nMerged) = [i ;k];
            data.contactsMerged.position(:,nMerged) = [data.position(:,i); data.position(:,k)];
            data.contactsMerged.mass(nMerged) = (data.mass(i)+data.mass(k)); mergedMass = data.contactsMerged.mass(nMerged);
            data.contactsMerged.positionMerged(:,nMerged) = 1/mergedMass*(data.position(:,i)*data.mass(i)+data.position(:,k)*data.mass(k));
            data.contactsMerged.velocityMerged(:,nMerged) =1/mergedMass*(data.velocity(:,i)*data.mass(i)+data.velocity(:,k)*data.mass(k));
            data.contactsMerged.relativePosition(:,nMerged) = [data.position(:,i)-data.contactsMerged.positionMerged(:,nMerged); data.position(:,k)-data.contactsMerged.positionMerged(:,nMerged)];
            
        end       
    %end 
end

end