function [fx,fz,torqY,data] = DEM2DinteractForce(par,data,c)
    r = data.radius;
    d = DEM2Ddist(data.position(1,:),data.position(2,:));
    N = par.N;
    Rsparse = sparse(N,N);
    for i = 1:length(c.contacts)
        k = c.contacts(i).a;
        if(k > 0)
            l = c.contacts(i).b;
            Rsparse(k,l) = r(k) + r(l);
            data.delta(k,l) =  Rsparse(k,l) - d(k,l);
        else
            continue;
        end
    end

    fx = sparse(N,N); fz = sparse(N,N); % particle interaction forces

    torqY = spalloc(N,N,2*N);
    Ft = zeros(2,N,N);
    F = sparse(N,N);
    frictionForce = zeros(2,N,N);
%for i=1:N-1
for j = 1:c.numParticleContacts
    i = c.contacts(j).a;
    k = c.contacts(j).b;
        if(data.contactsParticle.isInitialized(i,k))
            if(data.delta(i,k) < 0)
                data.contactsParticle.PassiveContactAge(i,k) = data.contactsParticle.PassiveContactAge(i,k) +1;
                if(data.contactsParticle.PassiveContactAge(i,k) > data.contactsParticle.maxContactAge)
                data.contactsParticle.ActiveContactAge(i,k) = 0;
                data.contactsParticle.isInitialized(i,k) = 0;
                data.contactsParticle.actuationPoint(:,i,k) = [0; 0];
                data.contactsParticle.accumulatedRollingDeformation(:,i,k) = [0;0];
                end
            end
    end
    normal = c.contacts(j).n;
    tangential = c.contacts(j).t;            
    %I = find( data.delta(i,(i+1):N) > 0 );
    %for j=1:length(I)
    if(data.delta(i,k) <0)
        continue;
    end
        if(data.contactsGlued(i,k))
            continue;
        end
        % normal direction
        data.contactsParticle.actuationPoint(:,i,k) = (r(i)*data.position(:,k) + r(k)*data.position(:,i))/Rsparse(i,k);
        data.contactsParticle.actuationPoint(:,k,i) = data.contactsParticle.actuationPoint(:,i,k); % redundant
        if(data.contactsParticle.isInitialized(i,k))
            data.contactsParticle.ActiveContactAge(i,k) = data.contactsParticle.ActiveContactAge(i,k) + 1;
           % data.contactsParticle.localContactPoint(:,i,k) = DEM2Drotation(data.angular(1,i))*data.contactsParticle.localContactPoint(:,i,k);
            data.contactsParticle.localContactPoint(:,i,k) = DEM2Drotation(data.angular(2,i)*par.dt)*data.contactsParticle.localContactPoint(:,i,k);
            data.contactsParticle.localContactPoint(:,k,i) = DEM2Drotation(data.angular(2,k)*par.dt)*data.contactsParticle.localContactPoint(:,k,i);
        else
            data.contactsParticle.isInitialized(i,k) = true;
            data.contactsParticle.PassiveContactAge(i,k) = 0;
            data.contactsParticle.localContactPoint(:,i,k) = data.contactsParticle.actuationPoint(:,i,k) - data.position(:,k);
            data.contactsParticle.localContactPoint(:,k,i) = data.contactsParticle.actuationPoint(:,i,k) - data.position(:,i);            
        end
        % effective mass
        effMass = data.mass(i)*data.mass(k)/(data.mass(i)+data.mass(k));
        % normal stiffness
        normalStiffness = par.Emodul*(r(i)+r(k))/2*pi/2;

        relVel_ik = data.velocity(:,i) - data.velocity(:,k);


        % cohesive force
        Fcohesion = 0;
        if(par.cohesion>0)
            Fcohesion = -par.cohesion*pi*((data.radius(i)+data.radius(k))/2)^2;
        end
        F_cons = normalStiffness*data.delta(i,k);
        F_diss = - par.dampN*2*sqrt(normalStiffness*effMass)*(relVel_ik'*normal);
        % normal interaction 
        F(i,k) = F_cons + F_diss + Fcohesion ;%ddeltadt(i,k); 
        F(k,i) = - F(i,k);
        % tangential interaction
        data.contactsParticle.globalContactPoint(:,i,k) = data.position(:,k) + data.contactsParticle.localContactPoint(:,i,k);
        data.contactsParticle.globalContactPoint(:,k,i) = data.position(:,i) + data.contactsParticle.localContactPoint(:,k,i);% check sign
        tangentialSpring = (data.contactsParticle.globalContactPoint(:,k,i) - data.contactsParticle.globalContactPoint(:,i,k));
        
        tangentialSpringLength = norm(tangentialSpring);
        tangentialStiffness = 1/1.2*par.Emodul*(r(i)+r(k))/2*pi/2;%par.kT;%
%%%%%%%%%%%%%%%%%%% Coulomb Friction %%%%%%%%%%%%%%%%%%%%%%%

        scaleSpringLength =  F(i,k)*par.mu/(tangentialStiffness*tangentialSpringLength+eps);
        isSticking = true;
        if( scaleSpringLength < 1 && tangentialSpringLength > 10*eps)
            %% slipping, sliding occurs
            %disp('slipping occurs')
            isSticking = false;
            tangentialSpring = scaleSpringLength*tangentialSpring;
            data.contactsParticle.globalContactPoint(:,i,k) = data.contactsParticle.actuationPoint(:,i,k) - 0.5*tangentialSpring;
            data.contactsParticle.globalContactPoint(:,k,i) = data.contactsParticle.actuationPoint(:,k,i) + 0.5*tangentialSpring;
            data.contactsParticle.localContactPoint(:,i,k) = data.contactsParticle.globalContactPoint(:,i,k) - data.position(:,k);
            data.contactsParticle.localContactPoint(:,k,i) = data.contactsParticle.globalContactPoint(:,k,i) - data.position(:,i);    
        end
        % add tangential dissipation force;
        Ft(:,i,k) =  - tangentialSpring*tangentialStiffness;
        Ft(:,k,i) =  + tangentialSpring*tangentialStiffness;
        %Ft_cons = tangentialSpring*tangentialStiffness%
%         disp(['Ft(1,k,i) no damping: ',num2str(Ft(1,k,i))])
        if(isSticking && par.dampTwall > 0)
               relativeTangentialVelocity = relVel_ik - relVel_ik'*normal*normal;
               %Ft_diss = par.dampT*2*sqrt(effMass*tangentialStiffness)*relativeTangentialVelocity
               Ft(:,i,k) = Ft(:,i,k) - par.dampT*2*sqrt(effMass*tangentialStiffness)*relativeTangentialVelocity;
               Ft(:,k,i) = Ft(:,k,i) + par.dampT*2*sqrt(effMass*tangentialStiffness)*relativeTangentialVelocity;
%                disp(['Ft(1,k,i) with damping: ',num2str(Ft(1,k,i))])
               
        end
        if(par.considerRotations)
            torqY(i,k) = det([(data.contactsParticle.actuationPoint(:,i,k)-data.position(:,i)) Ft(:,i,k)]);%Rsparse(i,k)/2; % check Obermayr S 29
            torqY(k,i) = det([(data.contactsParticle.actuationPoint(:,k,i)-data.position(:,k)) Ft(:,k,i)]);
            % Rolling resistence: 13.05.2020, compare to DEM2DwallForce
            data.contactsParticle.rollingDeformation(:,i,k) = (data.contactsParticle.globalContactPoint(:,k,i) + data.contactsParticle.globalContactPoint(:,i,k))/2 - data.contactsParticle.actuationPoint(:,i,k);% - ((globalContactPoint_1 - globalContactPoint_2)'*[nx;nz])*[nx;nz];
            data.contactsParticle.accumulatedRollingDeformation(:,i,k) = data.contactsParticle.accumulatedRollingDeformation(:,i,k) + data.contactsParticle.rollingDeformation(:,i,k); %rolling deformation
            data.contactsParticle.accumulatedRollingDeformation(:,i,k) = data.contactsParticle.accumulatedRollingDeformation(:,i,k) - (data.contactsParticle.accumulatedRollingDeformation(:,i,k)'*normal)*normal; %projection
            if(norm(data.contactsParticle.accumulatedRollingDeformation(:,i,k)) > abs(F(i,k))*par.mu/tangentialStiffness*par.Cr)
                data.contactsParticle.accumulatedRollingDeformation(:,i,k) = (data.contactsParticle.accumulatedRollingDeformation(:,i,k))/norm(data.contactsParticle.accumulatedRollingDeformation(:,i,k))*(F(i,k))*par.muWall/tangentialStiffness*par.Cr;
            end
% 
           % disp(['torqY(i,k) no rolling resistance ',num2str(torqY(i,k))])
            rollingResistance1 = det([(data.contactsParticle.actuationPoint(:,i,k)-data.position(:,i)) tangentialStiffness*data.contactsParticle.accumulatedRollingDeformation(:,i,k)]);
            rollingResistance2 = det([(data.contactsParticle.actuationPoint(:,i,k)-data.position(:,k)) tangentialStiffness*data.contactsParticle.accumulatedRollingDeformation(:,i,k)]);
            torqY(i,k) = torqY(i,k) + rollingResistance1; %2 d cross product
            torqY(k,i) = torqY(i,k) + rollingResistance2;
           % disp(['torqY(i,k) with rolling resistance ',num2str(torqY(i,k))])


        end


        % Euler Equations in 3D
        % omegadot1 = 1/I_1*(I_3 - I_2)omega2*omega3 = M_1
        % omegadot2 = 1/I_2*(I_1 - I_3)omega3*omega1 = M_2
        % omegadot3 = 1/I_3*(I_2 - I_1)omega1*omega2 = M_3
        %%%%%%%%%%%%%%%%%%% END: Coulomb Friction %%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%% Transformation to global coordinate system %%%%%%%%%%%%%%%%%%%%%%%
        f(:,i,k)= F(i,k)*normal + Ft(:,i,k);
        f(:,k,i)= F(k,i)*normal + Ft(:,k,i);   
        fx(i,k) = f(1,i,k); fz(i,k) = f(2,i,k);
        fx(k,i) = f(1,k,i); fz(k,i) = f(2,k,i);
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