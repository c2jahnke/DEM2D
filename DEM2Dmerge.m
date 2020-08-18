function [par,data] = DEM2Dmerge(par,data,i,k)
% merge particles
                data.contactsGlued(i,k) = true;
                data.contactsGlued(k,i) = true;

                data.contactsParticle.mergedParticles = true;
                % deactivate particles
                if(data.contactsParticle.deactivated(i) == 0 && data.contactsParticle.deactivated(k) == 0)
                    nMerged = data.contactsMerged.N + 1 ; data.contactsMerged.N = nMerged;
                    data.contactsMerged.aggregateSize(nMerged) = 2;
                    data.contactsMerged.timeFlag(nMerged) = true;
                    data.contactsMerged.index(nMerged,1:2) = [i ;k];
                    data.contactsMerged.position(:,nMerged) = [data.position(:,i); data.position(:,k)];
                    data.contactsMerged.mass(nMerged) = (data.mass(i)+data.mass(k)); 
                elseif(data.contactsParticle.deactivated(i) == 1 && data.contactsParticle.deactivated(k) == 1)
                    return;
                else
                    nMerged = data.contactsMerged.N;% + 1 ; data.contactsMerged.N = nMerged;
                    data.contactsMerged.aggregateSize(nMerged) = data.contactsMerged.aggregateSize(nMerged) +1;
                    if(data.contactsParticle.deactivated(i))
                        data.contactsMerged.index(nMerged,data.contactsMerged.aggregateSize(nMerged)+1) = k;
                        data.contactsMerged.mass(nMerged) = data.contactsMerged.mass(nMerged) + data.mass(k);
                    elseif(data.contactsParticle.deactivated(k))
                        data.contactsMerged.index(nMerged,data.contactsMerged.aggregateSize(nMerged)+1) = i;
                        data.contactsMerged.mass(nMerged) = data.contactsMerged.mass(nMerged) + data.mass(i);
                    else
                        warning('Particle merge not successfull.')
                    end
                end
                data.contactsParticle.deactivated(i) = 1;
                data.contactsParticle.deactivated(k) = 1;
                % append merged particles                
                mergedMass = data.contactsMerged.mass(nMerged);
                tmp = nonzeros(data.contactsMerged.index(nMerged,:));
                for jj = 1:data.contactsMerged.aggregateSize(nMerged) 
                    k = tmp(jj); 
                    data.contactsMerged.positionMerged(:,nMerged) = data.contactsMerged.positionMerged(:,nMerged) + data.position(:,k)*data.mass(k);

                    data.contactsMerged.velocityMerged(:,nMerged) =data.contactsMerged.velocityMerged(:,nMerged)+data.velocity(:,k)*data.mass(k);
                    %data.contactsMerged.relativePosition(:,nMerged) = [data.position(:,i)-data.contactsMerged.positionMerged(:,nMerged); data.position(:,k)-data.contactsMerged.positionMerged(:,nMerged)];
                    
                    %data.contactsMerged.relativePosition(2*jj-1:2*jj,nMerged) = data.position(:,k)-data.contactsMerged.positionMerged(:,nMerged);
                end                
                data.contactsMerged.positionMerged(:,nMerged) = 1/mergedMass*(data.contactsMerged.positionMerged(:,nMerged));
                data.contactsMerged.velocityMerged(:,nMerged) =1/mergedMass *data.contactsMerged.velocityMerged(:,nMerged);              
                data.contactsMerged.inertiaTensor(nMerged) = 0;
                for jj = 1:data.contactsMerged.aggregateSize(nMerged) 
                    k = tmp(jj);              
                    data.contactsMerged.relativePosition(2*jj-1:2*jj,nMerged) = data.position(:,k)-data.contactsMerged.positionMerged(:,nMerged);
                    data.contactsMerged.inertiaTensor(nMerged) = data.contactsMerged.inertiaTensor(nMerged) + 0.5*data.radius(k)^2*data.mass(k)+ data.mass(k)*norm(data.contactsMerged.relativePosition(2*jj-1:2*jj,nMerged));
                end 
end