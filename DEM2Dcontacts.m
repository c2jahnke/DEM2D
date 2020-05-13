classdef DEM2Dcontacts < handle
    properties
    contacts = [];
    numParticleContacts = 0;
    numWallContacts = 0; 
    numConstraints = 0;
    end
    methods
        function c =DEM2Dcontacts(data,par)
            X = data.position;
            
            if strcmp(par.software,'MATLAB')
               
                DT = delaunayTriangulation(X');
                E = edges(DT);
      
            elseif  strcmp(par.software,'GNU Octave') 
                    
              DT = delaunay(X(1,:),X(2,:));
              if length(DT) == 1
                E = [1,2];
              else
                E = [DT(:,1:2); DT(:,2:3); DT(:,[3,1])];
              end
            end
                DT = delaunayTriangulation(X');
                E = edges(DT);
                maxr=max([data.radius]);
                for k = 1:size(E,1)
                    i=E(k,1);
                    j=E(k,2);
                    Di=data.position(1:2,i) - data.position(1:2,j);
                    Di=Di';
                    nD=norm(Di);
                    d=nD-data.radius(i)-data.radius(j);
                    if d<par.collisionThreshold*maxr
                        normal=Di'/nD;
                        c.numConstraints=c.numConstraints+1;
                        c.numParticleContacts = c.numParticleContacts +1;
                        contacts(c.numConstraints) = DEM2Dcontact(i,j,data,normal,d);
                    end
                end
                if(isempty(E) && par.N > 1) % points collinear?
                    d = DEM2Ddist(X(1,:),X(2,:));
                    for k = 1:length(X)
                        for l = k+1:length(X)
                            d(k,l) = d(k,l) - data.radius(k) - data.radius(l);
                            if(d(k,l)<0)
                                Di=data.position(1:2,k) - data.position(1:2,l);
                                Di=Di';
                                nD=norm(Di);
                                normal = Di'/nD;
                                c.numConstraints=c.numConstraints+1;
                                c.numParticleContacts = c.numParticleContacts +1;
                                contacts(c.numConstraints) = DEM2Dcontact(k,l,data,normal,d(k,l));
                            end
                        end
                    end
                end
            elseif  strcmp(par.software,'GNU Octave') 
                %disp('Warning, currently not working')
%                 DT = delaunay(X(1,:),X(2,:)); %fix, currently not working
%                 triplot(DT,X(1,:),X(2,:));
% %               TRI = DelaunayTri(X');
% %                   E=edges(TRI);

                d = DEM2Ddist(X(1,:),X(2,:));
                for k = 1:par.N-1
                    for l = max(1,min(k+1,par.N)):par.N
                       d(k,l) = d(k,l) - data.radius(k) - data.radius(l);
                            if(d(k,l)<0)
                                Di=data.position(1:2,k) - data.position(1:2,l);
                                Di=Di';
                                nD=norm(Di);
                                normal = Di'/nD;
                                c.numConstraints=c.numConstraints+1;
                                c.numParticleContacts = c.numParticleContacts +1;
                                contacts(c.numConstraints) = DEM2Dcontact(k,l,data,normal,d(k,l));
                            end
                    end
                end 
            end

            
            % search for wall-contacts
            for i = 1:par.N
                x(i) = data.position(1,i); z(i) = data.position(2,i);
                if ((x(i)< par.bBox(1)+data.radius(i))) % left
                    normal = [1;0];
                    c.numConstraints=c.numConstraints+1;
                    c.numWallContacts = c.numWallContacts+1;
                    contacts(c.numConstraints) = DEM2Dcontact(-1,i,data,normal,0);
                end
                if ((x(i)> par.bBox(2)-data.radius(i))) % right
                    normal = [-1;0];
                    c.numConstraints=c.numConstraints+1;
                    c.numWallContacts = c.numWallContacts+1;
                    contacts(c.numConstraints) = DEM2Dcontact(-2,i,data,normal,0);
                end
                if ((z(i)< par.bBox(3)+data.radius(i))) % bottom
                    normal = [0;1];
                    c.numConstraints=c.numConstraints+1;
                    c.numWallContacts = c.numWallContacts+1;
                    contacts(c.numConstraints) = DEM2Dcontact(-3,i,data,normal,0);
                end
                if ((z(i)> par.bBox(4)-data.radius(i))) % top
                    normal = [0;-1];
                    c.numConstraints=c.numConstraints+1;
                    c.numWallContacts = c.numWallContacts+1;
                    contacts(c.numConstraints) = DEM2Dcontact(-4,i,data,normal,0);
                end
            end
            if ( c.numConstraints ~= c.numParticleContacts+c.numWallContacts)
                warning;
            end
            if c.numConstraints==0
                contacts=[];
            end
            c.contacts = contacts;
        end
    end
end