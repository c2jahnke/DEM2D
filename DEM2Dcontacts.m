classdef DEM2Dcontacts < handle
    properties
    contacts = [];
    numParticleContacts = 0;
    numWallContacts = 0; 
    numToolContacts = 0; 
    numConstraints = 0;
    end
    methods
        function c =DEM2Dcontacts(data,par)
            X = data.position;

             if strcmp(par.software,'MATLAB')
                DT = delaunayTriangulation(X');
                E = edges(DT);
                
              elseif  strcmp(par.software,'GNU Octave') 
                DT = delaunay(X(1,:),X(2,:)); %
                if length(DT) == 1
                  E = [];
                else
                  E = [DT(:,1:2); DT(:,2:3); DT(:,[3,1])];
                end
             end
             maxr=max([data.radius]);
                for l = 1:length(E')
                    i=E(l,2);
                    j=E(l,1);
                    Di=data.position(1:2,j) - data.position(1:2,i);
                    Di=Di';
                    nD=norm(Di);
                    d=nD-data.radius(j)-data.radius(i);
                    if d < (par.collisionThreshold-1)*maxr
                        normal=Di'/nD;
                        c.numConstraints=c.numConstraints+1;
                        c.numParticleContacts = c.numParticleContacts +1;
                        contacts(c.numConstraints) = DEM2Dcontact(i,j,data,par,normal,-d);
                    end

                end
            if(isempty(E) && par.N > 1) % points collinear?
                d = DEM2Ddist(X(1,:),X(2,:));
                for l = 1:length(X)
                    for k = l+1:length(X)
                        d(l,k) = d(l,k) - data.radius(l) - data.radius(k);
                        if(d(l,k)<maxr*(par.collisionThreshold-1))
                            Di=data.position(1:2,l) - data.position(1:2,k);
                            Di=Di';
                            nD=norm(Di);
                            normal = Di'/nD;
                            c.numConstraints=c.numConstraints+1;
                            c.numParticleContacts = c.numParticleContacts +1;
                            contacts(c.numConstraints) = DEM2Dcontact(k,l,data,par,normal,-d(l,k));
                        end
                    end
                end
            end
            
            % search for wall-contacts
            for j = 1:par.N
                x(j) = data.position(1,j); z(j) = data.position(2,j);
                if (x(j)- par.bBox(1) < data.radius(j)*par.collisionThreshold)%((x(i)< par.bBox(1)+data.radius(i))) % left
                    normal = [1;0]; d = data.radius(j) -abs(par.bBox(1)- x(j));
                    c.numConstraints=c.numConstraints+1;
                    c.numWallContacts = c.numWallContacts+1;
                    contacts(c.numConstraints) = DEM2Dcontact(-1,j,data,par,normal,d);
                end
                if (par.bBox(2)-x(j)<data.radius(j)*par.collisionThreshold)%((x(i)> par.bBox(2)-data.radius(i))) % right
                    normal = [-1;0]; d = data.radius(j) -abs(par.bBox(2) - x(j));
                    c.numConstraints=c.numConstraints+1;
                    c.numWallContacts = c.numWallContacts+1;
                    contacts(c.numConstraints) = DEM2Dcontact(-2,j,data,par,normal,d);
                end
                if z(j) - par.bBox(3)  < data.radius(j)*par.collisionThreshold%(z(i)< (par.bBox(3)+data.radius(i)))%z(i) - par.bBox(3)  < data.radius(i)*par.collisionThreshold% % bottom
                    normal = [0;1]; d = data.radius(j) - abs(par.bBox(3) - z(j));
                    c.numConstraints=c.numConstraints+1;
                    c.numWallContacts = c.numWallContacts+1;
                    contacts(c.numConstraints) = DEM2Dcontact(-3,j,data,par,normal,d);
                end
                if (par.bBox(4) - z(j) < data.radius(j)*par.collisionThreshold) %((z(i)> par.bBox(4)-data.radius(i))) % top
                    normal = [0;-1]; d = data.radius(j) -abs(par.bBox(4) - z(j));
                    c.numConstraints=c.numConstraints+1;
                    c.numWallContacts = c.numWallContacts+1;
                    contacts(c.numConstraints) = DEM2Dcontact(-4,j,data,par,normal,d);
                end
            end
            if(par.toolBool)
            % search for tool-contacts
            for j = 1:par.N
                x(j) = data.position(1,j); z(j) = data.position(2,j);
                if (x(j)- data.toolbBox(1) < data.radius(j)*par.collisionThreshold && (z(j) < data.toolbBox(4) && z(j) > data.toolbBox(3)))%abs(z(j)-data.toolbBox(4)) <  data.radius(j)|| abs(z(j)-data.toolbBox(3)) < data.radius(j)) % left
                    normal = [1;0]; d = data.radius(j) -abs(data.toolbBox(1)- x(j));
                    c.numConstraints=c.numConstraints+1;
                    c.numToolContacts = c.numToolContacts+1;
                    contacts(c.numConstraints) = DEM2Dcontact(-5,j,data,par,normal,d);
                end
                if (data.toolbBox(2)-x(j)<data.radius(j)*par.collisionThreshold && (z(j) < data.toolbBox(4) && z(j) > data.toolbBox(3)))%&& (abs(z(j)-data.radius(j)) < data.toolbBox(3) || abs(z(j)-data.radius(j)) > data.toolbBox(4))) % right
                    normal = [-1;0]; d = data.radius(j) -abs(data.toolbBox(2) - x(j));
                    c.numConstraints=c.numConstraints+1;
                    c.numToolContacts = c.numToolContacts+1;
                    contacts(c.numConstraints) = DEM2Dcontact(-6,j,data,par,normal,d);
                end
                if (z(j) - data.toolbBox(3)  < data.radius(j)*par.collisionThreshold && (x(j) < data.toolbBox(2) && x(j) > data.toolbBox(1)))%&& (abs(x(j)-data.radius(j)) < data.toolbBox(2) || abs(x(j)-data.radius(j)) > data.toolbBox(1)))% bottom
                    normal = [0;1]; d = data.radius(j) - abs(data.toolbBox(3) - z(j));
                    c.numConstraints=c.numConstraints+1;
                    c.numToolContacts = c.numToolContacts+1;
                    contacts(c.numConstraints) = DEM2Dcontact(-7,j,data,par,normal,d);
                end
                if (data.toolbBox(4) - z(j) < data.radius(j)*par.collisionThreshold && (x(j) < data.toolbBox(2) && x(j) > data.toolbBox(1)))%&& (abs(x(j)-data.radius(j)) < data.toolbBox(1) || abs(x(j)-data.radius(j)) > data.toolbBox(2))) % top
                    normal = [0;-1]; d = data.radius(j) -abs(data.toolbBox(4) - z(j));
                    c.numConstraints=c.numConstraints+1;
                    c.numToolContacts = c.numToolContacts+1;
                    contacts(c.numConstraints) = DEM2Dcontact(-8,j,data,par,normal,d);
                end % upper left corner
                if (norm(diag(data.toolbBox) - data.position(:,j)) < data.radius(j)*par.collisionThreshold && (z(j) > data.toolbBox(4) && x(j) < data.toolbBox(1)))%&& (abs(x(j)-data.radius(j)) < data.toolbBox(1) || abs(x(j)-data.radius(j)) > data.toolbBox(2))) % top
                    normal = (diag(data.toolbBox) - data.position(:,j))/norm(diag(data.toolbBox) - data.position(:,j)); 
                    d = data.radius(j) - norm(diag(data.toolbBox) - data.position(:,j));
                    c.numConstraints=c.numConstraints+1;
                    c.numToolContacts = c.numToolContacts+1;
                    contacts(c.numConstraints) = DEM2Dcontact(-9,j,data,par,normal,d);
                end % lower right corner
                if (norm([data.toolbBox(1); data.toolbBox(3)] - data.position(:,j))<data.radius(j)*par.collisionThreshold && (x(j) > data.toolbBox(1) && z(j) < data.toolbBox(3)))%&& (abs(z(j)-data.radius(j)) < data.toolbBox(3) || abs(z(j)-data.radius(j)) > data.toolbBox(4))) % right
                    normal = ([data.toolbBox(2); data.toolbBox(3)] - data.position(:,j))/norm([data.toolbBox(2);data.toolbBox(3)] - data.position(:,j)); 
                    d = data.radius(j) -norm([data.toolbBox(2); data.toolbBox(3)] - data.position(:,j));
                    c.numConstraints=c.numConstraints+1;
                    c.numToolContacts = c.numToolContacts+1;
                    contacts(c.numConstraints) = DEM2Dcontact(-10,j,data,par,normal,d);
                end % lower left corner
                if (norm(data.toolbBox(1,:)' - data.position(:,j))  < data.radius(j)*par.collisionThreshold && (x(j) < data.toolbBox(1) && z(j) < data.toolbBox(3)))%&& (abs(x(j)-data.radius(j)) < data.toolbBox(2) || abs(x(j)-data.radius(j)) > data.toolbBox(1)))% bottom
                    normal = (data.toolbBox(1,:)' - data.position(:,j))/norm(data.toolbBox(1,:)' - data.position(:,j)); 
                    d = data.radius(j) - norm(data.toolbBox(1,:)' - data.position(:,j));
                    c.numConstraints=c.numConstraints+1;
                    c.numToolContacts = c.numToolContacts+1;
                    contacts(c.numConstraints) = DEM2Dcontact(-11,j,data,par,normal,d);
                end % upper right corner
                if (norm(data.toolbBox(2,:)' - data.position(:,j)) < data.radius(j)*par.collisionThreshold && (x(j) > data.toolbBox(2) && z(j) > data.toolbBox(4)))%&& (abs(x(j)-data.radius(j)) < data.toolbBox(1) || abs(x(j)-data.radius(j)) > data.toolbBox(2))) % top
                    normal = (data.toolbBox(2,:)' - data.position(:,j))/norm(data.toolbBox(2,:)' - data.position(:,j)); 
                    d = data.radius(j) - norm(data.toolbBox(2,:)' - data.position(:,j));
                    c.numConstraints=c.numConstraints+1;
                    c.numToolContacts = c.numToolContacts+1;
                    contacts(c.numConstraints) = DEM2Dcontact(-12,j,data,par,normal,d);
                end
                
            end
            end
            if ( c.numConstraints ~= c.numParticleContacts+c.numWallContacts+c.numToolContacts)
                warning;
            end
            if c.numConstraints==0
                contacts=[];
            end
            c.contacts = contacts;
        end
    end
end