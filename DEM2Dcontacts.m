classdef DEM2Dcontacts < handle
    properties
    contacts = [];
    end
    methods
        function c =DEM2Dcontacts(data,par)
            X = data.position;
              disp('Warning')
             if strcmp(par.software,'MATLAB')
              DT = delaunayTriangulation(X');
              E = edges(DT);
            elseif  strcmp(par.software,'GNU Octave') 
                % DT = DelaunayTri(X',X'); fix, currently not working
                %E=edges(DT);
            end
            numconstraints = 0;
            maxr=max([data.radius]);
            for k = 1:length(E')
                i=E(k,1);
                j=E(k,2);
                Di=data.position(1:2,i) - data.position(1:2,j);
                Di=Di';
                nD=norm(Di);
                d=nD-data.radius(i)-data.radius(j);
                if d<par.collisionThreshold*maxr
                    normal=Di'/nD;
                    numconstraints=numconstraints+1;
                    contacts(numconstraints) = DEM2Dcontact(i,j,data,normal,d);
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
                            numconstraints=numconstraints+1;
                            contacts(numconstraints) = DEM2Dcontact(k,l,data,normal,d(k,l));
                        end
                    end
                end
            end
            if numconstraints==0
                contacts=[];
            end
            c.contacts = contacts;
        end
    end
end