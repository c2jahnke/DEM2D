classdef DEM2Dcontacts < handle
    properties
    contacts = [];
    end
    methods
        function c =DEM2Dcontacts(data,par)
            X = data.position;
            %DT = delaunayTriangulation(X);
            Tri = DelaunayTri(X');
            E = edges(Tri);
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
            if numconstraints==0
                contacts=[];
            end;
            c.contacts = contacts
        end
           
        
    end
end