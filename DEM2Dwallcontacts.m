classdef DEM2Dwallcontacts < handle
    properties
    contacts = [];
    end
    methods
        function c =DEM2Dwallcontacts(data,par)
            delta = sparse(0)
             x1 = par.bBox(1,:)+ data.radius(end); %lower bounding point
             x2 = par.bBox(2,:)- data.radius(end);%upper bounding point
              U = repmat(x1,par.N,1);
              V = repmat(x2,par.N,1);
              dist = data.position' -U; %distance from lower
              dist2 = data.position' -V; %distance from upper
              dist3 = dist - abs(dist);
              dist4 = dist2 + abs(dist2);
              indexLow = find(dist3);
              indexUp = find(dist4);
              numconstraints = length(indexLow) + length(indexUp);
              if numconstraints==0
                contacts=[];
                c.contacts = contacts;
                return;
              else
              for k = 1:length(indexLow)
                  iA = -1
                  deltaLow(indexLow) = dist(indexLow)
                  if(indexLow < par.N)
                     iA = -3
                  end
                  normal =  - data.position(indexLow)
                  contacts =  DEM2Dcontact(iA,iB,data,normal,d);
              end
              for l = 1:length(indexUp) 
                  deltaUp(indexLow) = dist2(indexLow)
              end
              end
            c.contacts = contacts;
        end
    end
end
