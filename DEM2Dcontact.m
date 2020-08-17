classdef DEM2Dcontact < handle
    properties
    a = 0;
    b = 0;
    n = zeros(2,1);
    t = zeros(2,1);
    
    % for PGJ
    Ga = zeros(2,2);
    Gb = zeros(2,2);
    l = zeros(2,1);
    rhs =zeros(2,1);
    distance = 0;
    eta = 0;
    compliance = 0;
    relVel = 0;
    redMass = 0;
    
    end
    methods
        function c = DEM2Dcontact(iA,iB,data,par,normal,d)
            c.a = iA;
            c.b = iB;
            c.n = -normal;
            c.t = [c.n(2) ;-c.n(1)];
            c.distance = d;
            % PGJ stuff
            
            c.Gb(1:2,1)=normal;
            c.Gb(1:2,2)=[-normal(2);normal(1)];%-c.t;
            c.compliance = 1/(par.dt^2*par.Emodul);
            if(iA>0)
                c.Ga(1:2,1)=-normal;
                c.Ga(1:2,2)=[normal(2);-normal(1)];%c.t;
                c.redMass = data.mass(iA)*data.mass(iB)/(data.mass(iA) + data.mass(iB));
                c.relVel= c.Gb(:,2)'*data.velocity(:,iB) + c.Ga(:,2)'*data.velocity(:,iA);%(data.velocity(:,iA) - data.velocity(:,iB))'*c.t;
            else
                c.redMass = data.mass(iB);
                c.relVel= c.Gb(:,2)'*data.velocity(:,iB);%data.velocity(:,iB)'*c.t;
            end
            c.eta = 2/( trace(c.Ga'*c.Ga + c.Gb'*c.Gb)+2*c.compliance)*c.redMass;
            c.rhs=[-d/par.dt + par.mu*abs(c.relVel);0];
        end
        function delta=ConeProj(c,par,delta) % from Jan Kleinerts code
            absDelta2=abs(delta(2));
            if absDelta2<=par.mu*delta(1)
                %sticking
            elseif absDelta2 <= - delta(1)/par.mu
                %breaking/no contact
                delta=[0;0];
            else
                %sliding
                temp=(delta(1)+par.mu*absDelta2)/(par.mu^2+1);
                delta=temp*[1;sign(delta(2))*par.mu];
            end
        end

    end
end