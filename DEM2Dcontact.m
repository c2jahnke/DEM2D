classdef DEM2Dcontact < handle
    properties
    a;
    b;
    n;
    t;
%     Ga;
%     Gb;
%     l;
%     rhs;
%     eta;
%     compliance
    end
    methods
        function c = DEM2Dcontact(iA,iB,data,normal,d)
            c.a = iA;
            c.b = iB;
            c.n = normal;
            c.t = [c.n(2) ;-c.n(1)];
%             c.Ga = [n t];
%             c.Gb = - c.Ga;
%             c.l=[0;0];
%             c.rhs=[0;0];
%             c.compliance = 0;
%             c.eta = 2/( trace(c.Ga'*c.Ga + c.Gb'*c.Gb)+2*c.compliance);
        end

    end
end