classdef DEM2Dcontact < handle
    properties
    a;
    b;
    n;
    t;
    end
    methods
        function c = DEM2Dcontact(iA,iB,data,normal,d)
            c.a = iA;
            c.b = iB;
            c.n = normal;
            c.t = [c.n(2) ;-c.n(1)];
        end

    end
end