function [data,par] =DEM2Dload(par)
    temp = load('showcases/AngleOfRepose200p/data.mat');
    data = temp.data;
    par.N = length(data.delta);
end
