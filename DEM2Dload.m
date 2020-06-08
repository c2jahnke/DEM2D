function [data,par] =DEM2Dload(par)
    temp = load('data200.mat');
    data = temp.data;
    par.N = length(data.delta);
end
