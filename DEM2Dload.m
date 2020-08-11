function [data,par] =DEM2Dload(par)
    temp = load('data.mat');
    data = temp.data;
    par.N = length(data.delta);
end
