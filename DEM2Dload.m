function [data,par] =DEM2Dload(par)
    temp = load('data_3pMerged.mat');
    data = temp.data;
    par.N = length(data.delta);
end
