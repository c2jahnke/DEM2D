function [data,par] =DEM2Dload(par,dataStr)
    temp = load(dataStr); % 'data.mat' 'data_50p.mat'
    data = temp.data;
    par.N = length(data.delta);
end
