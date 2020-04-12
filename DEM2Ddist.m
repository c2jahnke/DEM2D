function d = DEM2Ddist(x,y)

% Calculate the distance between particles

N = length(x);

Dx = zeros(N,N);
Dy = zeros(N,N);

for i=1:N
    for j=1:N
        Dx(i,j) = x(i)-x(j);
        Dy(i,j) = y(i)-y(j);
    end
end

d = sqrt(Dx.^2 + Dy.^2);

return