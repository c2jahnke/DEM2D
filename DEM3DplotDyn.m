function DEM3DplotDyn(data,par,P1,k)
r = data.radius; 

[Xs,Ys,Zs] = sphere(14);
axis([-0.03+par.bBox(1) par.bBox(2)+0.03 0.3 0.7 -0.03+par.bBox(3) par.bBox(4)+0.03 ])
% Partical position:
for i=1:par.N 
    surf(Xs.*r(i)+P1(k,1,i)*ones(15,15) ,Ys.*r(i)+0.5*ones(15,15), Zs.*r(i)+P1(k,2,i)*ones(15,15))
    hold on
    %axis equal
end

end