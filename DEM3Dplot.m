function DEM3Dplot(data,par)
p = data.position;
r = data.radius; 

[Xs,Ys,Zs] = sphere(14);
figure (1);
% Partical position:
for i=1:par.N 
    surf(Xs.*r(i)+p(1,i)*ones(15,15) ,Ys.*r(i)+0.5*ones(15,15), Zs.*r(i)+p(2,i)*ones(15,15))
    hold on
end
plot([par.bBox(1) par.bBox(1) par.bBox(1) par.bBox(2) par.bBox(2) par.bBox(1)],...
     [par.bBox(3) par.bBox(4) par.bBox(3) par.bBox(3) par.bBox(4) par.bBox(4)],'b-','LineWidth',2)

axis([-0.03+par.bBox(1) par.bBox(2)+0.03 0.3 0.7 -0.03+par.bBox(3) par.bBox(4)+0.03 ])
end