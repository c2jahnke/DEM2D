function DEM2Dplot(data,par)
p = data.position;
v = data.velocity;      
r = data.radius;
m = data.mass;               

% Draw the particles
theta = linspace(0,2*pi,15);
s = sin(theta);
c = cos(theta);


% Partical position:
for i=1:par.N
    plot([c*r(i)+p(1,i)*ones(1,15)],[s*r(i)+p(2,i)*ones(1,15)],'r-')
    hold on
end

% Walls: 
plot([par.bBox(1) par.bBox(1) par.bBox(1) par.bBox(2) par.bBox(2) par.bBox(1)],...
     [par.bBox(3) par.bBox(4) par.bBox(3) par.bBox(3) par.bBox(4) par.bBox(4)],'b-','LineWidth',2)


axis([-0.03+par.bBox(1) par.bBox(2)+0.03 -0.03+par.bBox(3) par.bBox(4)+0.03])
axis equal
hold off

end



