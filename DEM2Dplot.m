function DEM2Dplot(data,par)
p = data.position;
v = data.velocity;      
r = data.radius;
m = data.mass;               

AngleRes = 45; 
% Draw the particles
theta = linspace(0,2*pi,AngleRes);
scal = linspace(-1,1,AngleRes);
s = sin(theta);
c = cos(theta);


% Partical position:
for i=1:par.N
    plot([c*r(i)+p(1,i)*ones(1,AngleRes)],[s*r(i)+p(2,i)*ones(1,AngleRes)],'k-')
    hold on
    xCross = DEM2Drotation(data.angular(1,i))*[1;0];
    zCross = DEM2Drotation(data.angular(1,i))*[0;1];
    plot([xCross(1)*scal*r(i) + p(1,i)*ones(1,AngleRes)],[xCross(2)*scal*r(i) + p(2,i)*ones(1,AngleRes)],'k-')
    plot([zCross(1)*scal*r(i) + p(1,i)*ones(1,AngleRes)],[zCross(2)*scal*r(i) + p(2,i)*ones(1,AngleRes)],'k-')
end

% Walls: 
plot([par.bBox(1) par.bBox(1) par.bBox(1) par.bBox(2) par.bBox(2) par.bBox(1)],...
     [par.bBox(3) par.bBox(4) par.bBox(3) par.bBox(3) par.bBox(4) par.bBox(4)],'b-','LineWidth',2)


axis([-0.03+par.bBox(1) par.bBox(2)+0.03 -0.03+par.bBox(3) par.bBox(4)+0.03])
axis equal
hold off

end



