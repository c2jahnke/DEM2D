function DEM2DplotDyn(data,par,P1,A1,k)
p = data.position;
v = data.velocity;      
r = data.radius;
m = data.mass;               

% Draw the particles
theta = linspace(0,2*pi,15);
scal = linspace(-1,1,15);
s = sin(theta);
c = cos(theta);


% Partical position:
for i=1:par.N
    plot([c*r(i)+P1(k,1,i)*ones(1,15)],[s*r(i)+P1(k,2,i)*ones(1,15)],'k-')
    hold all
    xCross = DEM2Drotation(A1(k,1,i))*[1;0];
    zCross = DEM2Drotation(A1(k,1,i))*[0;1];
    plot([xCross(1)*scal*r(i) + P1(k,1,i)*ones(1,15)],[xCross(2)*scal*r(i) + P1(k,2,i)*ones(1,15)],'k-')
    plot([zCross(1)*scal*r(i) + P1(k,1,i)*ones(1,15)],[zCross(2)*scal*r(i) + P1(k,2,i)*ones(1,15)],'k-')

end

% Walls: 
plot([par.bBox(1) par.bBox(1) par.bBox(1) par.bBox(2) par.bBox(2) par.bBox(1)],...
     [par.bBox(3) par.bBox(4) par.bBox(3) par.bBox(3) par.bBox(4) par.bBox(4)],'b-','LineWidth',2)


end