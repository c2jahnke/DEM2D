function DEM2DplotSim(P1,A1,par,data,j)

for k=1:j    
r = data.radius;

% Draw the particles
theta = linspace(0,2*pi,15);
scal = linspace(-1,1,15);
s = sin(theta);
c = cos(theta);
DELAY = par.dt*par.step*0.04;

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
    hold    on
    title('Simulation')
    text(par.bBox(1)-0.15,par.bBox(3)-0.08,['t = ',num2str((k-1)*par.dt*par.step,'%10.2f') 's'])
    axis([-0.03+par.bBox(1) par.bBox(2)+0.03 -0.03+par.bBox(3) par.bBox(4)+0.03])
    axis equal
    hold off 
    F(j) = getframe();
    pause(DELAY)
end
end