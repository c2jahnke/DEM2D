function DEM2Dplot(data,par)
    drawArrow = @(x,y) quiver(x(1),y(1),x(2)-x(1),y(2)-y(1),0 );
    
    figure('units','normalized','outerposition',[0.3 0.0 0.5 1.0])
    h1 = gcf;
    
    set(gcf,'Color',[1 1 1]);
    hold on;
    p = data.position;
    v = data.velocity;     
    vnorm = norm(v);
    r = data.radius;              

    AngleRes = 45; 
    % Draw the particles
    theta = linspace(0,2*pi,AngleRes);
    scal = linspace(-1,1,AngleRes);
    s = sin(theta);
    c = cos(theta);


    % Partical position:
    for i=1:par.N
        x1 = [p(1,i) p(1,i)+v(1,i)/vnorm*r(i)];
        y1 = [p(2,i) p(2,i)+v(2,i)/vnorm*r(i)];
        plot([c*r(i)+p(1,i)*ones(1,AngleRes)],[s*r(i)+p(2,i)*ones(1,AngleRes)],'k-')
        hold on
        xCross = DEM2Drotation(data.angular(1,i))*[1;0];
        zCross = DEM2Drotation(data.angular(1,i))*[0;1];
        plot([xCross(1)*scal*r(i) + p(1,i)*ones(1,AngleRes)],[xCross(2)*scal*r(i) + p(2,i)*ones(1,AngleRes)],'k-')
        plot([zCross(1)*scal*r(i) + p(1,i)*ones(1,AngleRes)],[zCross(2)*scal*r(i) + p(2,i)*ones(1,AngleRes)],'k-')
        drawArrow(x1,y1); 
    end

    % Walls: 
    plot([par.bBox(1) par.bBox(1) par.bBox(1) par.bBox(2) par.bBox(2) par.bBox(1)],...
         [par.bBox(3) par.bBox(4) par.bBox(3) par.bBox(3) par.bBox(4) par.bBox(4)],'b-','LineWidth',2)
    if(par.toolBool)
        plot([par.toolbBox(1) par.toolbBox(1) par.toolbBox(1) par.toolbBox(2) par.toolbBox(2) par.toolbBox(1)],...
         [par.toolbBox(3) par.toolbBox(4) par.toolbBox(3) par.toolbBox(3) par.toolbBox(4) par.toolbBox(4)],'r-','LineWidth',2)
     end
    set(gca,'Visible','off','XTick',[],'YTick',[]);
    %axis([-0.95*par.r(2)+par.bBox(1) par.bBox(2)+0.95*par.r(2) -0.95*par.r(2)+par.bBox(3) par.bBox(4)+0.95*par.r(2)])
    axis equal
    hold off
end



