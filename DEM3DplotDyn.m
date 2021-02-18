function DEM3DplotDyn(P1,A1,data,par,j)
r = data.radius; 

[Xs,Ys,Zs] = sphere(14);
axis([-0.03+par.bBox(1) par.bBox(2)+0.03 0.3 0.7 -0.03+par.bBox(3) par.bBox(4)+0.03 ])
DELAY = par.dt*par.VisualizationStep*0.1/par.N*10;
CO = zeros(15,15,3);
CO(:,:,1) = zeros(15); % red
CO(:,:,2) = ones(15).*linspace(0.5,0.6,15); % green
CO(:,:,3) = ones(15).*linspace(0,1,15); % blue
% Partical position:
    for k=1:j 
        h1 = gcf;
        for i=1:par.N 
            
            xCross = DEM2Drotation(A1(k,1,i))*[1;0];
            zCross = DEM2Drotation(A1(k,1,i))*[0;1];
            s = surf(Xs.*r(i)+P1(k,1,i)*ones(15,15) ,Ys.*r(i)+0.5*ones(15,15), Zs.*r(i)+P1(k,2,i)*ones(15,15),CO);
            %rotate(s,[0 1 0],A1(k,1,i)*180/pi)
            hold on
            axis equal
        end
        title('Simulation','fontsize',par.videoFontsize)
        text(par.bBox(2)-0.55,par.bBox(4)-0.25,['t = ',num2str((k-1)*par.dt*par.VisualizationStep,'%10.2f') 's'],'fontsize',par.videoFontsize)
        axis([-0.15+par.bBox(1) par.bBox(2)+0.15 -0.15+par.bBox(3) par.bBox(4)+0.15])
        axis equal

        set(gca,'fontsize', par.videoFontsize);
        drawnow;
        F(j) = getframe(h1);
        pause(DELAY)
        hold off
    end
end