function DEM2DplotSimMerged(P1,V1,A1,PM,VM,par,data,j)
    drawArrow = @(x,y) quiver(x(1),y(1),x(2)-x(1),y(2)-y(1),0,'b-');
    figure('units','normalized','outerposition',[0 0 1 1])
    if(par.writeVid)
        video = VideoWriter(['videos/' par.videoname '.avi']);
        video.FrameRate = par.video_framerate;
        open(video);

    end
    h1 = subplot(1,2,1);

    r = data.radius;
    AngleRes = 45; 
    % Draw the particles
    theta = linspace(0,2*pi,AngleRes);
    scal = linspace(-1,1,AngleRes);
    s = sin(theta);
    c = cos(theta);
    DELAY = par.dt*par.VisualizationStep*0.4/par.N;
    for k=1:j    
        plot(reshape(P1(k,1,:),[1 par.N]),reshape(P1(k,2,:),[1 par.N]),'.','MarkerEdgeColor','w')
        hold all
        % Partical position:
        for i=1:par.N
            if(data.contactsParticle.deactivated(i))
                for kk = 1:data.contactsMerged.N
                    if(data.contactsMerged.index(1,kk) == i)
                        P1(k,:,i) = PM(k,:,kk); % wrong: check index
                    end % this part is not fully implemented
                end
            else
            text(P1(k,1,i),P1(k,2,i),num2str(i))
            vnorm = norm(V1(k,:,i));
            x1 = [P1(k,1,i) P1(k,1,i)+V1(k,1,i)/vnorm*r(i)];
            y1 = [P1(k,2,i) P1(k,2,i)+V1(k,2,i)/vnorm*r(i)];
            if(vnorm >0.01)
                    drawArrow(x1,y1); 
            end

            plot(c*r(i)+P1(k,1,i)*ones(1,AngleRes),s*r(i)+P1(k,2,i)*ones(1,AngleRes),'k-')
            
            xCross = DEM2Drotation(A1(k,1,i))*[1;0];
            zCross = DEM2Drotation(A1(k,1,i))*[0;1];
            plot(xCross(1)*scal*r(i) + P1(k,1,i)*ones(1,AngleRes),xCross(2)*scal*r(i) + P1(k,2,i)*ones(1,AngleRes),'k-')
            plot(zCross(1)*scal*r(i) + P1(k,1,i)*ones(1,AngleRes),zCross(2)*scal*r(i) + P1(k,2,i)*ones(1,AngleRes),'k-')

            end
        end

        % Walls: 
        plot([par.bBox(1) par.bBox(1) par.bBox(1) par.bBox(2) par.bBox(2) par.bBox(1)],...
         [par.bBox(3) par.bBox(4) par.bBox(3) par.bBox(3) par.bBox(4) par.bBox(4)],'b-','LineWidth',2)

        title('Simulation')
        text(par.bBox(2)-0.35,par.bBox(4)-0.25,['t = ',num2str((k-1)*par.dt*par.VisualizationStep,'%10.2f') 's'])
        axis([-0.15+par.bBox(1) par.bBox(2)+0.15 -0.15+par.bBox(3) par.bBox(4)+0.15])
        axis equal

        F(j) = getframe(h1);

        %pause(DELAY)
        %size(F(j).cdata)

        if(par.writeVid)
             writeVideo(video,F(j));
        end
        if(par.writePng)
            saveas(h1,['videos/' par.videoname '_' num2str(k) '.png'])
        end
        hold off
    end
    if(par.writeVid)
        close(video)
    end
end