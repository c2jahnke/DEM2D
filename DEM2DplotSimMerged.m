function DEM2DplotSimMerged(P1,V1,A1,PM,VM,par,data,j)
    drawArrow = @(x,y) quiver(x(1),y(1),x(2)-x(1),y(2)-y(1),0,'b-');
    figure('units','normalized','outerposition',[0.3 0.0 0.5 1.0])
    if(par.writeVid)
        video = VideoWriter(['videos/' par.videoname '.avi']);
        video.FrameRate = par.video_framerate;
        open(video);

    end
    h1 = gcf; %subplot(1,2,1);

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
                for kk = 1:data.contactsMerged.N % inefficient
                    if(data.contactsMerged.index(1,kk) == i)
                        ii = data.contactsMerged.index(2,kk);
                        if(data.contactsMerged.timePoint(i,ii) <= k)
                            text(P1(k,1,i),P1(k,2,i),[num2str(i) "M"],'Color','red','fontsize',par.videoFontsize)
                            text(P1(k,1,ii),P1(k,2,ii),[num2str(ii) "M"],'Color','red','fontsize',par.videoFontsize)
                        else
                            text(P1(k,1,i),P1(k,2,i),num2str(i),'fontsize',par.videoFontsize)
                            text(P1(k,1,ii),P1(k,2,ii),num2str(ii),'fontsize',par.videoFontsize)
                        end
                    end % this part is not fully implemented
                end
            else
                text(P1(k,1,i),P1(k,2,i),num2str(i),'fontsize',par.videoFontsize)
            end
            vnorm = norm(V1(k,:,i));
            x1 = [P1(k,1,i) P1(k,1,i)+V1(k,1,i)/vnorm*r(i)];
            y1 = [P1(k,2,i) P1(k,2,i)+V1(k,2,i)/vnorm*r(i)];
            if(vnorm >0.1)
                    drawArrow(x1,y1); 
            end

            plot(c*r(i)+P1(k,1,i)*ones(1,AngleRes),s*r(i)+P1(k,2,i)*ones(1,AngleRes),'k-')
            
            xCross = DEM2Drotation(A1(k,1,i))*[1;0];
            zCross = DEM2Drotation(A1(k,1,i))*[0;1];
            plot(xCross(1)*scal*r(i) + P1(k,1,i)*ones(1,AngleRes),xCross(2)*scal*r(i) + P1(k,2,i)*ones(1,AngleRes),'k-')
            plot(zCross(1)*scal*r(i) + P1(k,1,i)*ones(1,AngleRes),zCross(2)*scal*r(i) + P1(k,2,i)*ones(1,AngleRes),'k-')

            
        end

        % Walls: 
        plot([par.bBox(1) par.bBox(1) par.bBox(1) par.bBox(2) par.bBox(2) par.bBox(1)],...
         [par.bBox(3) par.bBox(4) par.bBox(3) par.bBox(3) par.bBox(4) par.bBox(4)],'b-','LineWidth',2)

        title('Simulation','fontsize',par.videoFontsize)
        text(par.bBox(2)-0.55,par.bBox(4)-0.25,['t = ',num2str((k-1)*par.dt*par.VisualizationStep,'%10.2f') 's'],'fontsize',par.videoFontsize)
        axis([-0.15+par.bBox(1) par.bBox(2)+0.15 -0.15+par.bBox(3) par.bBox(4)+0.15])
        axis equal

        set(gca,'fontsize', par.videoFontsize);
        drawnow;
        F(j) = getframe(h1);
        %pause(DELAY)
        %size(F(j).cdata)
        %hf=figure('Position', [100, 100, 1049, 895]);

        %pause(0.1)
        if(par.writeVid)
             writeVideo(video,F(j));
        end
        if(par.writePng)
            saveas(h1,['videos/' par.videoname '_' num2str(k) '.png'])
        end
        if(par.writeEps)
            saveas(h1,['videos/' par.videoname '_' num2str(k) '.eps'])
        end
        if(par.writePdf)
            
%             saveas(h1,['videos/' par.videoname '_' num2str(k) '.pdf'])
            screenposition = get(h1,'Position');
            set(h1,...
            'PaperPosition',[0 0 screenposition(3:4)],...
            'PaperSize',[screenposition(3:4)]);
            outPutFile = ['videos/' par.videoname '_' num2str(k)];
            print -dpdf outPutFile
            saveas(h1,['videos/' par.videoname '_' num2str(k) '.pdf'])
        end
        hold off
    end
    if(par.writeVid)
        close(video)
    end
end