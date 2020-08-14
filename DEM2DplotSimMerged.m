function DEM2DplotSimMerged(P1,V1,A1,PM,VM,par,data,j)
    drawArrow = @(x,y) quiver(x(1),y(1),x(2)-x(1),y(2)-y(1),0,'b-');
    figure('units','normalized','outerposition',[0.3 0.4 0.5 0.6])
    
    if(par.writeVid)
        video = VideoWriter(['videos/' par.videoname '.avi']);
        video.FrameRate = par.video_framerate;
        open(video);

    end
    set(gcf,'Color',[1 1 1]);
    h1 = gcf; %subplot(1,2,1);

    r = data.radius;
    AngleRes = 45; 
    % Draw the particles
    theta = linspace(0,2*pi,AngleRes);
    scal = linspace(-1,1,AngleRes);
    s = sin(theta);
    c = cos(theta);
    DELAY = par.dt*par.VisualizationStep/par.N;
    
    for k=1:j    
         % Walls: 
        plot([par.bBox(1) par.bBox(1) par.bBox(1) par.bBox(2) par.bBox(2) par.bBox(1)],...
         [par.bBox(3) par.bBox(4) par.bBox(3) par.bBox(3) par.bBox(4) par.bBox(4)],'b-','LineWidth',2)
        %plot(reshape(P1(k,1,:),[1 par.N]),reshape(P1(k,2,:),[1 par.N]),'.','MarkerEdgeColor','w')
        hold all
        if(par.toolBool)
        plot([par.toolbBox(1) par.toolbBox(1) par.toolbBox(1) par.toolbBox(2) par.toolbBox(2) par.toolbBox(1)],...
         [par.toolbBox(3) par.toolbBox(4) par.toolbBox(3) par.toolbBox(3) par.toolbBox(4) par.toolbBox(4)],'r-','LineWidth',2)
        end
        % Partical position:
        for i=1:par.N
            if(data.contactsParticle.deactivated(i))
                for kk = 1:data.contactsMerged.N % inefficient
                    if(data.contactsMerged.index(kk,1) == i)
                        ii = nonzeros(data.contactsMerged.index(kk,:));
                        for j = 1:data.contactsMerged.aggregateSize(kk) 
                        %ii = data.contactsMerged.index(kk,2);
                            if(any(data.contactsMerged.timePoint(i,ii(j))) && data.contactsMerged.timePoint(i,ii(j)) <= k && i < ii(j))
                                text(P1(k,1,i),P1(k,2,i),[num2str(i) "M" num2str(kk)],'Color','red','fontsize',par.videoPartFontsize)
                                text(P1(k,1,ii(j)),P1(k,2,ii(j)),[num2str(ii(j)) "M" num2str(kk)],'Color','red','fontsize',par.videoPartFontsize)
                            elseif(data.contactsMerged.timePoint(i,ii(j)) > k)
                                text(P1(k,1,i),P1(k,2,i),num2str(i),'fontsize',par.videoPartFontsize)
                                text(P1(k,1,ii(j)),P1(k,2,ii(j)),num2str(ii(j)),'fontsize',par.videoPartFontsize)
                            end
                        end
                    end % this part is not fully implemented
                end
            else
                text(P1(k,1,i),P1(k,2,i),num2str(i),'fontsize',par.videoPartFontsize)
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

       

        title([par.videoname '_merged'],'fontsize',par.videoFontsize)
        text(par.bBox(2)-2.95*par.r(2),par.bBox(4)-0.95*par.r(2),['t = ',num2str((k-1)*par.dt*par.VisualizationStep,'%10.2f') 's'],'fontsize',par.videoFontsize)
        
        axis([-0.15+par.bBox(1) par.bBox(2)+0.15 -0.15+par.bBox(3) par.bBox(4)+0.15])
        axis equal
        
        %set(gca,'Visible','off','XTick',[],'YTick',[]);
        set(gca,'fontsize', par.videoFontsize);
        title(par.videoname,'fontsize',par.videoFontsize)
        drawnow;
        F(j) = getframe(h1);
        pause(DELAY)

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