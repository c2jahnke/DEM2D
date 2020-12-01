function DEM2DplotAngularVelocity(data,par,output,visuIndex)

for i = 1:par.N
plot(output.timeInc,output.angular(output.timeInc,2,i)*180/pi)
hold on
title(['Angular Velocity of particle ' num2str(i)])
end
end
