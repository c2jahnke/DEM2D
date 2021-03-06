function DEM2DplotEnergy(data,par,output)
%% Plot energy in the particulate system
%% total energy should decrease due to friction and damping
%% only for par.gVert = 0;
V = zeros(par.N,output.finalVisuIndex); % speed
A = zeros(par.N,output.finalVisuIndex); % angular velocity
P = zeros(par.N,output.finalVisuIndex); % height
Ekin = zeros(1,output.finalVisuIndex);
Erot = zeros(1,output.finalVisuIndex);
Epot = zeros(1,output.finalVisuIndex);
Etotal = zeros(1,output.finalVisuIndex);
for k = 1 : output.finalVisuIndex
    for j = 1: par.N
       V(j,k) = norm( output.velocity(k,:,j));
       Ekin(k) = Ekin(k) + 1/2*data.mass(j)*V(j,k)^2;
       A(j,k) = norm(output.angular(k,2,j));
       Erot(k) = Erot(k) + 1/2*data.mass(j)*A(j,k)^2;
       P(j,k) = output.position(k,2,j) - par.bBox(3)-data.radius(j);
       Epot(k) = Epot(k) + data.mass(j)*(-par.g)*P(j,k);
    end
    Etotal(k) = Ekin(k) + Erot(k) + abs(Epot(k));
end
% kinetic energy
figure('Name','Total Energy')
subplot(2,1,1);
semilogy(output.timeInc*par.VisualResolution,Ekin(:))
hold on
semilogy(output.timeInc*par.VisualResolution,Erot(:))
hold on
semilogy(output.timeInc*par.VisualResolution,abs(Epot))
hold on
semilogy(output.timeInc*par.VisualResolution,Etotal)
ylabel("Energy [J = kg m/s^2]")
xlabel("Time [s]")
legend('Ekin', 'Erot','Epot','Etotal')
hold off
subplot(2,1,2)
plot(output.timeInc*par.VisualResolution,Etotal)
ylabel("Energy [J = kg m/s^2]")
xlabel("Time [s]")
legend('Etotal')
% potential energy
% rotational energy
% total energy
subdir = 'plots'
figName = [par.videoname '_totalEnergy'];
h2_fig = gcf;
if(par.writePdf)
    saveas(h2_fig,[subdir '\' figName '.fig']); saveas(h2_fig,[subdir '\' figName '.png']);
    print( h2_fig, [subdir '\' figName], ['-dpdf'] )
end


end
