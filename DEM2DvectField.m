function F = DEM2DvectField(vx,vz,fx,fz,fwx,fwz,m,par)

F(1) = vx;
F(2) = vz;
F(3) = (sum(fx) + fwx)/m ;
F(4) = (sum(fz) + fwz)/m - par.g;


end