function rot = DEM2Drotation(phi)
    rot = [cos(phi) -sin(phi); sin(phi) cos(phi)];
end