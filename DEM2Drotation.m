function rot = DEM2Drotation(rad)
    rot = [cos(rad) -sin(rad); sin(rad) cos(rad)];
end