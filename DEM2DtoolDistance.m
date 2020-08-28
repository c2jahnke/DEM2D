function dist = DEM2DtoolDistance(par,data)
    toolCenterMass = 0.5*(data.toolbBox(1,:) + data.toolbBox(2,:));
    distVector = [data.position(1,:) - toolCenterMass(1);data.position(2,:) - toolCenterMass(2)];
    dist = sqrt(distVector(1,:).^2 + distVector(2,:).^2);
end