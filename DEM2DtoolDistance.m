function [partDist, unfrozenIndex] = DEM2DtoolDistance(par,data)
    toolCenterMass = 0.5*(data.toolbBox(1,:) + data.toolbBox(2,:));
    distVector = [data.position(1,:) - toolCenterMass(1);data.position(2,:) - toolCenterMass(2)];
    partDist = sqrt(distVector(1,:).^2 + distVector(2,:).^2);
    unfrozenIndex = (sign(distVector(1,:)) == sign(par.toolSpeed(1,1)))';
end