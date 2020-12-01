function  DEM2Dsplit(par,data)

    data.contactsParticle.deactivated(:) = 0;
    data.contactsMerged.N = 0;
    data.contactsMerged.aggregateSize(:,:) = 0;
    data.contactsMerged.positionMerged(:,:) = 0;
    data.contactsMerged.velocityMerged(:,:) = 0;
    data.contactsMerged.angularMerged(:,:) = 0;
    data.contactsMerged.inertiaTensor (:,:) = 0;
    data.contactsMerged.mass = 0;
end