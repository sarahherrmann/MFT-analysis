FindTrackableTracks:

Trackable = track that left at least 1 cluster on at least 4 different MFT disks


Variables :
pt -> pt of the MC track (generated kinematic), //at the origin of the track ?
eta -> eta of the MC track (generated kinematic) //idem ?
phi -> phi of the MC track (generated kinematic), //at the origin of the track ?
zVtx -> position of the primary vertex (as used in the MC simulation) on the z axis


x -> origin of the MC generated track on the x axis
y -> origin of the MC generateds track on the y axis
z -> origin of the MC generateds track on the z axis
R=sqrt(x^2+y^2) -> radius of the production vertex point (only for MC)

Histograms:
  1D number of MClabels found with each occurence, as a function of this number of occurences

  2D pt vs Eta of the trackables
  2D phi vs eta of the trackables
  2D phi vs pt of the trackables
  2D zVtx vs Eta of the trackables
  2D R vs zVtx of the trackables
