# fuelCellChannelBypass
Set of codes developed to predict the effect of bypass channels on fuel cell performance. These codes were developed for the collaboration between Newcastle University (Dr Daniel Niblett) and University of New South Wales (Dr. Quentin Meyer).

1) PredictByPassMovement.m - This code evaluates an analytical prediction for the balance between viscous shear and intertial corner pressure drop for a single fuel cell serpentine channel. It compares the predictions from the analytical model to OpenFOAM, VOF CFD.

2) halfUnitCellDiffusivity.m - This code evaluates an analytical prediction for the reduction in effective diffusivity of a half-rib, half-channel unit cell due to the presence of a condensing water film under the rib.

3) fuelCellModel0D.m - This code predicts the performance of a fuel cell in the presence of a water film under the rib regions using code (2). The polarisation curve is reproduced in this scenario matching experimental conditions to predict the reason for the improved performance of flowfield channel bypass.
