# pCT

This will prepare the simulation proton data for the DROP-TVS pCT reconstruction algorithm!

Like with the pRad, you start with the **createTree.cc**, then you filter using either **filterSS.cc** or **filterDS.cc**, before building the binaries in **binary.cc** containing the proton data for the DROP-TVS to reconstruct.

**createTree.cc** <br />
Functions the same as in pRad, preparing the simulation output and forming a new root file to filter and build the binaries from.

**filterSS.cc** or **filterDS.cc**<br />
Only the head phantom and CTP phantoms are of interest in CT, so I have imported both the head phantom and CTP cylinder for use in hull-algorithm, just comment or uncomment the one you are interested in. Otherwise it functions pretty much in the same way as the pRad.

**binary.cc**<br />
The proton positions (from the hull and projected unto tracker plane positions) and wepl information is filled into some binary files and then DROP-TVS handles all the MLP, filtering, and reconstruction.
