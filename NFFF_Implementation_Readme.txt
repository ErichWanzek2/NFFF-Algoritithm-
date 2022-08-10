The NFFF paramters and arrays are set up in the MATLAB script titled:

Initialize_FieldCalculator.m

This is where I put arrays to define the radial distances, freuqnecies,theta, and phi angle values 
for the NFFF to pefrom aalong  with the bounds for the rectangular Huygen Surface.



The NFFF calcualtor that calcualtes the NFFF lx,ly,lz,Nx,Ny,NZ for every time step is performed by
the MATLAB script:

FieldCalculators.m



The NFFF post processing and plotting is found as addtional MATLAB code added at the end of the MATLAB
script titled:

PostProcess.m