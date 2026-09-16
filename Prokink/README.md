# Prokink
New python implementation of Prokink. Based on paper: https://doi.org/10.1093/protein/13.9.603
Fast and much easier to use than original Fortran code. I have used it a few times and results look reasonable, but IT IS NOT OFFICIALLY TESTED. 

Requires a python environment with numpy and mdtraj loaded. It is easy to build your own, but you can also use this one:
conda activate /home/agp2004/anaconda3/envs/tica_env

## To use:

# INPUTS. Either:
1) psf and dcd
2) pdb

# USER INPUT parameters:
proline_resid = resid of the proline. This is the resid (as in VMD atomselect) in the psf file NOT THE RESIDUE NUMBER

pre_pro_resid1 = resid of the first residue of the pre-proline helix part

pre_pro_resid2 = resid of the last residue of the pre-proline helix part

post_pro_resid1 = resid of the first residue of the post-proline helix part

post_pro_resid2 = resid of the last residue of the post-proline helix part

That's it. The program will take care of the rest for you. I implemented a checker that should stop you from using it incorrectly. You can turn it off if you want with checker=False.

Note that the original implemenatation required manual input of files including things like helix axis, etc. So the results may not exactly match the original, but they should be close. 

# OUTPUTS: Either:
1) Two vectors of length n_frames (number of frames in input dcd) containing the bend angles and wobble angles of each frame of the input dcd

2) Two number representing the bend angle and wobble angle of the input pdb

# EXAMPLE USE calculating the bend and wobble angles of TM6 over a trajectory of the 5HT2A:

import sys

sys.path.append('/athena/hwlab/scratch/agp2004/ProKink')

from prokink import *

psf = '/athena/hwlab/scratch/agp2004/NEW_5HT2A/no_palm_all_traj/all_traj/5HT/0/ionized.psf'

dcd = '/athena/hwlab/scratch/agp2004/NEW_5HT2A/no_palm_all_traj/all_traj/5HT/0/5HT_0.wrapped.aligned.dcd'

proline_resid = 338

pre_pro_resid1 = 315

pre_pro_resid2 = 336

post_pro_resid1 = 339

post_pro_resid2 = 349

bend_angles, wobble_angles = prokink(proline_resid, pre_pro_resid1, pre_pro_resid2, post_pro_resid1, post_pro_resid2, pdb=False, psf=psf, dcd=dcd)
