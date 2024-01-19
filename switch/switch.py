import numpy as np
from MDAnalysis import Universe
from groio import write_gro
import groio
import matplotlib.pyplot as plt
from MDAnalysis.analysis.distances import distance_array as get_dist
import sys



try:
    part = int(sys.argv[1])
except:
    print(f'ERROR: Could not load {sys.argv[1]}. Exiting.')
    exit(0)




fraction = 0.5
polyDispersity = 13

sig0   = 67.31e-9 # m
sig30S = (14e-9/sig0)
sig50S = (17e-9/sig0)
sig70S = (20e-9/sig0)


m0   = 649*5000 # Da
mass30S = 0.855e6/m0
mass50S = 1.39e6/m0
mass70S = 2.30e6/m0

file = open("count.log", "a")
factor = 10
cutoff_30_70 = 1.5*(sig30S + sig70S)
cutoff_50_70 = 1.5*(sig50S + sig70S)
cut_off = 0.12


traj = Universe(f'sd{part}.gro', f'sd{part}.gro')



pos = traj.trajectory[0].positions*0.1

n30S = traj.select_atoms('name R30 or name R31 or name R32').n_atoms
n50S = traj.select_atoms('name R50 or name R51 or name R52').n_atoms
n70S = int(traj.select_atoms('name R70').n_atoms/polyDispersity)

indices = traj.select_atoms('all').indices
idx_30S = traj.select_atoms('name R30 or name R31 or name R32').indices
idx_50S = traj.select_atoms('name R50 or name R51 or name R52').indices
idx_70S = traj.select_atoms('name R70').indices
idx_70S = idx_70S.reshape((n70S, polyDispersity))


select = np.random.random(size = int(n70S))


select[select<0.01] = 0
select[select>0.01] = 1
pos_poly = pos[traj.select_atoms('name R70').indices].reshape((n70S, polyDispersity, 3))
pos_30S = pos[traj.select_atoms('name R30 or name R31 or name R32').indices]
pos_50S = pos[traj.select_atoms('name R50 or name R51 or name R52').indices]

box = traj.dimensions*0.1

lx = box[0]
ly = box[1]
lz = box[2]



ctr = 0
for s, _poly, idx in zip(select, pos_poly, idx_70S):
    
    if(not s):
        continue
    if(_poly[:, 0].any() > lx or _poly[:, 0].any() < 0.0 or _poly[:, 1].any() > ly or _poly[:, 1].any() < 0.0 or _poly[:, 2].any() > lz or _poly[:, 2].any() < 0.0):
        continue

    sw_pos = np.copy(pos)
    dist30S = get_dist(pos_30S, _poly[0])[:, 0]
    dist50S = get_dist(pos_50S, _poly[0])[:, 0]
    #_poly[0] for the 0th monomer of the polysome
    
    if(dist30S.min() > factor*cutoff_30_70 or dist50S.min() > factor*cutoff_50_70):
        print('No ribsosome nearby')
        continue
    
    switch30S = idx_30S[dist30S.argmin()]
    switch50S = idx_50S[dist50S.argmin()]
    switch70S = idx[-1]#switch70S is the 12th monomer of the polysome
    
    com_mono = (mass30S*sw_pos[switch30S] + mass50S*sw_pos[switch50S])/(mass30S + mass50S)
    #center of mass of the selected 30s and 50s 
    com_poly = np.copy(sw_pos[switch70S])#copy the position of the 12th monomer of the switching polysome
    sw_pos[switch70S] = np.copy(com_mono)#12th monomer of the polysome replaced by the centre of mass of the 30s and 50s
    sw_pos[switch30S] = sw_pos[switch30S] - com_mono + com_poly
    sw_pos[switch50S] = sw_pos[switch50S] - com_mono + com_poly
    if(sw_pos[switch70S][0] > lx or sw_pos[switch70S][0] < 0.0 or sw_pos[switch30S][1] > ly or sw_pos[switch30S][1] < 0.0 or sw_pos[switch50S][2] > lz or sw_pos[switch50S][2] < 0.0):
    	continue
    
    _pos = np.copy(sw_pos[:idx[0]])
    _pos = np.append(_pos, sw_pos[switch70S].reshape((1, 3)), axis=0)
    _pos = np.append(_pos, sw_pos[idx[0]:switch70S], axis=0)
    _pos = np.append(_pos, sw_pos[switch70S+1:], axis=0)
    dummy_pos = np.copy(_pos)

    # dist_pos = get_dist(dummy_pos, dummy_pos,  backend='OpenMP')
    # np.fill_diagonal(dist_pos, dist_pos.max())
    # if np.any(dist_pos <= cut_off):
    #     continue

    
    dist_pos1 = get_dist(dummy_pos, dummy_pos[switch30S])
    dist_pos1[switch30S, 0] = dist_pos1.max()
    
    dist_pos2 = get_dist(dummy_pos, dummy_pos[switch50S])
    dist_pos2[switch50S, 0] = dist_pos2.max()
    
    index = int(switch70S - (polyDispersity-1))
    dist_pos3 = get_dist(dummy_pos, dummy_pos[index])
    dist_pos3[index, 0] = dist_pos3.max()

    dist_pos4 = get_dist(dummy_pos, dummy_pos[switch70S])
    dist_pos4[switch70S, 0] = dist_pos4.max()

    if(dist_pos1.min()<=cut_off or dist_pos2.min()<=cut_off or dist_pos3.min()<=cut_off or dist_pos4.min()<=cut_off):
        continue

    pos = np.copy(dummy_pos)
    ctr += 1

write_gro(f'next{part+1}.gro', pos, traj.select_atoms('all').resids.tolist(), traj.select_atoms('all').resnames.tolist(),
          traj.select_atoms('all').names.tolist(), traj.trajectory[0]._unitcell[:3]*0.1)


file.write("%d\t%d\n" % (len(select[select==1]), ctr))
    

