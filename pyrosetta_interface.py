#!/usr/bin/env python2.7

from pymol import cmd

import pyrosetta
from pyrosetta import rosetta

pyrosetta.init()
poses = {}


def load_pose(name, pdb_file):
    '''Load a pdb_file into a pose'''
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, pdb_file)
    poses[name] = pose

def score_pose(name):
    '''Score a pose'''
    sfxn = rosetta.core.scoring.get_score_function()
    sfxn(poses[name])

def color_pose_by_energy(name):
    '''Color the pose according to the energy.'''
    energies = [poses[name].energies().residue_total_energy(i) for i in range(1, poses[name].size() + 1)] ###DEBUG
    emax = 1.0 * max(energies)
    emin = 1.0 * min(energies)
    
    for i, e in enumerate(energies):
        cmd.set_color('res{0}_color'.format(i + 1), [ (e - emin) / (emax - emin), 0, (emax - e) / (emax - emin) ])
        cmd.color('res{0}_color'.format(i + 1), 'res {0}'.format(i + 1))


###TEST
load_pose('foo', 'demo/1arb.pdb')

pdb_os = rosetta.std.ostringstream()
poses['foo'].dump_pdb(pdb_os)
cmd.read_pdbstr(pdb_os.str(), 'foo')
score_pose('foo')
color_pose_by_energy('foo')
