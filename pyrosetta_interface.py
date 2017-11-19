#!/usr/bin/env python2.7

from pymol import cmd

import pyrosetta
from pyrosetta import rosetta

pyrosetta.init()
poses = {}

def update_structure(pose_name):
    '''Update the structure of the pose.'''
    pdb_os = rosetta.std.ostringstream()
    poses[pose_name].dump_pdb(pdb_os)
    cmd.delete(pose_name)
    cmd.read_pdbstr(pdb_os.str(), pose_name)

def load_pose(name, pdb_file):
    '''Load a pdb_file into a pose'''
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, pdb_file)
    poses[name] = pose

    update_structure(name)
    
def score_pose(name):
    '''Score a pose'''
    sfxn = rosetta.core.scoring.get_score_function()
    sfxn(poses[name])

def color_pose_by_energy(name):
    '''Color the pose according to the energy.'''
    score_pose(name)
    energies = [poses[name].energies().residue_total_energy(i) for i in range(1, poses[name].size() + 1)] ###DEBUG
    emax = 1.0 * max(energies)
    emin = 1.0 * min(energies)
    
    for i, e in enumerate(energies):
        cmd.set_color('res{0}_color'.format(i + 1), [ (e - emin) / (emax - emin), 0, (emax - e) / (emax - emin) ])
        cmd.color('res{0}_color'.format(i + 1), 'res {0}'.format(i + 1))

def selection_to_residue_selector(pose_name, selection):
    '''Get a residue selector corresponding to a pymol selection.'''
    myspace = {'l': set()}
    cmd.iterate(selection ,'l.add((model, resi + chain))', space=myspace)

    res_string = ','.join(r[1] for r in myspace['l'] if r[0] == pose_name)
    return rosetta.core.select.residue_selector.ResidueIndexSelector(res_string)

def repack_selected_residues(pose_name, selection):
    '''Repack the selected residues.'''
    res_selector = selection_to_residue_selector(pose_name, selection)
    not_res_selector = selection_to_residue_selector(pose_name, 'not ({0})'.format(selection))
  
    task_repack = rosetta.core.pack.task.operation.OperateOnResidueSubset(
            rosetta.core.pack.task.operation.RestrictToRepackingRLT(), res_selector)
    task_not_repack = rosetta.core.pack.task.operation.OperateOnResidueSubset(
            rosetta.core.pack.task.operation.PreventRepackingRLT(), not_res_selector)

    task_factory = rosetta.core.pack.task.TaskFactory()
    task_factory.push_back(task_repack)
    task_factory.push_back(task_not_repack)

    packer = rosetta.protocols.simple_moves.PackRotamersMover()
    packer.task_factory(task_factory)
    packer.apply(poses[pose_name])

    update_structure(pose_name)

cmd.extend("load_pose", load_pose)
cmd.extend("score_pose", score_pose)
cmd.extend("color_pose_by_energy", color_pose_by_energy)
cmd.extend("repack_selected_residues", repack_selected_residues)


###TEST
load_pose('foo', 'demo/1arb.pdb')
#load_pose('foo', 'demo/short.pdb')

#selection_to_residue_selector('(all)')
#repack_selected_residues('(all)')

#score_pose('foo')
#color_pose_by_energy('foo')
