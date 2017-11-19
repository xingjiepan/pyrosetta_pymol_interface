#!/usr/bin/env python2.7

from pymol import cmd

import pyrosetta
from pyrosetta import rosetta

pyrosetta.init()
poses = {}

def load_pose_into_object(pose_name):
    '''Load the pose into the object.'''
    pdb_os = rosetta.std.ostringstream()
    poses[pose_name].dump_pdb(pdb_os)
    cmd.read_pdbstr(pdb_os.str(), pose_name)
    cmd.color('green', '{0} and (n. c*)'.format(pose_name))

def update_structure(pose_name):
    '''Update the structure of the pose.'''
    pose = poses[pose_name]
    for i in range(1, pose.size() + 1):
        for j in range(1, pose.residue(i).natoms() + 1):
            name = pose.residue(i).atom_name(j)
            xyz = pose.residue(i).xyz(j)

            cmd.alter_state(1, '{0} and res {1} and n. {2}'.format(pose_name, i, name),
                    '(x, y, z)=({0}, {1}, {2})'.format(xyz.x, xyz.y, xyz.z))

def load_pose(name, pdb_file):
    '''Load a pdb_file into a pose'''
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, pdb_file)
    poses[name] = pose

    load_pose_into_object(name)
    
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
        cmd.label('res {0} and n. ca'.format(i + 1), '"{0:.2f}"'.format(energies[i]))

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
#load_pose('foo', 'demo/1arb.pdb')
#load_pose('foo', 'demo/short.pdb')
#cmd.alter_state(1, 'n. ca and res 2', '(x, y, z)=(0, 0, 0)')

#selection_to_residue_selector('(all)')
#repack_selected_residues('(all)')

#score_pose('foo')
#color_pose_by_energy('foo')
