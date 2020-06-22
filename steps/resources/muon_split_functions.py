from icecube import dataclasses, icetray, dataio
from I3Tray import *
import numpy as np

from egenerator.utils.angles import get_delta_psi_vector, get_delta_psi_vector_dir

from resources.geometry import distance_to_icecube_hull

from ic3_labels.labels.utils.general import particle_is_inside 

from scipy.spatial import ConvexHull


def get_max_energy_loss_id(tree):
    '''Find maximum energy loss in the tree.

    Paramters
    ---------
    tree : I3MCTree

    Returns
    -------
    [max_index, max_e_loss, hull_distance_of_max_e_loss] : list of     
        index of maximum energy loss, 
        value of maximum energy loss, 
        distance of that position to the hull (negative, if inside the detector)

    '''
    # Create energy array with two energies for the first two particles (unkown and muon)
    energies = [0, 0]
    distances = [0, 0]
    
    for d in tree[2:]:
        # Stashes the muons at every point and searches for the maximum energy inside the detector
        hull_distance = distance_to_icecube_hull(d.pos)
        distances.append(hull_distance)
        if d.type_string != 'MuMinus' and hull_distance < 0:
            energies.append(d.energy)
        else:
            energies.append(0)
    
    max_e_loss = max(energies)
    max_index = energies.index(max_e_loss)   
    # If max_e_loss = 0 --> muon was not going through the detector
    hull_distance_of_max_e_loss = distances[max_index] 
    return [max_index, max_e_loss, hull_distance_of_max_e_loss]


def build_tree_with_muon_split(frame, new_psi, random_seed):
    '''Change direction of particles after maximum energy loss happens. First save the old tree, the existing tree will 
    be deleted and replaced with a new tree.

    Paramters
    ---------
    frame : I3Frame
    new_psi : angle in degree
        Set new direction of following particles of the maximum energy loss.
    random_seed : int
        Random service for generating the new direction.
    '''

    tree = frame['I3MCTree']

    # Create a tree and initialize it with the old I3MCTree 
    frame.Put('OldTree', dataclasses.I3MCTree(tree))
    oldTree = frame['OldTree']

    # Set muon and primary
    primary = tree[0]
    muon = tree[1]

    # Empty new tree except unknown and muon
    # tree.erase_children(muon)

    # Delete tree
    frame.Delete('I3MCTree')

    # Create new I3MCTree
    frame.Put('I3MCTree', dataclasses.I3MCTree())
    
    # Initialize tree with the new I3MCTree and set primary and muon
    tree = frame['I3MCTree']
    tree.add_primary(primary)
    tree.append_child(primary.id, muon)


    # Get index of particle with highest energy loss
    # max_index = get_max_energy_loss_id(oldTree)[0]
    max_index = int(frame['I3MapSplit']['max_E_id'])
    
    # Create two variables to fill the tree
    found_max = False
    daughter_counter = 0
    daughter_direction_change = False

    # Loop through particles in oldTree except the first two particles because they already exist in tree
    for d in oldTree[2:]:

        # Check for maximum energy loss
        if d == oldTree[max_index]:
            # print('maximum energy loss: {}'.format(oldTree[max_index]))
            
            # Set parameter: maximum loss found, now muon2 has to be inserted, change direction of the following daughter particles
            found_max = True
            daughter_direction_change = True

            muon_split_helper = d

                       
            # Get new zenith and azimuth value for splitted particles with given delta_angle
            random_service = np.random.RandomState(random_seed)
            new_zenith_azimuth = get_delta_psi_vector(muon.dir.zenith, muon.dir.azimuth, random_service=random_service, delta_psi=new_psi, is_degree=True)

            # Insert key to frame
            # frame.Put('I3MapSplit', dataclasses.I3MapStringDouble())
            i3map = frame['I3MapSplit']
            # Save opening angle and new direction
            i3map['delta_psi_in_degree'] = new_psi 
            i3map['new_dir_x'] = new_zenith_azimuth[0][0]
            i3map['new_dir_y'] = new_zenith_azimuth[1][0]
            i3map['old_dir_x'] = d.dir.zenith
            i3map['old_dir_y'] = d.dir.azimuth
            i3map['time'] = d.time

            # Insert pos key to frame
            frame.Put('I3PosOfMaxEnergyLoss', dataclasses.I3Position(d.pos))

        # Insert daughter particles
        if daughter_counter != 0:
            # Check if this daughter particle is behind muon2
            if daughter_direction_change == True:
                 # Set new direction
                d.dir = dataclasses.I3Direction(new_zenith_azimuth[0][0], new_zenith_azimuth[1][0])
                # Set new position
                d.pos = muon_split_helper.pos + muon_split_helper.dir * (muon_split_helper.pos - d.pos).magnitude
                tree.append_child(daughter_id, d)
            else:
                tree.append_child(daughter_id, d)
            daughter_counter -= 1
            continue

        if found_max == False:
            # Append all particles until particle with max energy loss
            tree.append_child(muon.id, d)

            # Check if that particle has daughters
            daughter_counter = len(oldTree.get_daughters(d))
            if daughter_counter != 0:
                daughter_id = d.id

        else:
            # Append rest particles with new direction and position
            
            # Set new direction
            d.dir = dataclasses.I3Direction(new_zenith_azimuth[0][0], new_zenith_azimuth[1][0])

            # Set new position
            d.pos = muon_split_helper.pos + muon_split_helper.dir * (muon_split_helper.pos - d.pos).magnitude 
            
            tree.append_child(muon.id, d)
            

            # Check if that particle has daughters
            daughter_counter = len(oldTree.get_daughters(d))
            if daughter_counter != 0:
                daughter_id = d.id
                continue


    # print('OldTree: {}'.format(oldTree))
    # print('I3MCTree: {}'.format(tree))
    # print('----------- next one -------------')


def selection(self, frame):
    '''select events with a minimum energy loss in the detector.
    
    Parameters
    ----------
    frame : I3Frame
    percentage_energy_loss : float, between 0 and 1
        Defines the percentage of the maximum energy.

    Return
    ------
    Bool : False, if event does not fullfill the conditions
    '''

    tree = frame['I3MCTree']

    # Check if particle moves through detector
    if particle_is_inside(tree[1], self._convex_hull) != True:
        return False

    # Save output of get_max_energy_loss_id (id, energy, negative distance to hull)
    id_E_dist = get_max_energy_loss_id(tree)

    # Check if muon is going through detector --> about 60% of the events are rejected
    if id_E_dist[1] == 0:
        return False
    
    # Check for minimum energy loss: 30% loss -> 3% left, 50% loss -> 1.3% left
    min_energy_loss = self._percentage_energy_loss * tree[1].energy
    if min_energy_loss > id_E_dist[1]:
        return False
    
    # Energy loss should be in the middle of the detector
    hull_distance_of_max_e_loss = id_E_dist[2] # this value is negative
    min_distance_to_hull = -100
    if min_distance_to_hull < hull_distance_of_max_e_loss:
        return False
    
    # Save data in tree 
    frame.Put('I3MapSplit', dataclasses.I3MapStringDouble())
    frame['I3MapSplit']['max_E_id'] = id_E_dist[0]
    frame['I3MapSplit']['max_E_loss'] = id_E_dist[1]
    frame['I3MapSplit']['hull_dist'] = id_E_dist[2]
