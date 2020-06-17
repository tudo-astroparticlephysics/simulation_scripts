from icecube import dataclasses, icetray, dataio
from I3Tray import *
import numpy as np

from egenerator.utils.angles import get_delta_psi_vector, get_delta_psi_vector_dir
from resources.geometry import get_closest_point_on_edge, get_distance_to_edge, get_edge_intersection, distance_to_axis_aligned_Volume, distance_to_icecube_hull


### Function to get index of particle with the maximum energy loss ###

def get_max_energy_loss_id(tree):
    '''
        return: list of index of maximum energy loss, value of maximum energy loss and the distance of that position to the hull (negative, if inside the detector)

    '''
    # Create energy array with two energies for the first two particles (unkown and muon)
    energies = [0, 0]
    distances = [0, 0]
    
    for d in tree[2:]:
        # Stashes the muons at every point and searches for the maximum energy inside the detector
        hull_distance = distance_to_icecube_hull(d.pos)
        distances.append(hull_distance * (-1))
        if d.type_string != 'MuMinus' and hull_distance < 0:
            energies.append(d.energy)
        else:
            energies.append(0)
    
    max_index = energies.index(max(energies))   
    # If max_e_loss = 0 --> muon was not going through the detector
    max_e_loss = max(energies) 
    hull_distance_of_max_e_loss = max(distances) * (-1)
    return [max_index, max_e_loss, hull_distance_of_max_e_loss]


### Function to split muon track at highest energy loss and insert this particle in a the tree, but first save the old tree, so the existing I3MCTree will be replaced with a new tree      

def build_tree_with_muon_split(frame, new_psi, random_seed):
    # new_psi: angle in degree to set new direction of muon2 and the following particles
    # random_seed: int to set random service
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
    max_index = get_max_energy_loss_id(oldTree)[0]
    
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
            frame.Put('I3MapMuonChangeDir', dataclasses.I3MapStringDouble())
            i3map = frame['I3MapMuonChangeDir']
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
        # !!!!! Was passiert, wenn ein Teilchen mehrere Tochterteilchen hat?
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

### Function to filter events with special conditions ###

def selection(frame):
    tree = frame['I3MCTree']

    # Save output of get_max_energy_loss_id (id, energy, negative distance to hull)
    id_E_dist = get_max_energy_loss_id(tree)

    # Check if muon is going through detector --> about 60% of the events are rejected
    if id_E_dist[1] == 0:
        return False

    # Check for minimum energy loss: 30% loss -> 3% left, 50% loss -> 1.3% left, 40% loss -> 2% left
    min_energy_loss = 0.01 * tree[1].energy
    # print('min_energy_loss: {}'.format(min_energy_loss))
    # print('max_energy_loss: {}'.format(get_max_energy_loss_id(tree)[1]))
    if min_energy_loss > id_E_dist[1]:
        return False
    
    # Energy loss should be in the middle of the detector
    hull_distance_of_max_e_loss = id_E_dist[2] # this value is negative
    min_distance_to_hull = -100
    if min_distance_to_hull < hull_distance_of_max_e_loss:
        return False

    
