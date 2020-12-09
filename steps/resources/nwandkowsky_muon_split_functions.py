from icecube import dataclasses, icetray, dataio
from I3Tray import *
import numpy as np

from resources.geometry import distance_to_icecube_hull
from resources.muon_split_functions import get_delta_psi_vector, get_delta_psi_vector_dir, particle_is_inside
from ic3_labels.labels.utils import geometry
from scipy.spatial import ConvexHull




def build_tree_with_muon_split(frame, new_psi, random_service):
    '''Change direction of particles after maximum energy loss happens. First save the old tree, 
    the existing tree will 
    be deleted and replaced with a new tree.

    Paramters
    ---------
    frame : I3Frame
    new_psi : angle in degree
        Set new direction of following particles of the maximum energy loss.
    random_service : np.random.RandomState
        Random service for generating the new direction.
    '''

    
    tree = frame['I3MCTree']

    # Create a tree and initialize it with the old I3MCTree 
    frame.Put('OldTree', dataclasses.I3MCTree(tree))

    # Get muon
    muon_id = int(frame['I3MapSplit']['muon_id'])
    muon = tree[muon_id]

    # Get index of maximum energy loss    
    max_index = int(frame['I3MapSplit']['max_E_id'])

    # Create variables to fill the tree
    found_max = False

    
    # Loop through particles in oldTree
    tree_len = len(tree[muon_id:])
    for i in range(tree_len):
        # Skip particles until muon
        i += muon_id 
        # Get particle
        d = tree[i]
        # Check if muon track ends
        if d.time < 0:
            break
        # Check for maximum energy loss
        if d == tree[max_index]:
            # print('maximum energy loss: {}'.format(oldTree[max_index]))
            
            # Set parameter: maximum loss found, now muon2 has to be inserted, change direction of the following daughter particles
            found_max = True
            muon_split_helper = d
               
            # Get new zenith and azimuth value for splitted particles with given delta_angle
            # random_service = np.random.RandomState(random_seed)
            
            ### For linear analysis
            if new_psi != 0:    
                new_psi = np.exp(random_service.uniform(np.log(0.1), np.log(new_psi))) 
            new_zenith_azimuth = get_delta_psi_vector(muon.dir.zenith, muon.dir.azimuth, random_service=random_service, delta_psi=new_psi, is_degree=True)

            # Insert key to frame
            # frame.Put('I3MapSplit', dataclasses.I3MapStringDouble())
            i3map = frame['I3MapSplit']
            # Save opening angle and new direction
            i3map['delta_psi_in_degree'] = new_psi
            # i3map['incoming_muon_energy'] = muon.energy
            if d.type_string == 'NuclInt':    
                i3map['max_E_loss_type'] = -2000001004
            elif d.type_string == 'Brems':
                i3map['max_E_loss_type'] = -2000001001
            else:
                i3map['max_E_loss_type'] = 0
            i3map['new_zenith'] = new_zenith_azimuth[0][0]
            i3map['new_azimuth'] = new_zenith_azimuth[1][0]
            i3map['old_zenith'] = d.dir.zenith
            i3map['old_azimuth'] = d.dir.azimuth
            i3map['time'] = d.time
            i3map['pre_length'] = (tree[0].pos - tree[max_index].pos).magnitude
            i3map['post_length'] = (tree[max_index].pos - tree[-1].pos).magnitude

            # Insert pos key to frame
            # frame.Put('I3PosOfMaxEnergyLoss', dataclasses.I3Position(d.pos))

        if found_max == True:
            # Set new direction
            tree[i].dir = dataclasses.I3Direction(new_zenith_azimuth[0][0], new_zenith_azimuth[1][0])

            # Set new position
            tree[i].pos = muon_split_helper.pos + muon_split_helper.dir * (muon_split_helper.pos - d.pos).magnitude
            
            
            
def get_max_energy_loss_id_interaction(tree, muon_id, percentage_energy_loss, min_dist):
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
           
    # Get muon track               
    min_id = 1
    try:
        while tree[min_id+muon_id].minor_id < tree[muon_id].minor_id:
            min_id += 1
        if min_id == 1:
            # There is no muon track
            return [0, 0, 0] 
        muon_track = tree[muon_id+1:min_id+muon_id]   
    except:
        return [0, 0, 0]
    # distances = [distance_to_icecube_hull(d.pos) for d in muon_track]               
                   
    # energies = [d.energy if distances[l]<0 else 0 for d,l in zip(muon_track, range(len(distances)))]
    # max_E_loss = max(energies)
    # max_index_track = energies.index(max_E_loss)               
    # max_index = max_index_track + muon_id + 1 # Necessary to get access via tree               
    
    
    incoming_muon_energy = tree[muon_id].energy
    max_index_track = -1
    for i, d in enumerate(muon_track):
        # Check for daughter particles, otherwise the losses double
        if d.length != 0:
            continue
        # Check energy loss
        if (d.energy / incoming_muon_energy) > percentage_energy_loss:
            # Check position of energy loss
            if distance_to_icecube_hull(d.pos) > min_dist:
                break
            E_before_max_E_loss = incoming_muon_energy
            hull_distance_of_max_e_loss = distance_to_icecube_hull(d.pos)
            max_E_loss = d.energy
            max_index_track = i
            break
        incoming_muon_energy -= d.energy
    if max_index_track == -1:
        return [0, 0, 0]
    max_index = max_index_track + muon_id + 1 # Necessary to get access via tree
    
    
    '''Long loop instead of list comprehension
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
    '''
    
    
       
    # If max_e_loss = 0 --> muon was not going through the detector
    # hull_distance_of_max_e_loss = distances[max_index_track] 
    return [max_index, max_E_loss, hull_distance_of_max_e_loss, E_before_max_E_loss]



def get_max_energy_loss(tree, muon_id):
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
    # Get muon track               
    min_id = 1
#    try:
    tree_len = len(tree)
    # Check if muon is the last particle in tree
    if muon_id+1 >= tree_len:
        return [-1, 0, 0]
    while tree[min_id+muon_id].minor_id < tree[muon_id].minor_id:
        min_id += 1
        if min_id+muon_id >= tree_len:
            break
            
    if min_id == 1:
            # There is no muon track
        print('There is no muon track')
        return [-1, 0, 0] 
    
    muon_track = tree[muon_id+1:min_id+muon_id]   
    
#    except:
        # return [-1, 0, 0]
    
    distances = [distance_to_icecube_hull(d.pos) if d.length==0 else 0 for d in muon_track]               
                   
    energies = [d.energy if distances[l]<0 else 0 for d,l in zip(muon_track, range(len(distances)))]
    max_E_loss = max(energies)
    max_index_track = energies.index(max_E_loss)               
    max_index = max_index_track + muon_id + 1 # Necessary to get access via tree               
    # If max_e_loss = 0 --> muon was not going through the detector
    hull_distance_of_max_e_loss = distances[max_index_track]
    return [max_index, max_E_loss, hull_distance_of_max_e_loss]
    
'''    
    incoming_muon_energy = tree[muon_id].energy
    max_index_track = -1
    for i, d in enumerate(muon_track):
        # Check for daughter particles, otherwise the losses double
        if d.length != 0:
            continue
        # Check energy loss
        if (d.energy / incoming_muon_energy) > percentage_energy_loss:
            # Check position of energy loss
            if distance_to_icecube_hull(d.pos) > min_dist:
                break
            E_before_max_E_loss = incoming_muon_energy
            hull_distance_of_max_e_loss = distance_to_icecube_hull(d.pos)
            max_E_loss = d.energy
            max_index_track = i
            break
        incoming_muon_energy -= d.energy
    if max_index_track == -1:
        return [0, 0, 0]
    max_index = max_index_track + muon_id + 1 # Necessary to get access via tree
'''    
    
'''Long loop instead of list comprehension
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
'''



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
    
    # selection = 0: passes selection, 1: does not pass
    selection = 0

    # Get muon
    muon_id = 0
#    try:
    tree_len = len(tree)
    while ((tree[muon_id].type_string != 'MuMinus') & (tree[muon_id].type_string != 'MuPlus')):
        muon_id += 1
        # muon = tree[muon_id]
        if muon_id >= tree_len:
            return False
    muon = tree[muon_id]
        
#    except:
        # Muon not in tree
        # selection = 1
    # if (muon.type_string != 'MuMinus') & (muon.type_string != 'MuPlus'):
        # print('No muon found in tree (selection)-----ERROR')
        # return False
        
      
    # Check if particle moves through detector
    if particle_is_inside(muon, self._convex_hull) != True:
        selection = 1
        # return False

    # Save output of get_max_energy_loss_id (id, energy, negative distance to hull)
    # id_E_dist = get_max_energy_loss_id_interaction(tree, muon_id, self._percentage_energy_loss, self._min_dist)

    id_E_dist = get_max_energy_loss(tree, muon_id)
    
    # Check if muon is going through detector --> about 60% of the events are rejected
    if id_E_dist[1] == 0:
        # print('------you should never see this print!------')
        selection = 1
        # return False
    
    # Check for minimum energy loss: 30% loss -> 3% left, 50% loss -> 1.3% left
    min_energy_loss = self._percentage_energy_loss * muon.energy
    if min_energy_loss > id_E_dist[1]:
        selection = 1
        # return False
    
    # Energy loss should be in the middle of the detector
    hull_distance_of_max_e_loss = id_E_dist[2] # this value is negative
    if self._min_dist < hull_distance_of_max_e_loss:
        selection = 1
        # return False
        
        
    
    # Save data in tree 
    frame.Put('I3PosOfMaxEnergyLoss', dataclasses.I3Position(tree[id_E_dist[0]].pos))
    frame.Put('I3MapSplit', dataclasses.I3MapStringDouble())
    frame['I3MapSplit']['muon_id'] = muon_id
    frame['I3MapSplit']['incoming_muon_energy'] = muon.energy
    frame['I3MapSplit']['max_E_id'] = id_E_dist[0]
    frame['I3MapSplit']['max_E_loss'] = id_E_dist[1]
    frame['I3MapSplit']['hull_dist'] = id_E_dist[2]
    frame['I3MapSplit']['selection'] = selection # pass selection?
    
    return True


def selection_MC(self, frame):
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

    # Get muon
    muon_id = 0
    tree_len = len(tree)
#     try:
    while ((tree[muon_id].type_string != 'MuMinus') & (tree[muon_id].type_string != 'MuPlus')):
        muon_id += 1
        if muon_id >= tree_len:
            return False
    
    
    muon = tree[muon_id]
#     except:
        # Muon not in tree
    # if (muon.type_string != 'MuMinus') & (muon.type_string != 'MuPlus'):
        # print('No muon found in tree (selection)-----ERROR')
        # return False
        
      
    # Check if particle moves through detector, only checks the direction
    if particle_is_inside(muon, self._convex_hull) != True:
        return False

    # Save output of get_max_energy_loss_id (id, energy, negative distance to hull)
    # id_E_dist = get_max_energy_loss_id_interaction(tree, muon_id, self._percentage_energy_loss, self._min_dist)

    id_E_dist = get_max_energy_loss(tree, muon_id)
    
    # Check if muon is going through detector --> about 60% of the events are rejected, checks if particle looses energy in detector
    if id_E_dist[1] == 0:
        # print('------you should never see this print!------')
        return False
    
    # Check for minimum energy loss: 30% loss -> 3% left, 50% loss -> 1.3% left
    min_energy_loss = self._percentage_energy_loss * muon.energy
    if min_energy_loss > id_E_dist[1]:
        return False
    
    # Energy loss should be in the middle of the detector
    hull_distance_of_max_e_loss = id_E_dist[2] # this value is negative
    if self._min_dist < hull_distance_of_max_e_loss:
        return False
    
    # Save data in tree 
    frame.Put('I3PosOfMaxEnergyLoss', dataclasses.I3Position(tree[id_E_dist[0]].pos))
    frame.Put('I3MapSplit', dataclasses.I3MapStringDouble())
    frame['I3MapSplit']['muon_id'] = muon_id
    frame['I3MapSplit']['incoming_muon_energy'] = muon.energy
    frame['I3MapSplit']['max_E_id'] = id_E_dist[0]
    frame['I3MapSplit']['max_E_loss'] = id_E_dist[1]
    frame['I3MapSplit']['hull_dist'] = id_E_dist[2]
    frame['I3MapSplit']['selection'] = 0 # passes selection

   
    return True
    
    
def cut_millipede_out_of_detector(frame):
    '''Add new key to frame with list of energy losses only inside the detector.
    '''
    if 'SplineMPE_MillipedeHighEnergyMIE' not in frame.keys():
        print('No millipede in frame')
        return False
    millipede = frame['SplineMPE_MillipedeHighEnergyMIE']
    millipede_in_detector = [d for d in millipede if distance_to_icecube_hull(d.pos) < 0]
    if len(millipede_in_detector) == 0:
        p = dataclasses.I3Particle()
        p.energy = 0
        millipede_in_detector.append(p)
    frame.Put('SplineMPE_MillipedeInsideDetector', dataclasses.I3VectorI3Particle(millipede_in_detector))


