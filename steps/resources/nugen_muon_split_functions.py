from icecube import dataclasses, icetray, dataio
from I3Tray import *
import numpy as np

from resources.geometry import distance_to_icecube_hull
from resources.muon_split_functions import get_delta_psi_vector, get_delta_psi_vector_dir, particle_is_inside
from ic3_labels.labels.utils import geometry
from scipy.spatial import ConvexHull




def build_tree_with_muon_split_nugen(frame, new_psi, random_service, beta=1):
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
    beta : number to scale sampled angle
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
            
            d_type = d.type_string
            # Reject energy losses which are not NuclInt or Brems, ~2%
            if ((d_type != 'NuclInt') & (d_type != 'Brems')):
                # print('Type: ', d_type)
                return False 
            
            
            ### For linear analysis
            if new_psi != 0: 
                # Incoming muon energy
                E = frame['I3MapSplit']['E_before_max_E_loss']
                # Muon energy afert maximum energy loss
                # E_ = E - d.energy 
                E_ = E - frame['I3MapSplit']['max_E_loss']
                # Get angle in rad
                if d_type == 'NuclInt':
                    new_psi = get_new_psi_nuclint(E, E_, random_service) * beta
                elif d_type == 'Brems':
                    new_psi = get_new_psi_brems(E, E_, random_service) * beta
                
            new_zenith_azimuth = get_delta_psi_vector(muon.dir.zenith, muon.dir.azimuth, random_service=random_service, delta_psi=new_psi, is_degree=True)

            # Insert key to frame
            # frame.Put('I3MapSplit', dataclasses.I3MapStringDouble())
            i3map = frame['I3MapSplit']
            # Save opening angle and new direction
            i3map['delta_psi_in_degree'] = new_psi
            i3map['incoming_muon_energy'] = muon.energy
            i3map['E_after_max_E_loss'] = E_
            if d_type == 'NuclInt':    
                i3map['max_E_loss_type'] = 1 # -2000001004
            elif d_type == 'Brems':
                i3map['max_E_loss_type'] = 2 # -2000001001
            else:
                i3map['max_E_loss_type'] = 0
            i3map['delta_psi_in_degree'] = new_psi 
            i3map['new_zenith'] = new_zenith_azimuth[0][0]
            i3map['new_azimuth'] = new_zenith_azimuth[1][0]
            i3map['old_zenith'] = d.dir.zenith
            i3map['old_azimuth'] = d.dir.azimuth
            i3map['time'] = d.time
            i3map['pre_length'] = (tree[0].pos - tree[max_index].pos).magnitude
            i3map['post_length'] = (tree[max_index].pos - tree[-1].pos).magnitude

            # Insert pos key to frame
            frame.Put('I3PosOfMaxEnergyLoss', dataclasses.I3Position(d.pos))

        if found_max == True:
            # Set new direction
            tree[i].dir = dataclasses.I3Direction(new_zenith_azimuth[0][0], new_zenith_azimuth[1][0])

            # Set new position
            tree[i].pos = muon_split_helper.pos + muon_split_helper.dir * (muon_split_helper.pos - d.pos).magnitude
            
            
def get_new_psi_brems(E, E_, rnd_state, is_degree=True, theta_star=1):
    epsilon = E - E_
    mu = 0.1057 # muon mass
    p = rnd_state.uniform(0, 1)
    r_max = np.minimum(1, E_/epsilon) * E * theta_star / mu
    a = p * r_max**2 / (1+r_max**2)
    r = np.sqrt(a/(1-a))
    theta_photon = mu / E * r
    theta_mu = epsilon / E_ * theta_photon
    
    if is_degree:
        # print(np.rad2deg(theta_mu))
        return np.rad2deg(theta_mu)
    else:
        return theta_mu
    
def get_new_psi_nuclint(E, E_, rnd_state, is_degree=True):
    M = 0.9383 # Proton mass in GeV
    mu = 0.1057 # Muon mass in GeV
    m_0 = np.sqrt(0.4)
    p = rnd_state.uniform(0, 1)
    # nu = epsilon
    epsilon = E - E_
    y = epsilon / E
    t_max = 2 * M * epsilon
    t_min = (mu * y)**2 / (1 - y)
    t_1 = np.minimum(epsilon**2, m_0**2)
    t_p = (t_max * t_1) / ((t_max + t_1) * ((t_max * (t_min + t_1))\
                    / (t_min * (t_max + t_1)))**p - t_max)
    sin2 = (t_p - t_min) / (4 * (E * E_ - mu**2) - 2 * t_min)
    theta_mu = 2 * np.arcsin(np.sqrt(sin2))
    
    if is_degree:
        # print(np.rad2deg(theta_mu))
        return np.rad2deg(theta_mu)
    else:
        return theta_mu
            
            
def get_max_energy_loss_id_nugen(tree, muon_id, energy_loss, min_dist):
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
        while tree[min_id+muon_id].minor_id > (tree[muon_id].minor_id + 1) :
            min_id += 1
        if min_id == 1:
            # There is no muon track
            # print('no muon track')
            return [0, 0, 0]
        muon_track = tree[muon_id+1:min_id+muon_id]   
    except:
        # print('no muon track')
        return [0, 0, 0]
    
    max_index_track = -1
    for i,d in enumerate(muon_track):
        if d.type_string != tree[muon_id].type_string:
            continue
        Next = 2
        if i > len(muon_track) - 3:
            # print('no maximum energy loss of 40% found (100m inner)')
            return [0, 0, 0]
        if muon_track[i+Next].type_string != tree[muon_id].type_string:
            print('next particle type not equal muon type')
            Next = 1
        if d.energy < 0.000001:
            print('energy equal zero: ', d.energy)
            return [0, 0, 0]
        if ((d.energy - muon_track[i+Next].energy)/ d.energy) > energy_loss:
            if distance_to_icecube_hull(muon_track[i+1].pos) > min_dist:
                continue
            E_before_max_E_loss = d.energy
            max_E_loss = E_before_max_E_loss - muon_track[i+Next].energy
            hull_distance_of_max_e_loss = distance_to_icecube_hull(muon_track[i+1].pos)
            max_index_track = i + 1
    
    if max_index_track == -1:
        # print('no maximum energy loss of 40% found (100m inner)')
        return [0, 0, 0]
    
    
    # distances = [distance_to_icecube_hull(d.pos) for d in muon_track]             
    # energies = [d.energy if d.type_string == tree[muon_id].type_string and distances[l]<0 else 0 for d,l in zip(muon_track, range(len(distances)))]
    # max_E_loss = max(energies)
    # max_index_track = energies.index(max_E_loss)               
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



def selection_nugen(self, frame):
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
    try:
        while ((tree[muon_id].type_string != 'MuMinus') & (tree[muon_id].type_string != 'MuPlus')):
            muon_id += 1
            muon = tree[muon_id]
    except:
        # There is no muon in tree
        # print('no muon in tree')
        return False
        
      
    # Check if particle moves through detector
    if particle_is_inside(muon, self._convex_hull) != True:
        # print('particle is not going through detector')
        return False

    # Save output of get_max_energy_loss_id (id, energy, negative distance to hull)
    id_E_dist = get_max_energy_loss_id_nugen(tree, muon_id, self._percentage_energy_loss, self._min_dist)

    # Check if muon is going through detector --> about 60% of the events are rejected
    if id_E_dist[1] == 0:
        return False
    
    # Check for minimum energy loss: 30% loss -> 3% left, 50% loss -> 1.3% left
    # min_energy_loss = self._percentage_energy_loss * muon.energy
    # if min_energy_loss > id_E_dist[1]:
        # print('energy loss is too low')
        # return False
    
    # Energy loss should be in the middle of the detector
    # hull_distance_of_max_e_loss = id_E_dist[2] # this value is negative
    # if self._min_dist < hull_distance_of_max_e_loss:
        # print('energy loss not in the middle of detector')
        # return False
    
    # Save data in tree 
    frame.Put('I3MapSplit', dataclasses.I3MapStringDouble())
    frame['I3MapSplit']['muon_id'] = muon_id
    frame['I3MapSplit']['max_E_id'] = id_E_dist[0]
    frame['I3MapSplit']['max_E_loss'] = id_E_dist[1]
    frame['I3MapSplit']['hull_dist'] = id_E_dist[2]
    frame['I3MapSplit']['E_before_max_E_loss'] = id_E_dist[3]