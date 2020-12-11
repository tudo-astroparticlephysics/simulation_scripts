from icecube import dataclasses, icetray, dataio
from I3Tray import *
import numpy as np

from resources.geometry import distance_to_icecube_hull
from resources.muon_split_functions import get_delta_psi_vector, get_delta_psi_vector_dir, particle_is_inside
from ic3_labels.labels.utils import geometry
from ic3_labels.labels.utils.muon import get_muon_of_inice_neutrino
from ic3_labels.labels.utils import muon as mu_utils

from scipy.spatial import ConvexHull


########### NEW FUNCTIONS ###############################################
def selection(self, frame, perform_cut=False):
    '''Get primary muon and check:
        1. track hits detector
        2. maximum energy loss is higher than choosen E_loss
        3. check if maximum loss is at least x meters inside the detector
    '''
    
    # Get muon
    # tree = frame['I3MCTree']
    # muon_index = 0
    # while tree[muon_index].type_string not in ('MuMinus', 'MuPlus'):
        # muon_id += 1
        
    # Check selection, 1 = True
    selection = 1
    res = {
        'muon_minor_id': -1,
        'incoming_muon_energy': -1,
        'muon_energy_at_entry': -1,
        'max_p_minor_id': -1,
        'max_E_loss': -1,
        'hull_dist': -1,
        'E_before_max_loss': -1,
        'energy_before_max_loss': -1,
        'max_p_type': -1,
        'pre_length': -1,
        'post_length': -1,
    }
        
    muon = get_muon_of_inice_neutrino(frame)
    if muon == None:
        selection = 0
    else:
        # Check if particle moves through detector
        if particle_is_inside(muon, self._convex_hull) != True:
            selection = 0
        else:
            # Check maximum energy loss
            res_ = get_max_energy_loss(frame, muon, self._convex_hull)
            res.update(res_)
            entry, time, energy = mu_utils.get_muon_entry_info(frame, muon, self._convex_hull)
            res['muon_minor_id'] = muon.minor_id
            res['incoming_muon_energy'] = muon.energy
            res['muon_energy_at_entry'] = energy
            print('E at entry: ', res['muon_energy_at_entry'])
            if res['max_E_loss'] < res['muon_energy_at_entry'] * self._percentage_energy_loss:
                selection = 0
            elif res['hull_dist'] > self._min_dist:
                selection = 0
                
    res['selection'] = selection
    print(res)
    
    frame.Put('I3MapSplit', dataclasses.I3MapStringDouble(res))
    
    if perform_cut:
        return selection==1

    return True
    
    
    
    

def get_max_energy_loss(frame, muon, convex_hull):
    '''
        Return
        ------
        list of issues of maximum energy loss: minor_id, energy_loss, distance to hull
        and energy before maximum loss
    '''
    tree = frame['I3MCTree']
    muon_daughters = tree.get_daughters(muon)
    energy_losses_inside = [p.energy if (p.type_string not in ('MuMinus','MuPlus')) & (distance_to_icecube_hull(p.pos)<0) else 0 for p in muon_daughters]
    max_E_loss = float(np.max(energy_losses_inside))
    # if max_E_loss < muon.energy * min_energy_loss:
        # return None
        
    max_p_minor_id = -1
    hull_dist = -1
    E_before_max_loss = -1
    energy_before_max_loss = -1
    max_p_type = -1
    pre_length = -1
    post_length = -1 
    
    # Check if max_E_loss == 0 (max_E_loss is not in detector)
    if max_E_loss != 0:
        # Check min dist of max p
        max_index = energy_losses_inside.index(max_E_loss)
        max_p = muon_daughters[max_index]
        max_p_minor_id = max_p.minor_id
        hull_dist = distance_to_icecube_hull(max_p.pos)
        max_p_type = max_p.pdg_encoding
        # Get muon energy
        energy_before_max_loss = mu_utils.get_muon_energy_at_distance(frame, muon, (max_p.pos - muon.pos).magnitude - 1)

        if muon_daughters[max_index - 1].minor_id == max_p.minor_id-1:
            E_before_max_loss = muon_daughters[max_index-1].energy   
        elif muon_daughters[max_index - 2].minor_id == max_p.minor_id-1:
            E_before_max_loss = muon_daughters[max_index-2].energy
        elif muon_daughters[max_index - 3].minor_id == max_p.minor_id-1:
            E_before_max_loss = muon_daughters[max_index-3].energy
        else:        
            print('error in Energy before max loss')
            print(max_p.minor_id)
            print(muon_daughters[max_index - 1].minor_id)
            print(muon_daughters[max_index - 2].minor_id)
            print(muon_daughters[max_index - 3].minor_id)
            print(tree)
            
        # if distance_to_icecube_hull(max_p.pos) < min_dist:
    
        entry = mu_utils.get_muon_initial_point_inside(muon, convex_hull)
        exit = mu_utils.get_muon_exit_point(muon, convex_hull)
        if (entry == None) | (exit == None):
            pre_length = 0
            post_length = 0
        else:
            pre_length =  (entry - max_p.pos).magnitude
            post_length = (max_p.pos - exit).magnitude
         
                   
    return {
        'max_p_minor_id': max_p_minor_id,
        'max_E_loss': max_E_loss,
        'hull_dist': hull_dist,
        'E_before_max_loss': E_before_max_loss,
        'energy_before_max_loss': energy_before_max_loss,
        'max_p_type': max_p_type,
        'pre_length': pre_length,
        'post_length': post_length,
    }




def insert_deflection_angle(frame, new_psi, random_service, beta=1):

    max_p_minor_id = frame['I3MapSplit']['max_p_minor_id']
    if max_p_minor_id == -1:
        print('max_p_minor_id = -1')
        return True

    
    tree = frame['I3MCTree']
    
    frame.Put('OldTree', dataclasses.I3MCTree(tree))

    # muon_minor_id = frame['I3MapSplit']['muon_minor_id']
    
    found_max = False
    brems_nuclint = False
    
    
    # muon_daughters = tree.get_daughters(muon)
    
    
    # Problem: change all angles of tree, not only muon track angles
    for i in range(len(tree)):
        d = tree[i]
        if d.minor_id == max_p_minor_id:
            found_max = True
            muon_split_helper = d
               
            # Get new zenith and azimuth value for splitted particles with given delta_angle
            # random_service = np.random.RandomState(random_seed)
            
            d_type = d.type_string
            # Break energy losses which are not NuclInt or Brems, ~2%
            # if d_type in ('Brems', 'NuclInt'):            
                # brems_nuclint = True
            if new_psi != 0: 
                # Incoming muon energy
                E = frame['I3MapSplit']['E_before_max_loss']
                # Muon energy afert maximum energy loss
                # E_ = E - d.energy 
                E_ = E - frame['I3MapSplit']['max_E_loss']
                # Get angle in rad
                if d_type == 'NuclInt':
                    brems_nuclint = True
                    new_psi = get_new_psi_nuclint(E, E_, random_service) * beta
                elif d_type == 'Brems':
                    brems_nuclint = True
                    new_psi = get_new_psi_brems(E, E_, random_service) * beta
                else:
                    new_psi = 0
                
            new_zenith_azimuth = get_delta_psi_vector(d.dir.zenith, d.dir.azimuth, random_service=random_service, delta_psi=new_psi, is_degree=True)

            # Insert key to frame
            frame['I3MapSplit']['delta_psi_in_degree'] = new_psi
            frame['I3MapSplit']['new_zenith'] = new_zenith_azimuth[0][0]
            frame['I3MapSplit']['new_azimuth'] = new_zenith_azimuth[1][0]
            frame['I3MapSplit']['old_zenith'] = d.dir.zenith
            frame['I3MapSplit']['old_azimuth'] = d.dir.azimuth
            frame['I3MapSplit']['time'] = d.time
            
            # wrong! get lenght of track via daughters! check ditstance in detector
            # i3map['pre_length'] = (tree[0].pos - d.pos).magnitude
            # i3map['post_length'] = (d.pos - tree[-1].pos).magnitude

            # Insert pos key to frame
            frame.Put('I3PosOfMaxEnergyLoss', dataclasses.I3Position(d.pos))
            
            if brems_nuclint == False:
                break

        if found_max == True:
            
            ### FIX: CHANGE DIR ONLY OF MUON TRACK, NOT WHOLE TREE
            
            # Set new direction
            tree[i].dir = dataclasses.I3Direction(new_zenith_azimuth[0][0], new_zenith_azimuth[1][0])

            # Set new position
            tree[i].pos = muon_split_helper.pos + muon_split_helper.dir * (muon_split_helper.pos - d.pos).magnitude
        

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




##########################################################################





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
                print('Type: ', d_type)
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
            
            
def get_max_energy_loss_id_nugen_old(tree, muon_id, energy_loss, min_dist):
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


def get_max_energy_loss_(tree):
    
    distances = [distance_to_icecube_hull(d.pos) for d in tree]
    energies = [d.energy if d.type_string not in ('MuMinus', 'MuPlus') and distances[l] < 0 else 0 for d,l in zip(tree, range(len(distances)))]
    
    max_e_loss = max(energies)
    max_index = energies.index(max_e_loss)   
    # If max_e_loss = 0 --> muon was not going through the detector
    hull_distance_of_max_e_loss = distances[max_index] 
    return [max_index, max_e_loss, hull_distance_of_max_e_loss, E_before_max_E_loss]
    
    


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
    tree_len = len(tree)
    while ((tree[muon_id].type_string != 'MuMinus') & (tree[muon_id].type_string != 'MuPlus')):
        muon_id += 1
        if muon_id >= tree_len:
            return False
    
    muon = tree[muon_id]    
      
    # Check if particle moves through detector
    if particle_is_inside(muon, self._convex_hull) != True:
        # print('particle is not going through detector')
        return False

    # Save output of get_max_energy_loss_id (id, energy, negative distance to hull)
    id_E_dist = get_max_energy_loss(tree) # , muon_id, self._percentage_energy_loss, self._min_dist)

    # Check if muon is going through detector --> about 60% of the events are rejected
    if id_E_dist[1] == 0:
        return False
    
    # Check for minimum energy loss: 30% loss -> 3% left, 50% loss -> 1.3% left
    min_energy_loss = self._percentage_energy_loss * muon.energy
    if min_energy_loss > id_E_dist[1]:
        print('energy loss is too low')
        return False
    
    # Energy loss should be in the middle of the detector
    hull_distance_of_max_e_loss = id_E_dist[2] # this value is negative
    if self._min_dist < hull_distance_of_max_e_loss:
        print('energy loss not in the middle of detector')
        return False
    
    # Save data in tree 
    frame.Put('I3MapSplit', dataclasses.I3MapStringDouble())
    frame['I3MapSplit']['muon_id'] = muon_id
    frame['I3MapSplit']['max_E_id'] = id_E_dist[0]
    frame['I3MapSplit']['max_E_loss'] = id_E_dist[1]
    frame['I3MapSplit']['hull_dist'] = id_E_dist[2]
    frame['I3MapSplit']['E_before_max_E_loss'] = id_E_dist[3]