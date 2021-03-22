from icecube import dataclasses, icetray, dataio
from I3Tray import *
import numpy as np

from resources.geometry import distance_to_icecube_hull




# from ic3_labels.labels.utils.general import particle_is_inside 
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

    # Initialize tree with the new I3MCTree and set primary and muon
    muon = tree[1]

    max_index = int(frame['I3MapSplit']['max_E_id'])

    # Create variables to fill the tree
    found_max = False

    
    # Loop through particles in oldTree except the first two particles because they already exist in tree
    for i in range(len(tree)):
        d = tree[i]
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
    # distances = [0, 0]
    
    distances = [distance_to_icecube_hull(d.pos) for d in tree]
    energies = [d.energy if d.type_string != 'MuMinus' and distances[l] < 0 else 0 for d,l in zip(tree, range(len(distances)))]
    
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
    
    max_e_loss = max(energies)
    max_index = energies.index(max_e_loss)   
    # If max_e_loss = 0 --> muon was not going through the detector
    hull_distance_of_max_e_loss = distances[max_index] 
    return [max_index, max_e_loss, hull_distance_of_max_e_loss]


def build_tree_with_muon_split_oldversion(frame, new_psi, random_seed):
    '''Change direction of particles after maximum energy loss happens. First save the old tree, 
    the existing tree will 
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
    
    # Create variables to fill the tree
    found_max = False
    daughter_counter = 0
    daughter_direction_change = False
    daughter_daughter_counter = 0
    daughter_daughter_direction_change = False
    daughter_daughter_daughter_counter = 0
    daughter_daughter_daughter_direction_change = False

    # Loop through particles in oldTree except the first two particles because they already exist in tree
    for d in oldTree[2:]:

        # Check for maximum energy loss
        if d == oldTree[max_index]:
            # print('maximum energy loss: {}'.format(oldTree[max_index]))
            
            # Set parameter: maximum loss found, now muon2 has to be inserted, change direction of the following daughter particles
            found_max = True
            daughter_direction_change = True
            daughter_daughter_direction_change = True
            daughter_daughter_daughter_direction_change = True

            muon_split_helper = d


                       
            # Get new zenith and azimuth value for splitted particles with given delta_angle
            random_service = np.random.RandomState(random_seed)
            
            
            new_psi = np.exp(random_service.uniform(np.log(0.1), np.log(new_psi))) 
            new_zenith_azimuth = get_delta_psi_vector(muon.dir.zenith, muon.dir.azimuth, random_service=random_service, delta_psi=new_psi, is_degree=True)

            # Insert key to frame
            # frame.Put('I3MapSplit', dataclasses.I3MapStringDouble())
            i3map = frame['I3MapSplit']
            # Save opening angle and new direction
            i3map['delta_psi_in_degree'] = new_psi 
            i3map['new_zenith'] = new_zenith_azimuth[0][0]
            i3map['new_azimuth'] = new_zenith_azimuth[1][0]
            i3map['old_zenith'] = d.dir.zenith
            i3map['old_azimuth'] = d.dir.azimuth
            i3map['time'] = d.time
            i3map['pre_length'] = (oldTree[0].pos - oldTree[max_index].pos).magnitude
            i3map['post_length'] = (oldTree[max_index].pos - oldTree[-1].pos).magnitude

            # Insert pos key to frame
            frame.Put('I3PosOfMaxEnergyLoss', dataclasses.I3Position(d.pos))
                             
        # Inser daughter daughter daughter particles
        if daughter_daughter_daughter_counter != 0:
            if daughter_daughter_daughter_direction_change == True:
                d.dir = dataclasses.I3Direction(new_zenith_azimuth[0][0], new_zenith_azimuth[1][0])
                d.pos = muon_split_helper.pos + muon_split_helper.dir * (muon_split_helper.pos - d.pos).magnitude 
                tree.append_child(daughter_daughter_daughter_id, d)
            else:
                tree.append_child(daughter_daughter_daughter_id, d)
            daughter_daughter_daughter_counter -= 1
            continue

        # Insert daughter daughter particles
        if daughter_daughter_counter != 0:
            if daughter_daughter_direction_change == True:
                d.dir = dataclasses.I3Direction(new_zenith_azimuth[0][0], new_zenith_azimuth[1][0])
                d.pos = muon_split_helper.pos + muon_split_helper.dir * (muon_split_helper.pos - d.pos).magnitude 
                tree.append_child(daughter_daughter_id, d)
            else:
                tree.append_child(daughter_daughter_id, d)
            daughter_daughter_counter -= 1
            
            # Check for daughter daughter daughter
            daughter_daughter_daughter_counter = len(oldTree.get_daughters(d))
            if daughter_daughter_daughter_counter != 0:
                daughter_daughter_daughter_id = d.id
            continue                    
                          

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
            
            # Check for daughter daughter
            daughter_daughter_counter = len(oldTree.get_daughters(d))
            if daughter_daughter_counter != 0:
                daughter_daughter_id = d.id
                             
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
    if self._min_dist < hull_distance_of_max_e_loss:
        return False
    
    # Save data in tree 
    frame.Put('I3MapSplit', dataclasses.I3MapStringDouble())
    frame['I3MapSplit']['max_E_id'] = id_E_dist[0]
    frame['I3MapSplit']['max_E_loss'] = id_E_dist[1]
    frame['I3MapSplit']['hull_dist'] = id_E_dist[2]

def get_delta_psi_vector(zenith, azimuth, delta_psi,
                         random_service=None,
                         randomize_for_each_delta_psi=True,
                         is_degree=True,
                         return_angles=True):
    """Get new angles with an opening angle of delta_psi.
    Parameters
    ----------
    zenith : array_like
        The zenith angle of the input vector for which to compute a random
        new vector with an opening angle of delta_psi.
    azimuth : TYPE
        The azimuth angle of the input vector for which to compute a random
        new vector with an opening angle of delta_psi.
    delta_psi : float or array_like
        The opening angle. If 'is_degree' is True, then the unit is in degree,
        otherwise it is in radians.
        If an array is provided, broadcasting will be applied.
    random_service : None, optional
        An optional random number service to use for reproducibility.
    randomize_for_each_delta_psi : bool, optional
        If True, a random orthogonal vector is sampled for each specified
        delta_psi.
        If False, the direction vectors for the delta_psi opening angles
        are computed along along the same (random) geodesic.
    is_degree : bool, optional
        This specifies the input unit of 'delta_psi'.
        If True, the input unit of 'delta_psi' is degree.
        If False, it is radians.
    return_angles : bool, optional
        If True, the new random vector will be returned as zenith and azimuth
        angles: shape:  tuple([..., 1], [..., 1]).
        If False, it will be returned as a direction vector: shape: [..., 3].
    Returns
    -------
    array_like or tuple of array_like
        If return_angles is True:
            Return values are (zenith, azimuth)
            Shape:  tuple([..., 1], [..., 1])
        If return_angles is False:
            Return values are the new direction vectors in cartesian
            coordinates.
            Shape: [..., 3]
    """
    vec = np.array([np.sin(zenith) * np.cos(azimuth),
                    np.sin(zenith) * np.sin(azimuth),
                    np.cos(zenith)]).T
    vec = np.atleast_2d(vec)
    delta_vec = get_delta_psi_vector_dir(
        vec,
        delta_psi=delta_psi,
        random_service=random_service,
        randomize_for_each_delta_psi=randomize_for_each_delta_psi,
        is_degree=is_degree)
    if return_angles:
        # calculate zenith
        d_zenith = np.arccos(np.clip(delta_vec[..., 2], -1, 1))

        # calculate azimuth
        d_azimuth = (np.arctan2(delta_vec[..., 1], delta_vec[..., 0])
                     + 2 * np.pi) % (2 * np.pi)
        return d_zenith, d_azimuth
    else:
        return delta_vec


def get_delta_psi_vector_dir(vec, delta_psi,
                             randomize_for_each_delta_psi=True,
                             random_service=None,
                             is_degree=True):
    """Get a new direction vector with an opening angle of delta_psi to vec.
    Parameters
    ----------
    vec : array_like
        The vector for which to calculate a new random vector with an opening
        angle of delta_psi.
        Shape: [..., 3]
    delta_psi : float or array_like
        The opening angle. If 'is_degree' is True, then the unit is in degree,
        otherwise it is in radians.
        If an array is provided, broadcasting will be applied.
    randomize_for_each_delta_psi : bool, optional
        If True, a random orthogonal vector is sampled for each specified
        delta_psi.
        If False, the direction vectors for the delta_psi opening angles
        are computed along along the same (random) geodesic.
    random_service : None, optional
        An optional random number service to use for reproducibility.
    is_degree : bool, optional
        This specifies the input unit of 'delta_psi'.
        If True, the input unit of 'delta_psi' is degree.
        If False, it is radians.
    Returns
    -------
    array_like
        The new random vectors with an opening angle of 'delta_psi'.
        Shape: [..., 3]
    Raises
    ------
    ValueError
        If the specified opening angle is larger or equal to 90 degree.
        This calculation only supports angles up to 90 degree.
    """
    if random_service is None:
        random_service = np.random

    if is_degree:
        delta_psi = np.deg2rad(delta_psi)

    # allow broadcasting
    delta_psi = np.expand_dims(delta_psi, axis=-1)

    # This calculation is only valid if delta_psi < 90 degree
    if np.any(delta_psi > np.deg2rad(90)):
        msg = 'Delta Psi angle must be smaller than 90 degrees, but it is {!r}'
        raise ValueError(msg.format(np.rad2deg(delta_psi)))

    # get a random orthogonal vector
    if randomize_for_each_delta_psi:
        num_temp_vecs = max(len(vec), len(delta_psi))
    else:
        num_temp_vecs = len(vec)

    temp_vec = random_service.uniform(low=-1, high=1, size=(num_temp_vecs, 3))

    vec_orthogonal = np.cross(vec, temp_vec)
    vec_orthogonal /= np.linalg.norm(vec_orthogonal, axis=-1, keepdims=True)

    # calculate new vector with specified opening angle
    new_vec = vec + np.tan(delta_psi) * vec_orthogonal
    new_vec /= np.linalg.norm(new_vec, axis=-1, keepdims=True)
    return new_vec


def particle_is_inside(particle, convex_hull):
    '''Checks if a particle is inside the convex hull.
    The particle is considered inside if any part of its track is inside
    the convex hull. In the case of point like particles with length zero,
    the particle will be considered to be inside if the vertex is inside
    the convex hull.

    Parameters
    ----------
    particle : I3Particle
        The Particle to check.
    convex_hull : scipy.spatial.ConvexHull
        Defines the desired convex volume.

    Returns
    -------
    bool
        True if particle is inside, otherwise False.
    '''
    v_pos = (particle.pos.x, particle.pos.y, particle.pos.z)
    v_dir = (particle.dir.x, particle.dir.y, particle.dir.z)
    intersection_ts = geometry.get_intersections(convex_hull, v_pos, v_dir)

    # particle didn't hit convex_hull
    if intersection_ts.size == 0:
        return None

    # particle hit convex_hull:
    #   Expecting two intersections
    #   What happens if track is exactly along edge of hull?
    #   If only one ts: track exactly hit a corner of hull?
    # assert len(intersection_ts) == 2, 'Expected exactly 2 intersections'
    if len(intersection_ts) != 2:
        print('intersection_ts != 2')
        return False
    
    
    min_ts = min(intersection_ts)
    max_ts = max(intersection_ts)
    if min_ts <= 0 and max_ts >= 0:
        # starting event
        return True
    if max_ts < 0:
        # particle created after the convex hull
        return False
    if min_ts > particle.length + 1e-8:
        # particle stops before convex hull
        return False
    # everything else
    return True
