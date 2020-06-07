from icecube import dataclasses, icetray, dataio
from I3Tray import *
import numpy as np


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
    if np.any(delta_psi >= np.deg2rad(90)):
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


### Function to get index of particle with the maximum energy loss ###

def get_max_energy_loss_id(tree):
    # Create energy array with two energies for the first two particles (unkown and muon)
    energies = [0, 0]
    
    for d in tree[2:]:
        print(d.minor_id, distance_to_icecube_hull(d.pos))
        # Stashes the muons at every point and searches for the maximum energy inside the detector
        if d.type_string != 'MuMinus' and distance_to_icecube_hull(d.pos) < 0:
            energies.append(d.energy)
        else:
            energies.append(0)
    
    max_index = energies.index(max(energies))    
    return max_index

### Function to determine if particle is in ice ###

def get_closest_point_on_edge(edge_point1, edge_point2, point):
    '''Function to determine the closest point
            on an edge defined by the two points
            edge_point1 and edge_point2
    Parameters
    ----------
    edge_point1 : array-like shape=(3,)
        First edge point .
    edge_point2 : array-like shape=(3,)
        Second edge point .
    point : array-like shape=(3,)
        point of which to find the distance
        to the edge
    Returns
    -------
    distance: array-like shape=(3,)
        closest point on the edge
    '''
    if edge_point1 == edge_point2:
        return ValueError('Points do not define line.')
    A = np.array(edge_point1)
    B = np.array(edge_point2)
    P = np.array(point)
    vec_edge = B - A
    vec_point = P - A
    norm_edge = np.linalg.norm(vec_edge)
    t_projection = np.dot(vec_edge, vec_point) / (norm_edge**2)

    t_clipped = min(1, max(t_projection, 0))
    closest_point = A + t_clipped*vec_edge

    return closest_point


def get_distance_to_edge(edge_point1, edge_point2, point):
    '''Function to determine the closest distance of a point
            to an edge defined by the two points
            edge_point1 and edge_point2
    Parameters
    ----------
    edge_point1 : array-like shape=(3,)
        First edge point .
    edge_point2 : array-like shape=(3,)
        Second edge point .
    point : array-like shape=(3,)
        point of which to find the distance
        to the edge
    Returns
    -------
    distance: float
    '''
    closest_point = get_closest_point_on_edge(edge_point1,
                                              edge_point2, point)
    distance = np.linalg.norm(closest_point - point)
    return distance

def get_edge_intersection(edge_point1, edge_point2, point):
    '''Returns t:
        edge_point1 + u*(edge_point2-edge_point1)
        =
        point + t * (0, 1, 0)
        if u is within [0,1].
        [Helper Function to find out if point
         is inside the icecube 2D Polygon]
    Parameters
    ----------
    edge_point1 : array-like shape=(3,)
        First edge point .
    edge_point2 : array-like shape=(3,)
        Second edge point .
    point : array-like shape=(3,)
        point of which to find the distance
        to the edge
    Returns
    -------
    t: float.
        If intersection is within edge
        othwise returns nan.
    '''
    if edge_point1 == edge_point2:
        return ValueError('Points do not define line.')
    A = np.array(edge_point1)
    B = np.array(edge_point2)
    P = np.array(point)
    vec_edge = B - A
    vec_point = P - A

    u = vec_point[0] / vec_edge[0]
    t = u * vec_edge[1] - vec_point[1]

    if u > -1e-8 and u < 1 + 1e-8:
        return t
    return float('nan')

def distance_to_axis_aligned_Volume(pos, points, z_min, z_max):
    '''Function to determine the closest distance of a point
       to the edge of a Volume defined by z_zmin,z_max and a
       2D-Polygon described through a List of counterclockwise
       points.
    Parameters
    ----------
    pos :I3Position
        Position.
    points : array-like shape=(?,3)
        List of counterclockwise points
        describing the polygon of the volume
        in the x-y-plane
    z_max : float
        Top layer of IceCube-Doms
    z_min : float
        Bottom layer of IceCube-Doms
    Returns
    -------
    distance: float
        closest distance from the point
        to the edge of the volume
        negativ if point is inside,
        positiv if point is outside
    '''
    no_of_points = len(points)
    edges = [(points[i], points[(i + 1) % (no_of_points)])
             for i in range(no_of_points)]
    xy_distance = float('inf')
    list_of_ts = []

    for edge in edges:
        x = (edge[0][0], edge[1][0])
        y = (edge[0][1], edge[1][1])
        distance = get_distance_to_edge(edge[0], edge[1],
                                        [pos[0], pos[1], 0])
        t = get_edge_intersection(edge[0], edge[1],
                                  [pos[0], pos[1], 0])
        if not np.isnan(t):
            list_of_ts.append(t)
        if distance < xy_distance:
            xy_distance = distance
    is_inside_xy = False
    if len(list_of_ts) == 2:
        # u's are pos and negativ
        if list_of_ts[0]*list_of_ts[1] < 0:
            is_inside_xy = True
        # point is exactly on border
        elif len([t for t in list_of_ts if t == 0]) == 1:
            is_inside_xy = True

    # ---- Calculate z_distance
    is_inside_z = False
    if pos[2] < z_min:
        # underneath detector
        z_distance = z_min - pos[2]
    elif pos[2] < z_max:
        # same height
        is_inside_z = True
        z_distance = min(pos[2] - z_min, z_max - pos[2])
    else:
        # above detector
        z_distance = pos[2] - z_max

    # ---- Combine distances
    if is_inside_z:
        if is_inside_xy:
            # inside detector
            distance = - min(xy_distance, z_distance)
        else:
            distance = xy_distance
    else:
        if is_inside_xy:
            distance = z_distance
        else:
            distance = np.sqrt(z_distance**2 + xy_distance**2)

    return distance

def distance_to_icecube_hull(pos, z_min=-502, z_max=501):
    '''Function to determine the closest distance of a point
            to the icecube hull. This is only
            an approximate distance.
    Parameters
    ----------
    pos :I3Position
        Position.
    z_max : float
        Top layer of IceCube-Doms
    z_min : float
        Bottom layer of IceCube-Doms
    Returns
    -------
    distance: float
        closest distance from the point
        to the icecube hull
        negativ if point is inside,
        positiv if point is outside
    '''
    points = [
           [-570.90002441, -125.13999939, 0],  # string 31
           [-256.14001465, -521.08001709, 0],  # string 1
           [ 361.        , -422.82998657, 0],  # string 6
           [ 576.36999512,  170.91999817, 0],  # string 50
           [ 338.44000244,  463.72000122, 0],  # string 74
           [ 101.04000092,  412.79000854, 0],  # string 72
           [  22.11000061,  509.5       , 0],  # string 78
           [-347.88000488,  451.51998901, 0],  # string 75
            ]
    return distance_to_axis_aligned_Volume(pos, points, z_min, z_max)




# Function to split muon track at highest energy loss and insert this particle in a the tree, but first save the old tree, so the existing I3MCTree will be replaced with a new tree      
def build_tree_with_muon_split(frame, new_psi, random_seed):
    # new_psi: angle in degree to set new direction of muon2 and the following particles
    # random_seed: int to set random service
    tree = frame['I3MCTree']

    # Create a tree and initialize it with the old I3MCTree 
    frame.Put('OldTree', dataclasses.I3MCTree(tree))
    oldTree = frame['OldTree']

    # Set muon
    muon = tree[1]

    # Empty new tree except unknown and muon
    tree.erase_children(muon)

    # Get index of particle with highest energy loss
    max_index = get_max_energy_loss_id(oldTree)
    
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
            # insert_muon2 = True
            daughter_direction_change = True

            muon_split_helper = d

                       
            # Get new zenith and azimuth value for splitted particles with given delta_angle
            random_service = np.random.RandomState(random_seed)
            new_zenith_azimuth = get_delta_psi_vector(muon.dir.zenith, muon.dir.azimuth, random_service=random_service, delta_psi=new_psi, is_degree=True)

            # Insert key to frame
            frame.Put('I3MapMuonSplit', dataclasses.I3MapStringDouble())
            i3map = frame['I3MapMuonSplit']
            # Save opening angle and new direction
            i3map['delta_psi_in_degree'] = new_psi 
            i3map['new_dir_x'] = new_zenith_azimuth[0][0]
            i3map['new_dir_y'] = new_zenith_azimuth[1][0]

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


    print('OldTree: {}'.format(oldTree))
    print('I3MCTree: {}'.format(tree))
    print('----------- next one -------------')
