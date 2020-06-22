from icecube import dataclasses, icetray, dataio
from I3Tray import *
import numpy as np

from egenerator.utils.angles import get_delta_psi_vector, get_delta_psi_vector_dir

from resources.geometry import get_closest_point_on_edge, get_distance_to_edge, get_edge_intersection, distance_to_axis_aligned_Volume, distance_to_icecube_hull

from ic3_labels.labels.utils.general import particle_is_inside 

from scipy.spatial import ConvexHull


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

def selection(self, frame):

    tree = frame['I3MCTree']
    if particle_is_inside(tree[1], self._convex_hull) != True:
        return False

    # Save output of get_max_energy_loss_id (id, energy, negative distance to hull)
    id_E_dist = get_max_energy_loss_id(tree)

    # Check if muon is going through detector --> about 60% of the events are rejected
    if id_E_dist[1] == 0:
        return False
    
    # Check for minimum energy loss: 30% loss -> 3% left, 50% loss -> 1.3% left
    min_energy_loss = self._E_loss * tree[1].energy
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
    
    return True

class SelectionModule(icetray.I3ConditionalModule):

    """Class to split muon track

    """

    def __init__(self, context):
        """Initializes class

        Parameters
        ----------
        context : TYPE
            Description
        """
        icetray.I3ConditionalModule.__init__(self, context)
        # self.AddOutBox('OutBox') ?? whats this ??
        self.AddParameter('MinDist', 
                        'minimum distance of highest energy loss in detector', -100)
        self.AddParameter('percentage_energy_loss',
                        'percentage energy loss',
                        0.3)

    def Configure(self):
        """Set 
        """
        self._min_dist = self.GetParameter('MinDist')
        self._E_loss = self.GetParameter('percentage_energy_loss')

    def Geometry(self, frame):
        """Summary

        Parameters
        ----------
        frame : I3Frame
            Current geometry frame
        """
        print('---------Geometry method was used---------')

        geoMap = frame['I3Geometry'].omgeo
        domPosDict = {(i[0][0], i[0][1]): (i[1].position.x,
                                           i[1].position.y,
                                           i[1].position.z)
                      for i in geoMap if i[1].omtype.name == 'IceCube'}
        points = [
            domPosDict[(31, 1)], domPosDict[(1, 1)],
            domPosDict[(6, 1)], domPosDict[(50, 1)],
            domPosDict[(74, 1)], domPosDict[(72, 1)],
            domPosDict[(78, 1)], domPosDict[(75, 1)],

            domPosDict[(31, 60)], domPosDict[(1, 60)],
            domPosDict[(6, 60)], domPosDict[(50, 60)],
            domPosDict[(74, 60)], domPosDict[(72, 60)],
            domPosDict[(78, 60)], domPosDict[(75, 60)]
            ]
        self._convex_hull = ConvexHull(points)
        self._dom_pos_dict = domPosDict
        self.PushFrame(frame)

    def DAQ(self, frame):
        """Process DAQ frame

        Parameters
        ----------
        frame : I3Frame
            The current Q-frame.
        """
        
        if selection(self, frame) == False:
            return False
        
        build_tree_with_muon_split(frame, self.new_psi, self.random_seed)

        self.PushFrame(frame)

    def Physics(self, frame):
        """Process Physics frame

        Parameters
        ----------
        frame : I3Frame
            The current P-Frame.
        """
        self.PushFrame(frame)
