from icecube import dataclasses, icetray, dataio
from I3Tray import *
import numpy as np
from scipy.spatial import ConvexHull

# from egenerator.utils.angles import get_delta_psi_vector, get_delta_psi_vector_dir

# from resources.geometry import get_closest_point_on_edge, get_distance_to_edge, get_edge_intersection, distance_to_axis_aligned_Volume, distance_to_icecube_hull

from resources.muon_split_functions import build_tree_with_muon_split

class SplitModule(icetray.I3ConditionalModule):

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
        self.AddParameter('NewPsi', 
                        'angle in degree to set new direction after split')
        self.AddParameter('RandomSeed',
                        'int for random service',
                        42)

    def Configure(self):
        """Set new angle and random_seed
        """
        self.new_psi = self.GetParameter('NewPsi')
        self.random_seed = self.GetParameter('RandomSeed')

    def Geometry(self, frame):
        """Summary

        Parameters
        ----------
        frame : I3Frame
            Current geometry frame
        """
        print('---------Geometry method was used---------')

        #geoMap = frame['I3Geometry'].omgeo
        #domPosDict = {(i[0][0], i[0][1]): (i[1].position.x,
        #                                   i[1].position.y,
        #                                   i[1].position.z)
        #              for i in geoMap if i[1].omtype.name == 'IceCube'}
        #points = [
        #    domPosDict[(31, 1)], domPosDict[(1, 1)],
        #    domPosDict[(6, 1)], domPosDict[(50, 1)],
        #    domPosDict[(74, 1)], domPosDict[(72, 1)],
        #    domPosDict[(78, 1)], domPosDict[(75, 1)],

        #    domPosDict[(31, 60)], domPosDict[(1, 60)],
        #    domPosDict[(6, 60)], domPosDict[(50, 60)],
        #    domPosDict[(74, 60)], domPosDict[(72, 60)],
        #    domPosDict[(78, 60)], domPosDict[(75, 60)]
        #    ]
        #self._convex_hull = ConvexHull(points)
        #self._dom_pos_dict = domPosDict
        self.PushFrame(frame)

    def DAQ(self, frame):
        """Process DAQ frame

        Parameters
        ----------
        frame : I3Frame
            The current Q-frame.
        """

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
