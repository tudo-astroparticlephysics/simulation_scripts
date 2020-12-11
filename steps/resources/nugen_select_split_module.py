from icecube import dataclasses, icetray, dataio
from I3Tray import *
import numpy as np
from scipy.spatial import ConvexHull

from resources.nugen_muon_split_functions import insert_deflection_angle, selection, cut_millipede_out_of_detector

class NuGenSelectSplitModule(icetray.I3ConditionalModule):

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
                        'minimum distance of highest energy loss in detector',
                        -100)
        self.AddParameter('percentage_energy_loss',
                        'percentage energy loss',
                        0.3)
        self.AddParameter('NewPsi', 
                        'angle in degree to set new direction after split')
        self.AddParameter('RandomSeed',
                        'int for random service',
                        42)
        self.AddParameter('beta', 
                        'int to scale sampling angle',
                          1)
        self.AddParameter('perform_cut',
                         'bool to write down all events and only check selection',
                         False)

    def Configure(self):
        """Set 
        """
        self._min_dist = self.GetParameter('MinDist')
        self._percentage_energy_loss = self.GetParameter('percentage_energy_loss')
        self._new_psi = self.GetParameter('NewPsi')
        # self._random_seed = self.GetParameter('RandomSeed')
        self._random_service = np.random.RandomState(self.GetParameter('RandomSeed'))
        self._beta = self.GetParameter('beta')
        self._perform_cut = self.GetParameter('perform_cut')

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
        
        # if selection_nugen(self, frame) == False:
            # return False
        
        # if build_tree_with_muon_split_nugen(frame, self._new_psi, self._random_service, self._beta) == False:
            # return False

        self.PushFrame(frame)

    def Physics(self, frame):
        """Process Physics frame

        Parameters
        ----------
        frame : I3Frame
            The current P-Frame.
        """
        
        if selection(self, frame, self._perform_cut) == False:
            return False
        
        insert_deflection_angle(frame, self._new_psi, self._random_service, self._beta)
        
        cut_millipede_out_of_detector(frame)
        
        
        # if selection_nugen(self, frame) == False:
            # return False
        
        # if build_tree_with_muon_split_nugen(frame, self._new_psi, self._random_service, self._beta) == False:
            # return False
        
        self.PushFrame(frame)
