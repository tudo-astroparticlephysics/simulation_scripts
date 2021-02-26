from __future__ import division
import numpy as np

from I3Tray import I3Tray
from icecube import icetray, dataclasses

from ic3_labels.labels.utils import muon as mu_utils
from ic3_labels.labels.utils import detector


class InjectSingleVetoMuon(icetray.I3ConditionalModule):

    """Class to inject an accompanying muon for a provided neutrino event

    This is intended to run after neutrino injection. It is expected that the
    first primary in the provided I3MCTree is the injected neutrino event.
    For this particle, an accompanying muon will be injected at the convex hull
    of the detector.

    Attributes
    ----------
    mctree_name : str
        The name of the I3MCTree key.
    """

    def __init__(self, context):
        """Class to inject accompanying muons

        Parameters
        ----------
        context : TYPE
            Description
        """
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddOutBox('OutBox')
        self.AddParameter(
            'mctree_name', 'The name of the I3MCTree to use.', 'I3MCTree')
        self.AddParameter(
            'n_frames_per_neutrino',
            'The number of coincidence events to inject per imported event.',
            1)
        self.AddParameter(
            'sampling_settings',
            'Settings specifying how the muon energy is sampled.',
            {'method': 'power_law', 'range': [10, 1e7], 'gamma': 2})
        self.AddParameter(
            'random_service',
            'The random service or seed to use. If this is an '
            'integer, a numpy random state will be created with '
            'the seed set to `random_service`',
            42)
        self.AddParameter(
            'output_key',
            'The key to which the sampling information is written to.',
            'MCVetoMuonInjectionInfo')

    def Configure(self):
        """Configures MuonLossProfileFilter.
        """
        self.mctree_name = self.GetParameter('n_frames_per_neutrino')
        self.n_frames_per_neutrino = self.GetParameter('n_frames_per_neutrino')
        self.sampling_settings = self.GetParameter('sampling_settings')
        self.random_service = self.GetParameter('random_service')
        self.sampling_method = self.sampling_settings['method']
        self.output_key = self.GetParameter('output_key')

        if isinstance(self.random_service, int):
            self.random_service = np.random.RandomState(self.random_service)

    def DAQ(self, frame):
        """Inject accompanying muons
        """

        for i in range(self.n_frames_per_neutrino):

            # get I3MCTree
            mc_tree = dataclasses.I3MCTree(frame[self.mctree_name])

            # primary particle
            primary = mc_tree.get_primaries()[0]

            # compute entry point
            entry, time, nu_energy = mu_utils.get_muon_entry_info(
                frame, primary, convex_hull=detector.icecube_hull)

            # in order to not land exactly on the conved hull and potentially
            # cause issues for label generation, we will walk back a little further
            # for primary and muon injection
            dt = time - 10 - primary.time
            inj_pos = primary.pos + dt * primary.dir
            inj_time = time - 10
            inj_dir = dataclasses.I3Direction(primary.dir)

            # inject new muon
            mc_tree, injection_info = self._inject_muon(
                mc_tree=mc_tree,
                inj_pos=inj_pos,
                inj_time=inj_time,
                inj_dir=inj_dir,
            )

            # copy frame
            frame_copy = icetray.I3Frame(frame)

            # replace I3MCTree
            del frame_copy[self.mctree_name]
            frame_copy[self.mctree_name] = mc_tree

            # add info to frame
            injection_info['injection_counter'] = float(i)
            frame[self.output_key] = injection_info

            # push frame on to subsequent modules
            self.PushFrame()

    def _powerlaw_sampler(e_min, e_max, gamma, num=1):
        """Sample from Powerlaw Distribution

        Sample `num` events from a power law with index gamma between x
        low and xhigh by using the analytic inversion method.
        The power law pdf is given by
        .. math::
           \mathrm{pdf}(\gamma) = x^{-\gamma} / \mathrm{norm}
        where norm ensures an area under curve of one. Positive spectral index
        gamma means a falling spectrum.
        Note: When :math:`\gamma=1` the integral is
        .. math::
           \int 1/x \mathrm{d}x = ln(x) + c
        This case is also handled.

        Parameters
        ----------
        e_min : float
            The lower bound of the PDF, needed for proper normalization.
        e_max : float
            The upper bound of the PDF, needed for proper normalization.
        gamma : float, optional
            Power law index.
        num : int
            Number of random numbers to generate.

        Returns
        -------
        np.ndarray
            The random numbers sampled from a powerlaw.
        """
        u = self.random_service.uniform(size=int(num))

        if gamma == 1:
            return np.exp(u * np.log(xhigh / xlow)) * xlow
        else:
            radicant = (u * (xhigh**(1. - gamma) - xlow**(1. - gamma))
                        + xlow**(1. - gamma))
            return radicant**(1. / (1. - gamma))

    def _sample_energy(self):
        """Sample Energy of the injected Muon

        Returns
        -------
        double
            The sampled energy
        """
        if self.sampling_method == 'power_law':
            energy = self._powerlaw_sampler(
                e_min=self.sampling_settings['range'][0],
                e_max=self.sampling_settings['range'][1],
                gamma=self.sampling_settings['gamma'],
            )

        else:
            raise ValueError('Unknown smapling method: {}'.format(
                self.sampling_method))
        return energy

    def _inject_muon(self, mc_tree, inj_pos, inj_time, inj_dir):
        """Inject accompanying muon in provided I3MCTree

        Parameters
        ----------
        mc_tree : I3MCTree
            The I3MCTree into which the muon will be injected.
        inj_pos : I3Position
            The position at which to inject the muon.
        inj_time : double
            The time at which to inject the muon.
        inj_dir : I3Direction
            The direction of the injected muon.

        Returns
        -------
        I3MCTree
            The modified I3MCTree with the injected Muon.
        I3MapStringDouble
            Information on the injected muon.
        """
        muon_primary = dataclasses.I3Particle()
        muon_primary.shape = dataclasses.I3Particle.ParticleShape.Primary
        muon_primary.dir = dataclasses.I3Direction(inj_dir)
        muon_primary.pos = dataclasses.I3Position(inj_pos)
        muon_primary.time = inj_time

        muon = dataclasses.I3Particle()
        muon.dir = dataclasses.I3Direction(inj_dir)
        muon.pos = dataclasses.I3Position(inj_pos)
        muon.time = inj_time
        muon.location_type = dataclasses.I3Particle.LocationType.InIce

        # sample type: MuPlus or MuMinus
        pdg_encoding = self.random_service.choice([13, -13])
        muon.pdg_encoding = pdg_encoding

        # sample energy
        muon.energy = self._sample_energy()

        # add muon primary to I3MCTree
        mc_tree.add_primary(muon_primary)

        # add muon as new child
        mc_tree.append_child(muon_primary, muon)

        # add info
        injection_info = dataclasses.I3MapStringDouble({
            'muon_energy': muon.energy,
            'muon_pdg_encoding': muon.pdg_encoding,
        })

        return mc_tree, injection_info
