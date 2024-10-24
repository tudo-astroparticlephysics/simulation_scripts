import numpy as np

from icecube import icetray, dataclasses
from icecube.icetray import I3Units

from . import geometry


class SphereTrackFactory(icetray.I3ConditionalModule):
    def __init__(self, context):
        """Class to create and inject Tracks.

        Tracks are created at a sampled anchor point within the specified
        bounds. The specified time and energy range is defined at the point
        of entry of the track in the sphere around the detector center.

        Steps to simualte:
            1. Sample anchor point within bounds
            2. Sample direction, time and energy
            3. Compute point of entry from anchor point
            4. Create track starting at the entry point in the sphere


        Parameters
        ----------
        context : TYPE
            Description
        """
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddOutBox('OutBox')
        self.AddParameter('azimuth_range',
                          '[min, max] of track azimuth in degree.', [0, 360])
        self.AddParameter('zenith_range',
                          '[min, max] of track zenith in degree.', [0, 180])
        self.AddParameter('sample_uniformly_on_sphere',
                          'If True, zenith is sampled uniformly in cos(zenith)'
                          ' which results in a uniform distribution on the '
                          'sphere. If False, zenith is sampled uniformly in '
                          'zenith which leads to more densly sampled points '
                          'at the poles of the sphere', True)
        self.AddParameter('track_energy_range', '', [10000, 10000])
        self.AddParameter('time_range', '[min, max] of entry time in ns.',
                          [9000, 12000])
        self.AddParameter('x_range',
                          '[min, max] of anchor x-coordinate in meters.',
                          [-500, 500])
        self.AddParameter('y_range',
                          '[min, max] of anchor y-coordinate in meters.',
                          [-500, 500])
        self.AddParameter('z_range',
                          '[min, max] of anchor z-coordinate in meters.',
                          [-500, 500])
        self.AddParameter('max_anchor_distance',
                          'Maximum distance of anchor outside of convex hull '
                          'around IceCube. If the drawn anchor is further '
                          'outside of the convex hull than the specified '
                          'amount, a new anchor position will be drawn.'
                          'If max_anchor_distance is None, the sampled anchor '
                          'position will be accepted regardless of its '
                          'distance to the convex hull.',
                          None)
        self.AddParameter('injection_sphere_radius',
                          'The radius of the sphere around the detector center '
                          'where the track is injected in meters.',
                          750)
        self.AddParameter('random_state', '', 1337)
        self.AddParameter('random_service', '', None)
        self.AddParameter('num_events', '', 1)
        self.AddParameter('oversampling_factor',
                          'Oversampling Factor to be used. Simulation is '
                          'averaged over these many simulations.',
                          None)
        self.AddParameter('constant_vars',
                          'These variables are only sampled once when the '
                          'module is being configured. They are kept constant '
                          'afterwards. This can for instance be used to keep '
                          'certain parameters such as the direction constant '
                          'for events of a run. Allowed options are: '
                          'anchor, zenith, azimuth, track_energy, time. '
                          'The variables must be passed as a list of str',
                          None)

    def Configure(self):
        """Configures Track Factory.

        Raises
        ------
        ValueError
            If interaction type or flavor is unkown.
        """
        self.azimuth_range = self.GetParameter('azimuth_range')
        self.zenith_range = self.GetParameter('zenith_range')
        self.sample_in_cos = self.GetParameter('sample_uniformly_on_sphere')
        self.cos_zenith_range = [np.cos(np.deg2rad(self.zenith_range[1])),
                                 np.cos(np.deg2rad(self.zenith_range[0]))]
        self.track_energy_range = self.GetParameter('track_energy_range')
        self.log_track_energy_range = [
                                    np.log10(self.track_energy_range[0]),
                                    np.log10(self.track_energy_range[1])]
        self.time_range = self.GetParameter('time_range')
        self.x_range = self.GetParameter('x_range')
        self.y_range = self.GetParameter('y_range')
        self.z_range = self.GetParameter('z_range')
        self.max_anchor_distance = self.GetParameter('max_anchor_distance')
        self.injection_sphere_radius = self.GetParameter('injection_sphere_radius')
        self.random_state = self.GetParameter('random_state')
        self.random_service = self.GetParameter('random_service')
        if not isinstance(self.random_state, np.random.RandomState):
            self.random_state = np.random.RandomState(self.random_state)
        self.num_events = self.GetParameter('num_events')
        self.oversampling_factor = self.GetParameter('oversampling_factor')
        if self.oversampling_factor is None:
            self.oversampling_factor = 1
        if self.max_anchor_distance is None:
            self.max_anchor_distance = float('inf')
        self.constant_vars = self.GetParameter('constant_vars')
        if self.constant_vars is None:
            self.constant_vars = []
        self.events_done = 0

        # make lowercase
        self.constant_vars = [f.lower() for f in self.constant_vars]

        # --------------
        # sanity checks:
        # --------------
        for const_var in self.constant_vars:
            if const_var not in ['anchor', 'zenith', 'azimuth', 'time',
                                 'track_energy']:
                raise ValueError('Var unknown: {!r}'.format(const_var))

        if self.oversampling_factor < 1:
            raise ValueError('Oversampling must be set to "None" or integer'
                             ' greater than 1. It is currently set to: '
                             '{!r}'.format(self.oversampling_factor))

        if self.injection_sphere_radius <= 600 + self.max_anchor_distance:
            raise ValueError(
                f'Injection sphere radius ({self.injection_sphere_radius})'
                'should be larger than 600m +  max_anchor_distance. '
                f'Currently set to: {self.injection_sphere_radius}'
            )

        # --------------------
        # sample constant vars
        # --------------------
        # anchor
        if 'anchor' in self.constant_vars:
            self.anchor = self._sample_anchor()

        if 'time' in self.constant_vars:
            self.entry_time = \
                self.random_service.uniform(*self.time_range)*I3Units.ns

        # direction
        if 'azimuth' in self.constant_vars:
            self.azimuth = \
                self.random_service.uniform(*self.azimuth_range)*I3Units.deg
        if 'zenith' in self.constant_vars:
            if self.sample_in_cos:
                zenith = np.rad2deg(np.arccos(
                    self.random_service.uniform(*self.cos_zenith_range)))
            else:
                zenith = self.random_service.uniform(*self.zenith_range)
            self.zenith = zenith*I3Units.deg

        # energy
        if 'track_energy' in self.constant_vars:
            self.log_track_energy = self.random_service.uniform(
                                *self.log_track_energy_range) * I3Units.GeV
        # --------------------

    def _sample_anchor(self):
        """Sample an anchor within allowd distance of IceCube Convex Hull.

        Returns
        -------
        I3Position
            The sampled anchor point.
        """
        # anchor
        point_is_inside = False
        while not point_is_inside:
            anchor_x = self.random_service.uniform(*self.x_range) * I3Units.m
            anchor_y = self.random_service.uniform(*self.y_range) * I3Units.m
            anchor_z = self.random_service.uniform(*self.z_range) * I3Units.m
            anchor = dataclasses.I3Position(
                            anchor_x * I3Units.m,
                            anchor_y * I3Units.m,
                            anchor_z * I3Units.m)
            dist = geometry.distance_to_icecube_hull(anchor)
            point_is_inside = (
                dist < self.max_anchor_distance and
                anchor.magnitude < self.injection_sphere_radius -1e-6
            )

        return anchor

    def _get_sphere_intersection(self, anchor, direction):
        """Get intersection distances of track with sphere around detector center.

        Parameters
        ----------
        anchor : I3Position
            Anchor point of track.
        direction : I3Direction
            Direction of track.

        Returns
        -------
        tuple
            The distance to the entry and exit point of the track with the
            sphere around the detector center.
        """
        if anchor.magnitude >= self.injection_sphere_radius - 1e-6:
            raise ValueError('Anchor point is outside of injection sphere.')

        p = (direction.x*anchor.x + direction.y*anchor.y + direction.z*anchor.z)
        sqrt_term = np.sqrt(p**2 - (anchor.magnitude**2 - self.injection_sphere_radius**2))

        dist_entry = -p - sqrt_term
        dist_exit = -p + sqrt_term
        return dist_entry, dist_exit

    def _generate_daughters(self, primary, track_energy, **kwargs):
        """Generate daughter particles.

        Parameters
        ----------
        primary : I3Particle
            The primary particle.
        track_energy : float
            The energy of the track.
        **kwargs
            Additional keyword arguments

        Returns
        -------
        list
            A list of daughter particles.
        """
        daughter = dataclasses.I3Particle()
        daughter.time = primary.time
        daughter.dir = primary.dir
        daughter.speed = primary.speed
        daughter.pos = primary.pos
        daughter.energy = track_energy * I3Units.GeV
        daughter.shape = dataclasses.I3Particle.InfiniteTrack
        daughter.type = dataclasses.I3Particle.ParticleType.MuMinus
        daughter.location_type = dataclasses.I3Particle.LocationType.InIce
        return [daughter]

    def DAQ(self, frame):
        """Inject track into I3MCtree.

        Parameters
        ----------
        frame : icetray.I3Frame.DAQ
            An I3 q-frame.

        Raises
        ------
        ValueError
            If interaction type is unknown.
        """
        # --------------
        # sample track
        # --------------
        # anchor
        if 'anchor' in self.constant_vars:
            anchor = self.anchor
        else:
            anchor = self._sample_anchor()

        # entry time
        if 'time' in self.constant_vars:
            entry_time = self.entrytime
        else:
            entry_time = \
                self.random_service.uniform(*self.time_range)*I3Units.ns

        # direction
        if 'azimuth' in self.constant_vars:
            azimuth = self.azimuth
        else:
            azimuth = \
                self.random_service.uniform(*self.azimuth_range)*I3Units.deg
        if 'zenith' in self.constant_vars:
            zenith = self.zenith
        else:
            if self.sample_in_cos:
                zenith = np.rad2deg(np.arccos(
                    self.random_service.uniform(*self.cos_zenith_range)))
            else:
                zenith = self.random_service.uniform(*self.zenith_range)
            zenith = zenith*I3Units.deg
        direction = dataclasses.I3Direction(zenith, azimuth)

        # energy
        if 'track_energy' in self.constant_vars:
            log_track_energy = self.log_track_energy
        else:
            log_track_energy = self.random_service.uniform(
                                *self.log_track_energy_range) * I3Units.GeV
        track_energy = 10**log_track_energy

        # compute entry point
        dist_entry, _ = self._get_sphere_intersection(anchor, direction)
        entry = anchor + dist_entry * direction

        assert np.allclose(entry.magnitude, self.injection_sphere_radius, atol=1e-3), entry

        # create particle
        primary = dataclasses.I3Particle()

        primary.time = entry_time * I3Units.ns
        primary.dir = direction
        primary.pos = entry
        primary.speed = dataclasses.I3Constants.c
        # Assume the entry position in range is in ice, so the primary is the
        # in ice neutrino that interacts
        primary.location_type = dataclasses.I3Particle.LocationType.InIce

        # create daughters
        daughters = self._generate_daughters(primary, track_energy)

        # oversampling
        for i in range(self.oversampling_factor):
            if i > 0:
                # create a new frame
                frame = icetray.I3Frame(frame)
                del frame['I3MCTree_preMuonProp']
                del frame['oversampling']

            # Fill primary and daughter particles into a MCTree
            primary_copy = dataclasses.I3Particle(primary)
            mctree = dataclasses.I3MCTree()
            mctree.add_primary(primary_copy)
            for daughter in daughters:
                mctree.append_child(primary_copy, dataclasses.I3Particle(daughter))

            frame['I3MCTree_preMuonProp'] = mctree
            if self.oversampling_factor > 1:
                frame['oversampling'] = dataclasses.I3MapStringInt({
                                        'event_num_in_run': self.events_done,
                                        'oversampling_num': i,
                                    })
            self.PushFrame(frame)

        self.events_done += 1
        if self.events_done >= self.num_events:
            self.RequestSuspension()


class SphereBundleFactory(SphereTrackFactory):
    def __init__(self, context):
        raise NotImplementedError

    def _generate_daughters(self, primary, track_energy, **kwargs):
        """Generate daughter particles.

        Parameters
        ----------
        primary : I3Particle
            The primary particle.
        track_energy : float
            The energy of the track.
        **kwargs
            Additional keyword arguments

        Returns
        -------
        list
            A list of daughter particles.
        """
        raise NotImplementedError