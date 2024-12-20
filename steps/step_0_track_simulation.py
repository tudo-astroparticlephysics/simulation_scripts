#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/icetray-start
#METAPROJECT icetray/v1.10.0
from __future__ import division
import os
import sys
if "ENV_SITE_PACKAGES" in os.environ:
    sys.path.insert(1, os.environ["ENV_SITE_PACKAGES"])

import time
import click
import yaml
import numpy as np

from icecube.simprod import segments

from icecube import icetray, dataclasses

from utils import create_random_services, get_run_folder
from resources.track_factory import SphereTrackFactory
from resources.oversampling import DAQFrameMultiplier


class DummyMCTreeRenaming(icetray.I3ConditionalModule):
    def __init__(self, context):
        """Class to add dummy I3MCTree to frame from I3MCTree_preMuonProp

        Parameters
        ----------
        context : TYPE
            Description
        """
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddOutBox('OutBox')

    def DAQ(self, frame):
        """Inject casacdes into I3MCtree.

        Parameters
        ----------
        frame : icetray.I3Frame.DAQ
            An I3 q-frame.
        """

        pre_tree = frame['I3MCTree_preMuonProp']
        frame['I3MCTree'] = dataclasses.I3MCTree(pre_tree)
        self.PushFrame(frame)


@click.command()
@click.argument('cfg', type=click.Path(exists=True))
@click.argument('run_number', type=int)
@click.option('--scratch/--no-scratch', default=True)
def main(cfg, run_number, scratch):
    start_time = time.time()

    with open(cfg, 'r') as stream:
        cfg = yaml.full_load(stream)
    cfg['run_number'] = run_number
    cfg['run_folder'] = get_run_folder(run_number)
    if scratch:
        outfile = cfg['scratchfile_pattern'].format(**cfg)
    else:
        outfile = cfg['outfile_pattern'].format(**cfg)
    outfile = outfile.replace(' ', '0')

    click.echo('Run: {}'.format(run_number))
    click.echo('Outfile: {}'.format(outfile))
    click.echo('Azimuth: [{},{}]'.format(*cfg['azimuth_range']))
    click.echo('Zenith: [{},{}]'.format(*cfg['zenith_range']))
    click.echo('Energy: [{},{}]'.format(*cfg['track_energy_range']))
    click.echo('Anchor x: [{},{}]'.format(*cfg['x_range']))
    click.echo('Anchor y: [{},{}]'.format(*cfg['y_range']))
    click.echo('Anchor z: [{},{}]'.format(*cfg['z_range']))
    click.echo('Max anchor dist: {}'.format(cfg['max_anchor_distance']))
    click.echo('Inj. radius: {}'.format(cfg['injection_sphere_radius']))
    click.echo('Sample uniformly: {}'.format(cfg['sample_uniformly_on_sphere']))

    # crate random services
    if 'random_service_use_gslrng' not in cfg:
        cfg['random_service_use_gslrng'] = False
    random_services, _ = create_random_services(
        dataset_number=cfg['dataset_number'],
        run_number=cfg['run_number'],
        seed=cfg['seed'],
        n_services=2,
        use_gslrng=cfg['random_service_use_gslrng'])

    # --------------------------------------
    # Build IceTray
    # --------------------------------------
    tray = icetray.I3Tray()
    tray.AddModule('I3InfiniteSource', 'source',
                   # Prefix=gcdfile,
                   Stream=icetray.I3Frame.DAQ)

    if 'track_injection_kwargs' not in cfg:
        cfg['track_injection_kwargs'] = {}
    if 'constant_vars' not in cfg:
        cfg['constant_vars'] = None
    if 'oversample_after_proposal' in cfg and \
            cfg['oversample_after_proposal']:
        oversampling_factor_injection = None
        oversampling_factor_photon = cfg['oversampling_factor']
    else:
        oversampling_factor_injection = cfg['oversampling_factor']
        oversampling_factor_photon = None

    tray.AddModule(SphereTrackFactory,
                   'make_tracks',
                   azimuth_range=cfg['azimuth_range'],
                   zenith_range=cfg['zenith_range'],
                   sample_uniformly_on_sphere=cfg[
                                        'sample_uniformly_on_sphere'],
                   track_energy_range=cfg['track_energy_range'],
                   time_range=cfg['time_range'],
                   x_range=cfg['x_range'],
                   y_range=cfg['y_range'],
                   z_range=cfg['z_range'],
                   max_anchor_distance=cfg['max_anchor_distance'],
                   injection_sphere_radius=cfg['injection_sphere_radius'],
                   num_events=cfg['n_events_per_run'],
                   oversampling_factor=oversampling_factor_injection,
                   random_state=cfg['seed'],
                   random_service=random_services[0],
                   constant_vars=cfg['constant_vars'],
                   **cfg['track_injection_kwargs']
                   )

    # propagate muons if config exists in config
    # Note: Snowstorm may perform muon propagation internally
    if 'muon_propagation_config' in cfg:
        tray.AddSegment(segments.PropagateMuons,
                        'propagate_muons',
                        RandomService=random_services[1],
                        **cfg['muon_propagation_config'])
    else:
        # In this case we are not propagating the I3MCTree yet, but
        # are letting this be done by snowstorm propagation
        # We need to add a key named 'I3MCTree', since snowstorm expects this
        # It will propagate the particles for us.
        tray.AddModule(DummyMCTreeRenaming, 'DummyMCTreeRenaming')

    tray.AddModule(DAQFrameMultiplier, 'DAQFrameMultiplier',
                   oversampling_factor=oversampling_factor_photon)

    # --------------------------------------
    # Distance Splits
    # --------------------------------------
    if cfg['distance_splits'] is not None:
        click.echo('SplittingDistance: {}'.format(
            cfg['distance_splits']))
        distance_splits = np.atleast_1d(cfg['distance_splits'])
        dom_limits = np.atleast_1d(cfg['threshold_doms'])
        if len(dom_limits) == 1:
            dom_limits = np.ones_like(distance_splits) * cfg['threshold_doms']
        oversize_factors = np.atleast_1d(cfg['oversize_factors'])
        order = np.argsort(distance_splits)

        distance_splits = distance_splits[order]
        dom_limits = dom_limits[order]
        oversize_factors = oversize_factors[order]

        stream_objects = generate_stream_object(distance_splits,
                                                dom_limits,
                                                oversize_factors)
        tray.AddModule(OversizeSplitterNSplits,
                       "OversizeSplitterNSplits",
                       thresholds=distance_splits,
                       thresholds_doms=dom_limits,
                       oversize_factors=oversize_factors)
        for stream_i in stream_objects:
            outfile_i = stream_i.transform_filepath(outfile)
            tray.AddModule("I3Writer",
                           "writer_{}".format(stream_i.stream_name),
                           Filename=outfile_i,
                           Streams=[icetray.I3Frame.DAQ,
                                    icetray.I3Frame.Physics,
                                    icetray.I3Frame.Stream('S'),
                                    icetray.I3Frame.Stream('M')],
                           If=stream_i)
            click.echo('Output ({}): {}'.format(stream_i.stream_name,
                                                outfile_i))
    else:
        click.echo('Output: {}'.format(outfile))
        tray.AddModule("I3Writer", "writer",
                       Filename=outfile,
                       Streams=[icetray.I3Frame.DAQ,
                                icetray.I3Frame.Physics,
                                icetray.I3Frame.Stream('S'),
                                icetray.I3Frame.Stream('M')])
    # --------------------------------------

    click.echo('Scratch: {}'.format(scratch))
    tray.AddModule("TrashCan", "the can")
    tray.Execute()
    tray.Finish()

    end_time = time.time()
    print("That took "+str(end_time - start_time)+" seconds.")


if __name__ == '__main__':
    main()
