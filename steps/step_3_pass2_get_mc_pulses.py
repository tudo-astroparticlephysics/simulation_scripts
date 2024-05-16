#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT icerec/V05-01-06
import os
import sys
import time

import click
import yaml

from I3Tray import I3Tray, I3Units
from icecube import icetray, dataclasses, simclasses
from icecube.filterscripts import filter_globals
from icecube.filterscripts.baseproc import BaseProcessing
from icecube.STTools.seededRT.configuration_services import \
    I3DOMLinkSeededRTConfigurationService
from icecube import filter_tools

from utils import get_run_folder
from step_3_pass2_get_pulses import MergeOversampledEvents

from resources.pulses import GetMCPulses


@click.command()
@click.argument('cfg', type=click.Path(exists=True))
@click.argument('run_number', type=int)
@click.option('--scratch/--no-scratch', default=True)
@click.argument('do_merging_if_necessary', type=bool, default=True)
def main(cfg, run_number, scratch, do_merging_if_necessary):
    with open(cfg, 'r') as stream:
        if int(yaml.__version__[0]) < 5:
            # backwards compatibility for yaml versions before version 5
            cfg = yaml.load(stream)
        else:
            cfg = yaml.full_load(stream)
    cfg['run_number'] = run_number
    cfg['run_folder'] = get_run_folder(run_number)

    infile = cfg['infile_pattern'].format(**cfg)
    infile = infile.replace(' ', '0')
    infile = infile.replace('Level0.{}'.format(cfg['previous_step']),
                            'Level0.{}'.format(cfg['previous_step'] % 10))
    infile = infile.replace('2012_pass2', 'pass2')

    if scratch:
        outfile = cfg['scratchfile_pattern'].format(**cfg)
    else:
        outfile = cfg['outfile_pattern'].format(**cfg)
    outfile = outfile.replace('Level0.{}'.format(cfg['step']),
                              'Level0.{}'.format(cfg['step'] % 10))
    outfile = outfile.replace(' ', '0')
    outfile = outfile.replace('2012_pass2', 'pass2')
    print('Outfile != $FINAL_OUT clean up for crashed scripts not possible!')

    tray = I3Tray()
    tray.AddModule('I3Reader',
                   'i3 reader',
                   FilenameList=[cfg['gcd_pass2'], infile])

    # get MC pulses
    tray.AddModule(GetMCPulses, "GetMCPulses",
                   I3MCPEorPulseSeriesMap='I3MCPESeriesMapWithoutNoise',
                   OutputKey='MCPulses',
                   CreatePFrames=True)

    # merge oversampled events: calculate average hits
    if cfg['oversampling_factor'] is not None and do_merging_if_necessary:
        if 'oversampling_merge_events' in cfg:
            merge_events = cfg['oversampling_merge_events']
        else:
            # backward compability
            merge_events = True

        if merge_events:
            tray.AddModule(MergeOversampledEvents, 'MergeOversampledEvents',
                           OversamplingFactor=cfg['oversampling_factor'],
                           PulseKey='MCPulses')

    # Make space and delete uneeded keys
    keys_to_delete = [
        'I3MCPESeriesMap',
        'I3MCPulseSeriesMap',
        'I3MCPESeriesMapWithoutNoise',
        'I3MCPulseSeriesMapParticleIDMap',
        'I3MCPulseSeriesMapPrimaryIDMap',
        'InIceRawData',
        'IceTopRawData',
        ]
    tray.AddModule('Delete', 'DeleteKeys',
                   keys=keys_to_delete)

    tray.AddModule("I3Writer", "EventWriter",
                   filename=outfile,
                   Streams=[icetray.I3Frame.DAQ,
                            icetray.I3Frame.Physics,
                            icetray.I3Frame.TrayInfo,
                            icetray.I3Frame.Simulation,
                            icetray.I3Frame.Stream('M')])
    tray.AddModule("TrashCan", "the can")
    tray.Execute()
    tray.Finish()


if __name__ == '__main__':
    main()
