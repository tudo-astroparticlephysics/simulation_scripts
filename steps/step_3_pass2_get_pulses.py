#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT icerec/V05-01-06
import os
import sys
import time

import click
import yaml

from I3Tray import I3Tray, I3Units
from icecube import icetray, dataclasses
from icecube.filterscripts import filter_globals
from icecube.filterscripts.baseproc import BaseProcessing
from icecube.STTools.seededRT.configuration_services import \
    I3DOMLinkSeededRTConfigurationService
from icecube import filter_tools

from utils import get_run_folder

from resources.pulses import GetPulses
from resources.pulses import MergeOversampledEvents


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

    # get pulses
    tray.AddSegment(GetPulses, "GetPulses",
                    decode=False,
                    simulation=True,
                    )

    # Throw out unneeded streams and keys
    if 'oversampling_keep_keys' not in cfg:
        cfg['oversampling_keep_keys'] = []
    elif cfg['oversampling_keep_keys'] is None:
        cfg['oversampling_keep_keys'] = []

    if cfg['L1_keep_untriggered']:
        stream_name = filter_globals.InIceSplitter
    else:
        stream_name = filter_globals.NullSplitter
    tray.AddModule("KeepFromSubstream", "DeleteSubstream",
                   StreamName=stream_name,
                   KeepKeys=['do_not_keep_anything'])

    # merge oversampled events: calculate average hits
    if cfg['oversampling_factor'] is not None and do_merging_if_necessary:
        if 'oversampling_merge_events' in cfg:
            merge_events = cfg['oversampling_merge_events']
        else:
            # backward compability
            merge_events = True

        if merge_events:
            tray.AddModule(MergeOversampledEvents, 'MergeOversampledEvents',
                           OversamplingFactor=cfg['oversampling_factor'])
    keys_to_keep = [
        'TimeShift',
        'I3MCTree_preMuonProp',
        'I3MCTree',
        'I3MCWeightDict',
        'CorsikaWeightMap',
        'MMCTrackList',
        'I3EventHeader',
        'I3SuperDST',
        'RNGState',
        'oversampling',
        'AggregatedPulses',
        'InIceDSTPulses',
        'InIceDSTPulsesTimeRange',
        'CalibrationErrata',
        'SaturationWindows',
        'DOMDeadTimesMC',
        'SplitUncleanedInIcePulses',
        'SplitUncleanedInIcePulsesTimeRange',
        'SplitUncleanedInIceDSTPulsesTimeRange',
        'I3TriggerHierarchy',
        'GCFilter_GCFilterMJD',
        ]
    keys_to_keep += filter_globals.inice_split_keeps + \
        filter_globals.onlinel2filter_keeps

    if 'event_import_settings' in cfg:
        import_keys = [cfg['event_import_settings']['rename_dict'].get(k, k)
                       for k in cfg['event_import_settings']['keys_to_import']]
        keys_to_keep += import_keys

    tray.AddModule("Keep", "keep_before_merge",
                   keys=keys_to_keep + cfg['oversampling_keep_keys'])

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
