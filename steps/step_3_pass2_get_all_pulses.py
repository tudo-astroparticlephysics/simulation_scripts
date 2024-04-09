#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/icetray-start
#METAPROJECT icetray/v1.10.0
import time
import click
import yaml

from icecube import icetray
from icecube.filterscripts import filter_globals

from utils import get_run_folder

from resources.pulses import GetPulses
from resources.pulses import MergeOversampledEvents


@click.command()
@click.argument('cfg', type=click.Path(exists=True))
@click.argument('run_number', type=int)
@click.option('--scratch/--no-scratch', default=True)
def main(cfg, run_number, scratch):
    start_time = time.time()
    
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

    # get reco pulses
    if cfg['get_reco_pulses']:
        tray.AddSegment(
            GetPulses, "GetPulses", decode=False, simulation=True,
        )
    
    if cfg['get_mc_pulses']:
        default_get_mc_pulses_kwargs = {
            'I3MCPESeriesMap': 'I3MCPESeriesMap',
            'OutputKey': 'MCPulses',
            'CreatePFrames': True,
        }
        if ('get_mc_pulses_kwargs' in cfg['get_mc_pulses_kwargs'] and 
                cfg['get_mc_pulses_kwargs'] is not None):
            default_get_mc_pulses_kwargs.update(cfg['get_mc_pulses_kwargs'])
            
        tray.AddModule(
            GetMCPulses, "GetMCPulses", **default_get_mc_pulses_kwargs
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
    if cfg['oversampling_factor'] is not None:
        if cfg['oversampling_merge_events']:
            
            if cfg['get_reco_pulses'] and cfg['get_mc_pulses']:
                raise NotImplementedError(
                    'Merging of multiple pulses at the same time is not implemented yet'
            )
                
            if cfg['get_reco_pulses']:
                tray.AddModule(
                    MergeOversampledEvents, 'MergeOversampledEvents',
                    OversamplingFactor=cfg['oversampling_factor'],
                    PulseKey='InIceDSTPulses',
                )
            if cfg['get_mc_pulses']:
                tray.AddModule(
                    MergeOversampledEvents, 'MergeOversampledEvents',
                    OversamplingFactor=cfg['oversampling_factor'],
                    PulseKey=default_get_mc_pulses_kwargs['OutputKey'],
                )

            
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
    
    # Make space and delete uneeded keys
    if 'keys_to_delete' in cfg and cfg["keys_to_delete"] is not None:
        tray.AddModule('Delete', 'DeleteKeys', keys=cfg["keys_to_delete"])

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

    end_time = time.time()
    print("That took "+str(end_time - start_time)+" seconds.")

if __name__ == '__main__':
    main()
