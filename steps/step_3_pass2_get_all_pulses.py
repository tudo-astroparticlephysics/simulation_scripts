#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/icetray-start
#METAPROJECT icetray/v1.10.0
import time
import click
import yaml

from icecube import icetray
from icecube import sim_services  # for 'I3MCPEMerger'
from icecube.filterscripts import filter_globals

from utils import get_run_folder

from resources.pulses import (
    GetMCPulses,
    GetPulses,
    MergeOversampledEvents,
)


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

    tray = icetray.I3Tray()
    tray.AddModule('I3Reader',
                   'i3 reader',
                   FilenameList=[cfg['gcd_pass2'], infile])

    # get reco pulses
    if cfg['get_reco_pulses']:
        tray.AddSegment(
            GetPulses, "GetPulses", decode=False, simulation=True,
        )
    
    # get mc pulses
    output_keys = []
    if cfg['get_mc_pulses']:
        for i, mc_pulses_cfg in enumerate(cfg['get_mc_pulses_kwargs_list']):
            default_get_mc_pulses_kwargs = {
                'I3MCPESeriesMap': 'I3MCPESeriesMap',
                'TimeWindowCompression': 1,
                'I3MCPEorPulseSeriesMap': 'I3MCPESeriesMap',
                'OutputKey': 'MCPEs',
                'CreatePFrames': False,
                'WriteToQFrame': True,
            }
            default_get_mc_pulses_kwargs.update(mc_pulses_cfg)
            output_keys.append(default_get_mc_pulses_kwargs['OutputKey'])
            
            # compress I3MCPESeriesMap
            time_window = default_get_mc_pulses_kwargs.pop("TimeWindowCompression")
            mcpe_series = default_get_mc_pulses_kwargs.pop('I3MCPESeriesMap')
            if time_window is not None:
                tray.AddModule(
                    "I3MCPEMerger", f"I3MCPEMerger_{i:03d}",
                    Input=mcpe_series, 
                    Output=mcpe_series,
                    timeWindow=time_window * icetray.I3Units.ns,
                )
            
            tray.AddModule(
                GetMCPulses, f"GetMCPulses_{i:03d}", **default_get_mc_pulses_kwargs
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
        'DOMDeadTimesMC',
        'SplitUncleanedInIcePulses',
        'SplitUncleanedInIcePulsesTimeRange',
        'SplitUncleanedInIceDSTPulsesTimeRange',
        'I3TriggerHierarchy',
        'GCFilter_GCFilterMJD',
        ]
    keys_to_keep += filter_globals.inice_split_keeps
    keys_to_keep += filter_globals.onlinel2filter_keeps
    keys_to_keep += output_keys

    if 'event_import_settings' in cfg:
        import_keys = [cfg['event_import_settings']['rename_dict'].get(k, k)
                       for k in cfg['event_import_settings']['keys_to_import']]
        keys_to_keep += import_keys

    tray.AddModule("Keep", "keep_before_merge",
                   keys=keys_to_keep + cfg['oversampling_keep_keys'])
    
    # Make space and delete uneeded keys
    if 'keys_to_delete' in cfg and cfg["keys_to_delete"] is not None:
        tray.AddModule('Delete', 'DeleteKeys', keys=cfg["keys_to_delete"])
    
    
    if "i3_streams" in cfg and cfg["i3_streams"] is not None:
        i3_streams = [
            icetray.I3Frame.Stream(s) for s in cfg["i3_streams"]
        ]
    else:
        i3_streams = [
            icetray.I3Frame.DAQ,
            icetray.I3Frame.Physics,
            icetray.I3Frame.TrayInfo,
            icetray.I3Frame.Simulation,
            icetray.I3Frame.Stream("M"),
        ]

    tray.AddModule("I3Writer", "EventWriter",
                   filename=outfile,
                   Streams=i3_streams)
    tray.AddModule("TrashCan", "the can")
    tray.Execute()
    tray.Finish()

    end_time = time.time()
    print("That took "+str(end_time - start_time)+" seconds.")

if __name__ == '__main__':
    main()
