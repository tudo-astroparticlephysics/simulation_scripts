#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/icetray-start
#METAPROJECT icetray/v1.10.0
import time
import click
import yaml

from icecube import icetray, dataclasses
from icecube import sim_services  # for 'I3MCPEMerger'
from icecube.filterscripts import filter_globals

from utils import get_run_folder

from resources.pulses import (
    GetMCPulses,
    GetPulses,
    MoveSuperDST,
    MergeOversampledEvents,
    MergePulsesNearbyInTime,
    CompressPulses,
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

        if "add_no_noise_pulses" in cfg and cfg["add_no_noise_pulses"]:

            # save the original MC and reco pulses, rename MCPulses to utilize
            tray.Add(
                "Rename", "RenameToRunWithoutNoise",
                Keys=[
                    # temporarily save original keys somwhere else
                    "I3MCPulseSeriesMap", "temp_I3MCPulseSeriesMap",
                    "I3MCPulseSeriesMapParticleIDMap", "temp_I3MCPulseSeriesMapParticleIDMap",
                    "IceTopRawData", "temp_IceTopRawData",
                    "InIceRawData", "temp_InIceRawData",
                    "BeaconLaunches", "temp_BeaconLaunches",
                    "I3TriggerHierarchy", "temp_I3TriggerHierarchy",
                    "I3Triggers", "temp_I3Triggers",

                    # Move WithoutNoise to the original key
                    "I3MCPulseSeriesMapWithoutNoise", "I3MCPulseSeriesMap",
                    "I3MCPulseSeriesMapWithoutNoiseParticleIDMap", "I3MCPulseSeriesMapParticleIDMap",
                    "IceTopRawDataWithoutNoise", "IceTopRawData",
                    "InIceRawDataWithoutNoise", "InIceRawData",
                    "BeaconLaunchesWithoutNoise", "BeaconLaunches",
                    "I3TriggerHierarchyWithoutNoise", "I3TriggerHierarchy",
                    "I3TriggersWithoutNoise", "I3Triggers",
                ],
            )

            # run reco pulses without noise
            tray.AddSegment(
                GetPulses, "GetPulsesWithoutNoise",
                decode=False, simulation=True,
            )

            # move the created Pulses
            tray.AddModule(
                MoveSuperDST, "MoveSuperDST",
                InputKey="I3SuperDST",
                OutputKeyPattern="{}WithoutNoise",
            )

            # rename the created reco pulses, and revert previous ones
            tray.Add(
                "Rename", "RenameToRevert",
                Keys=[
                    # rename created Pulses
                    "MCPulses", "MCPulsesWithoutNoise",
                    "BadDomsList", "BadDomsListWithoutNoise",
                    "BadDomsListSLC", "BadDomsListSLCWithoutNoise",
                    # move without noise keys back
                    "I3MCPulseSeriesMap", "I3MCPulseSeriesMapWithoutNoise",
                    "I3MCPulseSeriesMapParticleIDMap", "I3MCPulseSeriesMapWithoutNoiseParticleIDMap",
                    "IceTopRawData", "IceTopRawDataWithoutNoise",
                    "InIceRawData", "InIceRawDataWithoutNoise",
                    "BeaconLaunches", "BeaconLaunchesWithoutNoise",
                    "I3TriggerHierarchy", "I3TriggerHierarchyWithoutNoise",
                    "I3Triggers", "I3TriggersWithoutNoise",
                    # revert original keys
                    "temp_I3MCPulseSeriesMap", "I3MCPulseSeriesMap",
                    "temp_I3MCPulseSeriesMapParticleIDMap", "I3MCPulseSeriesMapParticleIDMap",
                    "temp_IceTopRawData", "IceTopRawData",
                    "temp_InIceRawData", "InIceRawData",
                    "temp_BeaconLaunches", "BeaconLaunches",
                    "temp_I3TriggerHierarchy", "I3TriggerHierarchy",
                    "temp_I3Triggers", "I3Triggers",
                ],
            )

            # delete created keys that we don't need anymore
            tray.Add(
                "Delete", "DeleteNoNoisePulses",
                Keys=[
                    "CalibratedIceTopATWD_HLC",
                    "CalibratedIceTopATWD_SLC",
                    "CalibratedIceTopFADC_HLC",
                    "CalibratedWaveformRange",
                    "CalibratedWaveforms",
                    "CleanIceTopRawData",
                    "CleanInIceRawData",
                    "ClusterCleaningExcludedTanks",
                    "DSTTriggers",
                    "GetPulsesBaseProc_SimTrimmer_HighCharge",
                    "GetPulsesBaseProc_SimTrimmer_I3SuperDST_CalibratedWaveforms_Borked",
                    "GetPulsesBaseProc_SimTrimmer_I3SuperDST_CalibratedWaveforms_Chi",
                    "HLCTankPulses",
                    "I3SuperDST",
                    "IceTopCalibratedWaveformRange",
                    "IceTopCalibratedWaveforms",
                    "IceTopDSTPulses",
                    "IceTopHLCPulseInfo",
                    "IceTopHLCVEMPulses",
                    "IceTopPulses",
                    "IceTopPulses_HLC",
                    "IceTopPulses_SLC",
                    "IceTopRawData",
                    "IceTopSLCVEMPulses",
                    "InIceDSTPulses",
                    "QTriggerHierarchy",
                    "RawDSTPulses",
                    "SLCTankPulses",
                    "SimTrimmer",
                    "TankPulseMergerExcludedTanks",
                    "TankPulseMergerExcludedTanksSLC",
                    "UncleanedInIcePulses",
                    "UncleanedInIcePulsesTimeRange",
                ],
            )

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
            if time_window is not None and mcpe_series is not None:
                tray.AddModule(
                    "I3MCPEMerger", f"I3MCPEMerger_{i:03d}",
                    Input=mcpe_series,
                    Output=mcpe_series,
                    timeWindow=time_window * icetray.I3Units.ns,
                )

            if time_window is not None:
                default_get_mc_pulses_kwargs["PulseWidth"] = time_window * icetray.I3Units.ns

            tray.AddModule(
                GetMCPulses, f"GetMCPulses_{i:03d}", **default_get_mc_pulses_kwargs
            )

    # merge reco and mc pulses in time
    if "merge_pulses_cfg" in cfg and cfg["merge_pulses_cfg"]:
        print("Merging pulses:", cfg["merge_pulses_cfg"])
        tray.AddModule(
            MergePulsesNearbyInTime, "MergePulsesNearbyInTime",
            **cfg["merge_pulses_cfg"]
        )
        if "OutputKeys" in cfg["merge_pulses_cfg"]:
            output_keys += cfg["merge_pulses_cfg"]["OutputKeys"]

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
            for pulse_key in output_keys:
                tray.AddModule(
                    MergeOversampledEvents, 'MergeOversampledEvents',
                    OversamplingFactor=cfg['oversampling_factor'],
                    PulseKey=pulse_key,
                )

    if "compress_pulses_cfg" in cfg and cfg["compress_pulses_cfg"]:
        tray.AddModule(
            CompressPulses, "CompressPulses",
            **cfg["compress_pulses_cfg"]
        )
        if "OutputKeys" in cfg["compress_pulses_cfg"]:
            output_keys += cfg["compress_pulses_cfg"]["OutputKeys"]

    keys_to_keep = [
        'TimeShift',
        'I3MCTree_preMuonProp',
        'I3MCTree',
        'I3MCWeightDict',
        'CorsikaWeightMap',
        'MMCTrackList',
        'I3EventHeader',
        'I3SuperDST',
        'I3SuperDSTWithoutNoise',
        'RNGState',
        'oversampling',
        'AggregatedPulses',
        'InIceDSTPulses',
        'InIceDSTPulsesWithowtNoise',
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
