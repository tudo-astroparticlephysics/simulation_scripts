#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/icetray-start
#METAPROJECT /data/user/mhuennefeld/software/icecube/py3-v4.1.0/combo_V01-00-00-RC0/build
#--METAPROJECT /data/user/mhuennefeld/software/icecube/py3-v4.1.0/combo_V01-00-00/build
#--METAPROJECT combo/V01-00-00 # <-- Causes segfaults, therefore use RC0
import os
import sys
if "ENV_SITE_PACKAGES" in os.environ:
    sys.path.insert(1, os.environ["ENV_SITE_PACKAGES"])

import click
import yaml

from icecube.simprod import segments
from I3Tray import I3Tray
from icecube import icetray, dataclasses, dataio, phys_services
from utils import create_random_services, get_run_folder

from resources.pulses import RemoveEarlyPulses

# Load libraries
from icecube import clsim
from icecube import phys_services
from icecube import sim_services
from icecube import vuvuzela
from icecube import DOMLauncher
from icecube import trigger_sim
import time

MCPE_SERIES_MAP = 'I3MCPESeriesMap'


@click.command()
@click.argument('cfg', type=click.Path(exists=True))
@click.argument('run_number', type=int)
@click.option('--scratch/--no-scratch', default=True)
def main(cfg, run_number, scratch):
    with open(cfg, 'r') as stream:
        cfg = yaml.full_load(stream)
    cfg['run_number'] = run_number
    cfg['run_folder'] = get_run_folder(run_number)
    infile = cfg['infile_pattern'].format(**cfg)
    infile = infile.replace(' ', '0')
    infile = infile.replace('Level0.{}'.format(cfg['previous_step']),
                            'Level0.{}'.format(cfg['previous_step'] % 10))

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

    start_time = time.time()

    tray.context['I3FileStager'] = dataio.get_stagers()

    random_services, run_id = create_random_services(
        dataset_number=cfg['dataset_number'],
        run_number=cfg['run_number'],
        seed=cfg['seed'],
        n_services=2,
        use_gslrng=cfg['random_service_use_gslrng'])
    tray.context['I3RandomService'] = random_services[0]
    tray.context['I3RandomServiceForNoiseFree'] = random_services[1]

    tray.Add('I3Reader', FilenameList=[cfg['gcd_pass2'], infile])

    """
    Perform Detector simulation:
        https://code.icecube.wisc.edu/projects/icecube/browser/IceCube/
        meta-projects/combo/stable/simprod-scripts/python/segments/
        DetectorSim.py
    """

    # Combine MCPEs from both detectors
    if cfg['det_is_genie_simulation']:
        tray.Add("Rename", Keys=[MCPE_SERIES_MAP, 'GenieMCPEs'])
        tray.Add("I3CombineMCPE",
                 InputResponses=["GenieMCPEs", "BackgroundMCPEs"],
                 OutputResponse=MCPE_SERIES_MAP)
        tray.Add("Delete", Keys=['BackgroundMCPEs', 'GenieMCPEs'])
    if cfg['det_is_icetop_simulation']:
        tray.Add("Rename", Keys=[MCPE_SERIES_MAP, 'InIceMCPEs'])
        tray.Add("I3CombineMCPE",
                 InputResponses=["IceTopMCPEs", "InIceMCPEs"],
                 OutputResponse=MCPE_SERIES_MAP)
        tray.Add("Delete", Keys=['InIceMCPEs', 'IceTopMCPEs'])

    # Sample a different efficiency
    sample_eff = cfg['det_dom_eff_resmapling_sample_efficiency']
    generated_eff = cfg['det_dom_eff_resmapling_generated_efficiency']
    if sample_eff > 0.0:
        if sample_eff > generated_eff:
            msg = 'Cannot upscale from GeneratedEfficiency %s to '
            msg += 'SampleEfficiency %s'
            icecube.icetray.logging.log_fatal(
               msg % (sample_eff, generated_eff))
        tray.AddSegment(segments.MultiDomEffSample, "resample",
                        GeneratedEfficiency=generated_eff,
                        SampleEfficiencies=[sample_eff],
                        InputSeriesName=MCPE_SERIES_MAP,
                        DeleteOriginalSeries=True,
                        OverwriteOriginalSeries=True)

    if run_number < cfg['det_keep_all_upto']:
        cfg['det_keep_mc_hits'] = True
        cfg['det_keep_propagated_mc_tree'] = True
        cfg['det_keep_mc_pulses'] = True

    TimeShiftSkipKeys = [
        "SnowstormParameterRanges",
        "SnowstormParameters",
        "SnowstormParametrizations",
        "SnowstormProposalDistribution",
        "WavelengthAcceptance",
        "WavelengthGenerationBias",
        "LeptonInjectorProperties",
        "EventProperties",
        "MediumProperties",
    ]

    if "add_no_noise_pulses" in cfg and cfg["add_no_noise_pulses"]:
        if cfg["det_skip_noise_generation"]:
            raise ValueError(
                "Do not use add_no_noise_pulses with "
                "det_skip_noise_generation"
            )

        # copy original MCPEs
        tray.Add("Copy", Keys=[MCPE_SERIES_MAP, "MCPEsOriginal"])

        tray.AddSegment(segments.DetectorSim, "DetectorSimWithoutNoise",
            RandomService='I3RandomServiceForNoiseFree',
            RunID=run_id,
            GCDFile=cfg['gcd_pass2'],
            KeepMCHits=cfg['det_keep_mc_hits'],
            KeepPropagatedMCTree=cfg['det_keep_propagated_mc_tree'],
            KeepMCPulses=cfg['det_keep_mc_pulses'],
            SkipNoiseGenerator=True,
            LowMem=cfg['det_low_mem'],
            InputPESeriesMapName=MCPE_SERIES_MAP,
            BeaconLaunches=cfg['det_add_beacon_launches'],
            FilterTrigger=cfg['det_filter_trigger'],
            TimeShiftSkipKeys=TimeShiftSkipKeys,
        )
        tray.Add(
            "Rename", "RenameNoNoisePulses",
            Keys=[
                "I3MCPulseSeriesMap", "I3MCPulseSeriesMapWithoutNoise",
                "I3MCPulseSeriesMapParticleIDMap", "I3MCPulseSeriesMapWithoutNoiseParticleIDMap",
                "IceTopRawData", "IceTopRawDataWithoutNoise",
                "InIceRawData", "InIceRawDataWithoutNoise",
                "BeaconLaunches", "BeaconLaunchesWithoutNoise",
                "I3EventHeader", "I3EventHeaderWithoutNoise",
                "I3TriggerHierarchy", "I3TriggerHierarchyWithoutNoise",
                "I3Triggers", "I3TriggersWithoutNoise",
                "TimeShift", "TimeShiftWithoutNoise",
                "DOMSets", "DOMSetsWithoutNoise",
            ],
        )
        # restore original MCPEs
        tray.Add("Delete", Keys=[MCPE_SERIES_MAP, MCPE_SERIES_MAP + "WithoutNoise"])
        tray.Add("Rename", Keys=["MCPEsOriginal", MCPE_SERIES_MAP])

    tray.AddSegment(segments.DetectorSim, "DetectorSim",
                    RandomService='I3RandomService',
                    RunID=run_id,
                    GCDFile=cfg['gcd_pass2'],
                    KeepMCHits=cfg['det_keep_mc_hits'],
                    KeepPropagatedMCTree=cfg['det_keep_propagated_mc_tree'],
                    KeepMCPulses=cfg['det_keep_mc_pulses'],
                    SkipNoiseGenerator=cfg['det_skip_noise_generation'],
                    LowMem=cfg['det_low_mem'],
                    InputPESeriesMapName=MCPE_SERIES_MAP,
                    BeaconLaunches=cfg['det_add_beacon_launches'],
                    FilterTrigger=cfg['det_filter_trigger'],
                    TimeShiftSkipKeys=TimeShiftSkipKeys,
                    )

    if "add_no_noise_pulses" in cfg and cfg["add_no_noise_pulses"]:
        class UpdateTimeShift(icetray.I3Module):
            def DAQ(self, frame):

                # combine time shift and remove the old keys
                time_shift1 = frame["TimeShift"].value
                time_shift2 = frame["TimeShiftWithoutNoise"].value
                del frame["TimeShift"]
                del frame["TimeShiftWithoutNoise"]
                frame["TimeShift"] = dataclasses.I3Double(
                    time_shift1 + time_shift2
                )

                # update event header
                if "I3EventHeaderWithoutNoise" in frame:
                    event_header = dataclasses.I3EventHeader(
                        frame["I3EventHeaderWithoutNoise"]
                    )
                    event_header.start_time -= time_shift1
                    event_header.end_time -= time_shift1
                    event_header.event_id = frame["I3EventHeader"].event_id
                    del frame["I3EventHeaderWithoutNoise"]
                else:
                    event_header = dataclasses.I3EventHeader(
                        frame["I3EventHeader"]
                    )
                frame["I3EventHeaderWithoutNoise"] = event_header

                self.PushFrame(frame)
        tray.AddModule(UpdateTimeShift, "UpdateTimeShift")

        # it seems that adding the noise-free pulses in parallel
        # has some effect on the time shift and injected noise.
        # We remove the early pulses (<-512ns) here to avoid this.
        # (I3SuperDST can only express times >-512ns)
        tray.AddModule(
            RemoveEarlyPulses, "RemoveEarlyPulses",
            InputKeys=[
                "InIceRawData",
                "InIceRawDataWithoutNoise",
            ],
        )

    if not cfg['det_keep_mc_pulses']:
        tray.Add("Delete", Keys=['I3MCPulseSeriesMapPrimaryIDMap'])

    if not cfg['det_keep_mc_hits']:
        tray.Add("Delete", Keys=['I3MCPESeriesMapParticleIDMap'])

    if cfg['det_remove_keys_from_m_frame']:
        tray.Add("Delete", Keys=cfg['det_remove_keys_from_m_frame'])

    if cfg['det_convert_to_linear_tree']:
        tray.AddModule(segments.ConvertToLinearizedMCTree,
                       "lineartree", streams=[icetray.I3Frame.DAQ])

    tray.AddModule("I3Writer", "EventWriter",
                   filename=outfile,
                   Streams=[icetray.I3Frame.DAQ,
                            icetray.I3Frame.Physics,
                            icetray.I3Frame.TrayInfo,
                            icetray.I3Frame.Simulation,
                            icetray.I3Frame.Stream('m'),
                            icetray.I3Frame.Stream('M')])
    tray.AddModule("TrashCan", "the can")
    tray.Execute()
    tray.Finish()

    end_time = time.time()
    print("That took "+str(end_time - start_time)+" seconds.")


if __name__ == '__main__':
    main()
