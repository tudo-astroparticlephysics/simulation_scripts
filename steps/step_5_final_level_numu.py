#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/eganster/combo.releases.V01-00-02.py3-v4.0.1.RHEL_7_x86_64
import click
import yaml
import time 
import numpy as np 

from I3Tray import I3Tray
from icecube import icetray, dataclasses, dataio

# L3 processing
from icecube.level3_filter_muon.MuonL3TraySegment import MuonL3
# L4 processing
from icecube.finallevel_filter_diffusenumu import level4
# L5 processing
from icecube.finallevel_filter_diffusenumu import level5
# Post-L5 processing
from icecube.finallevel_filter_diffusenumu import post_level5
# hdf writer
from icecube.finallevel_filter_diffusenumu.write_hdf import write_hdf
# import mapping of SnowStorm parameters
from icecube.finallevel_filter_diffusenumu import snobo_parameters


from utils import get_run_folder


@click.command()
@click.argument('cfg', type=click.Path(exists=True))
@click.argument('run_number', type=int)
@click.option('--scratch/--no-scratch', default=True)
def main(cfg, run_number, scratch):
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
    infile = infile.replace('Level0.{}'.format(cfg['previous_step'] % 10),
                            'Level2')
    infile = infile.replace('2012_pass2', 'pass2')

    if scratch:
        outfile = cfg['scratchfile_pattern'].format(**cfg)
    else:
        outfile = cfg['outfile_pattern'].format(**cfg)
    outfile = outfile.replace('Level0.{}'.format(cfg['step']), 'Level5')
    outfile = outfile.replace(' ', '0')
    outfile = outfile.replace('2012_pass2', 'pass2')
    print('Outfile != $FINAL_OUT clean up for crashed scripts not possible!')

    tray = I3Tray()
    
    start_time = time.time()
    
    """The main L1 script"""
    tray.AddModule('I3Reader',
                   'i3 reader',
                   FilenameList=[cfg['gcd_pass2'], infile])

    # Default options
    options = {
        'photonicsdir': '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/SPICEMie/',
        'photonicsdriverdir': '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/SPICEMie/driverfiles',
        'photonicsdriverfile': 'mu_photorec.list',
        'infmuonampsplinepath': '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfBareMu_mie_abs_z20a10_V2.fits',
        'infmuonprobsplinepath': '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfBareMu_mie_prob_z20a10_V2.fits',
        'cascadeampsplinepath': '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/ems_mie_z20_a10.abs.fits',
        'cascadeprobsplinepath': '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/ems_mie_z20_a10.prob.fits',
        'restoretwformc': False,
        'stochampsplinepath': '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfHighEStoch_mie_abs_z20a10.fits',
        'stochprobsplinepath': '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfHighEStoch_mie_prob_z20a10.fits',
        'do_postL5': False,
        'is_MC': True,
    }
    if 'DiffuseNuMuFinalLevelSettings' in cfg:
        options.update(cfg['DiffuseNuMuFinalLevelSettings'])

    # L3 processing
    tray.AddSegment(MuonL3,
                    gcdfile=None,
                    infiles=None,
                    output_i3="",
                    output_hd5="",
                    output_root="",
                    photonicsdir=options['photonicsdir'],
                    photonicsdriverdir=options['photonicsdriverdir'],
                    photonicsdriverfile=options['photonicsdriverfile'],
                    infmuonampsplinepath=options['infmuonampsplinepath'],
                    infmuonprobsplinepath=options['infmuonprobsplinepath'],
                    cascadeampsplinepath=options['cascadeampsplinepath'],
                    cascadeprobsplinepath=options['cascadeprobsplinepath'],
                    restore_timewindow_forMC=options['restoretwformc']
                    )

    # L4 processing
    tray.AddSegment(level4.IC12L4,
                    gcdfile=None,
                    infiles=None,
                    table_paths=options,
                    is_numu=False,
                    pulsemap_name="TWSRTHVInIcePulsesIC")
    # L5 processing
    tray.Add(level5.segments.Scorer,
             CutFunc=level5.segments.CutFunc,
             CascCut=0.5)
    # millipede
    tray.Add(level5.segments.millipede_segment, "MillipedeLosses",
             table_paths=options)
    # paraboloid
    tray.Add(level5.segments.paraboloid_segment, "Paraboloid",
             table_paths=options,
             pulses="TWSRTHVInIcePulsesIC")

    # postL5 (Renes PS-recos and PassedHESE bool)
    if options["do_postL5"]:
        tray.Add(post_level5.pass1_ps_reco_paraboloid, "PostDiffusePSRecoSplineMPE",
                 PulsesName          = "TWSRTHVInIcePulsesIC",
                 configuration       = "max",
                 TrackSeedList       = ["SplineMPEIC"],
                 EnergyEstimators    = ["SplineMPEICTruncatedEnergySPICEMie_AllDOMS_Muon"],
                )
        if options["is_MC"]:
            tray.Add(post_level5.muon_energies, "getMuonEnergies")
        tray.Add(post_level5.add_hese_tag, "AddHESETag")

    # add the Snowstorm parameters as I3Doubles to the frame:
    tray.AddModule(snobo_parameters.map_parameters, "SnowstormParameterMapper",
                   Streams=[icetray.I3Frame.Stream('M')])

    # Write output
    tray.AddModule("I3Writer", "EventWriter",
                   filename=outfile,
                   Streams=[icetray.I3Frame.DAQ,
                            icetray.I3Frame.Physics,
                            icetray.I3Frame.TrayInfo,
                            icetray.I3Frame.Simulation,
                            icetray.I3Frame.Stream('M')],
                   DropOrphanStreams=[icetray.I3Frame.DAQ])
    tray.AddModule("TrashCan", "the can")
    tray.Execute()
    tray.Finish()
    del tray

    end_time = time.time()
    print('That took {} seconds'.format(np.round(end_time - start_time, 1)))

if __name__ == '__main__':
    main()
