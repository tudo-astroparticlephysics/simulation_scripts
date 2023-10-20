#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/icetray-start
#METAPROJECT icetray/v1.8.1
from __future__ import division
import click
import yaml
import numpy as np
import glob

from icecube.simprod import segments
from icecube.sim_services.propagation import get_propagators
from icecube.icetray import I3Tray
from icecube import icetray, dataclasses, neutrino_generator, simclasses,dataio

from utils import create_random_services, get_run_folder
from resources import geometry
from resources.oversampling import DAQFrameMultiplier
from resources.import_events import ImportEvents


@click.command()
@click.argument('cfg', type=click.Path(exists=True))
@click.argument('run_number', type=int)
@click.option('--scratch/--no-scratch', default=True)
def main(cfg, run_number, scratch):
    with open(cfg, 'r') as stream:
        cfg = yaml.full_load(stream)
    cfg['run_number'] = run_number
    cfg['run_folder'] = get_run_folder(run_number)
    if scratch:
        outfile = cfg['scratchfile_pattern'].format(**cfg)
    else:
        outfile = cfg['outfile_pattern'].format(**cfg)
    outfile = outfile.replace(' ', '0')

    # ------------------------------
    # get list of files for this run
    # ------------------------------
    import_cfg = cfg['event_import_settings']
    glob_files = import_cfg['input_file_glob_list']
    if isinstance(glob_files, str):
        # single string provided
        files = glob.glob(glob_files.format(run_number=run_number))
    else:
        # list of file globs provided
        files = []
        for file_pattern in glob_files:
            files.extend(glob.glob(file_pattern.format(run_number=run_number)))
    # sort files
    files = sorted(files)
    # ------------------------------

    click.echo('Run: {}'.format(run_number))
    click.echo('Outfile: {}'.format(outfile))
    click.echo('input Files:')
    for file in files:
        click.echo('\t{}'.format(file))

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
    tray = I3Tray()
    tray.context['I3RandomService'] = random_services[0]

    # import events from another I3-file
    tray.AddModule("I3Reader", "reader", FilenameList=files)

    # TODO: Oversampling
    # tray.AddModule(DAQFrameMultiplier, 'DAQFrameMultiplier',
    #                oversampling_factor=cfg['oversampling_factor'],
    #                mctree_keys=[import_cfg['mctree_name']])

    #mctree_keys=[import_cfg['mctree_name']]

    tray.Add('Delete', Keys=["I3MCTree", "MMCTrackList","I3MCTree_preMuonProp_RNGState"])
    tray.Add('Rename', Keys=['I3MCTree_preMuonProp','I3MCTree'])

    nugen_config = cfg["neutrino_generator_config"]
    propmode = neutrino_generator.to_propagation_mode(nugen_config["propmode"])


    if "CylinderRadiusInM" in nugen_config:
        detcylrad = nugen_config["CylinderRadiusInM"]*icetray.I3Units.m
    else:
        detcylrad = 950*icetray.I3Units.m
    if "CylinderHeightInM" in nugen_config:
        detcyllen = nugen_config["CylinderHeightInM"]*icetray.I3Units.m
    else:
        detcyllen = 1900*icetray.I3Units.m

    origin_x = 0.*icetray.I3Units.m
    origin_y = 0.*icetray.I3Units.m
    origin_z = 0.*icetray.I3Units.m
    cylinderparams = [detcylrad,detcyllen,origin_x,origin_y,origin_z]

    tray.AddService("I3EarthModelServiceFactory", "EarthModelService",
                    EarthModels = nugen_config["earth"],
                    MaterialModels = nugen_config["material"],
                    IceCapType = nugen_config["icecapmodel"])

    tray.AddService("I3NuGSteeringFactory", "NuGSteer",
                    EarthModelName = "EarthModelService",
                    NEvents = -1,
                    CylinderParams = cylinderparams
                    )

    selectNeutrinoCfg = cfg["select_neutrino_cfg"]

    tray.Add('I3NuGSourceSelector', **selectNeutrinoCfg if selectNeutrinoCfg else {})
    
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


if __name__ == '__main__':
    main()
