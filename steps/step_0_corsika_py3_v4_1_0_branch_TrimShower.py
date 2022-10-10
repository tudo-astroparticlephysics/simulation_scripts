#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/icetray-start
#METAPROJECT /data/user/mhuennefeld/software/icecube/py3-v4.1.0/icetray_trimshower_branch/build
import os
import click
import yaml

import numpy as np

from icecube.simprod import segments

from I3Tray import I3Tray
from icecube import icetray, dataclasses
from icecube import sim_services, MuonGun

from utils import (
    create_random_services,
    create_random_services_settings,
    get_run_folder,
    load_class,
)
from resources.biased_simulation import BaseSimulationBias
from dom_distance_cut import OversizeSplitterNSplits, generate_stream_object


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
    if scratch:
        outfile = cfg['scratchfile_pattern'].format(**cfg)
    else:
        outfile = cfg['outfile_pattern'].format(**cfg)
    outfile = outfile.replace(' ', '0')

    # Print out CORSIKA information
    click.echo('CORSIKA module: {}'.format(cfg['corsika_module']))
    corsika_settings = cfg['corsika_settings']
    click.echo('CORSIKA settings:')
    for name, value in corsika_settings.items():
        click.echo('  {}: {}'.format(name, value))

    if 'ApplyBaseSimulationBias' in cfg and cfg['ApplyBaseSimulationBias']:
        click.echo('Apply simulation bias: True')
    else:
        click.echo('Apply simulation bias: False')

    tray = I3Tray()

    if 'ApplyBaseSimulationBias' in cfg and cfg['ApplyBaseSimulationBias']:
        n_services = 3
    else:
        n_services = 2

    # --------------------------------------------
    # Create random number generators and settings
    # --------------------------------------------
    random_services_settings, _ = create_random_services_settings(
        dataset_number=cfg['dataset_number'],
        run_number=cfg['run_number'],
        seed=cfg['seed'],
        n_services=n_services,
        use_gslrng=cfg['random_service_use_gslrng'],
    )

    random_services, _ = create_random_services(
        dataset_number=cfg['dataset_number'],
        run_number=cfg['run_number'],
        seed=cfg['seed'],
        n_services=n_services,
        use_gslrng=cfg['random_service_use_gslrng'],
    )
    random_service, random_service_prop = random_services[:2]
    tray.context['I3RandomService'] = random_service
    corsika_rng_cfg = random_services_settings[0]

    # -----------------------------------
    # Run script to generate CORSIKA file
    # -----------------------------------
    corsika_file = os.path.abspath(outfile + '.corsika_temp.i3')
    module_class = load_class(cfg['corsika_module'])
    module = module_class()

    # add overwrite necessary settings
    corsika_settings.update({

        # random number generator
        'usegslrng': cfg['random_service_use_gslrng'],
        'seed': corsika_rng_cfg['seed'],
        'corsikaseed': corsika_rng_cfg['seed'],

        'gcdfile': cfg['gcd'],
        'outputfile': corsika_file,
        'compress': False,
        'runnum': cfg['run_number'],
        'nshowers': cfg['n_events_per_run'],
    })

    # add additional rng info if not using gslrng
    if not cfg['random_service_use_gslrng']:
        corsika_settings.update({
            # random number generator
            'procnum': corsika_rng_cfg['streamnum'],
            'nproc': corsika_rng_cfg['nstreams'],
        })

    # configure module
    for name, value in corsika_settings.items():
        module.SetParameter(name, value)

    # execute module
    stats = {}
    return_val = module.Execute(stats)
    if return_val != 0:
        raise ValueError(
            'CORSIKA module exited with non-zero return value:', return_val)
    print('Stats')
    print(stats)
    # -----------------------------------

    # read previously created corsika file
    tray.Add('I3Reader', FilenameList=[cfg['gcd'], corsika_file])

    # propagate muons if config exists in config
    # Note: Snowstorm may perform muon propagation internally
    if 'muon_propagation_config' in cfg:
        tray.AddModule('Rename', keys=['I3MCTree', 'I3MCTree_preMuonProp'])
        tray.AddSegment(segments.PropagateMuons,
                        'propagate_muons',
                        RandomService=random_service_prop,
                        **cfg['muon_propagation_config'])

    # Bias simulation if desired
    if 'ApplyBaseSimulationBias' in cfg and cfg['ApplyBaseSimulationBias']:
        tray.AddModule(
            BaseSimulationBias,
            'BaseSimulationBias',
            random_service=random_services[2],
            **cfg['BaseSimulationBiasSettings']
        )

    if cfg['distance_splits'] is not None:
        import dom_distance_cut as dom_cut
        click.echo('Oversizestreams')
        stream_objects = dom_cut.generate_stream_object(
            cut_distances=cfg['distance_splits'],
            dom_limits=cfg['threshold_doms'],
            oversize_factors=cfg['oversize_factors'])
        tray.AddModule(dom_cut.OversizeSplitterNSplits,
                       "OversizeSplitterNSplits",
                       thresholds=cfg['distance_splits'],
                       thresholds_doms=cfg['threshold_doms'],
                       oversize_factors=cfg['oversize_factors'],
                       simulaton_type=cfg['neutrino_flavor'].lower())
        for stream_i in stream_objects:
            outfile_i = stream_i.transform_filepath(outfile)
            click.echo('\t{}'.format(stream_i))
            click.echo('\tOutfile: {}'.format(outfile_i))
            tray.AddModule("I3Writer",
                           "writer_{}".format(stream_i.stream_name),
                           Filename=outfile_i,
                           Streams=[icetray.I3Frame.DAQ,
                                    icetray.I3Frame.Physics,
                                    icetray.I3Frame.Stream('S'),
                                    icetray.I3Frame.Stream('M')],
                           If=stream_i)
    else:
        click.echo('Output: {}'.format(outfile))
        tray.AddModule("I3Writer", "writer",
                       Filename=outfile,
                       Streams=[icetray.I3Frame.DAQ,
                                icetray.I3Frame.Physics,
                                icetray.I3Frame.Stream('S'),
                                icetray.I3Frame.Stream('M')])
    click.echo('Scratch: {}'.format(scratch))
    tray.AddModule("TrashCan", "the can")
    tray.Execute()
    tray.Finish()

    # remove temporarily created CORSIKA file
    print('Cleaning up temp CORSIKA file: {}'.format(corsika_file))
    os.remove(corsika_file)


if __name__ == '__main__':
    main()
