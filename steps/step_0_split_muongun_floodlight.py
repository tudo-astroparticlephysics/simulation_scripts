#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/icetray-start
#METAPROJECT /mnt/lfs7/user/mhuennefeld/software/icecube/py3-v4.1.0/combo_V01-00-00-RC0/build
import click
import yaml

import time

import numpy as np

from icecube.simprod import segments

from I3Tray import I3Tray
from icecube import icetray, dataclasses
from icecube import sim_services, MuonGun

from utils import create_random_services, get_run_folder
from dom_distance_cut import OversizeSplitterNSplits, generate_stream_object

# import muon split
from resources.muon_split_module import SplitModule
from resources.muon_split_functions import selection
# from resources.server_new_tree_hull import SelectionModule

def start_timer(frame):
    frame['starttime'] = dataclasses.I3Double(time.time())

def end_timer(frame):
    frame['endtime'] = dataclasses.I3Double(time.time())
    print('selection time {}'.format(frame['endtime'].value - frame['starttime'].value))

def start_timer2(frame):
    frame['starttime2'] = dataclasses.I3Double(time.time())

def end_timer2(frame):
    frame['endtime2'] = dataclasses.I3Double(time.time())
    print('split time {}'.format(frame['endtime2'].value - frame['starttime2'].value))


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
    if cfg['distance_splits'] is not None:
        click.echo('SplittingDistances: {}'.format(cfg['distance_splits']))
        click.echo('Oversizefactors: {}'.format(cfg['oversize_factors']))
    click.echo('NEvents: {}'.format(cfg['n_events_per_run']))
    click.echo('EMin: {}'.format(cfg['e_min']))
    click.echo('EMax: {}'.format(cfg['e_max']))
    click.echo('EBreak: {}'.format(cfg['muongun_e_break']))
    click.echo('Gamma: {}'.format(cfg['gamma']))
    click.echo('ZenithMin: {}'.format(cfg['zenith_min']))
    click.echo('ZenithMax: {}'.format(cfg['zenith_max']))
    click.echo('NewPsi: {}'.format(cfg['new_psi']))
    click.echo('PercentagEnergyLoss: {}'.format(cfg['percentage_energy_loss']))

    starttime = time.time()

    tray = I3Tray()

    random_services, _ = create_random_services(
        dataset_number=cfg['dataset_number'],
        run_number=cfg['run_number'],
        seed=cfg['seed'],
        n_services=2,
        use_gslrng=True)

    random_service, random_service_prop = random_services
    tray.context['I3RandomService'] = random_service

    tray.AddModule("I3InfiniteSource",
                   "TheSource",
                   Prefix=cfg['gcd'],
                   Stream=icetray.I3Frame.DAQ)

    tray.AddSegment(
        segments.GenerateSingleMuons,
        "GenerateCosmicRayMuons",
        NumEvents=cfg['n_events_per_run'],
        FromEnergy=cfg['e_min'] * icetray.I3Units.GeV,
        ToEnergy=cfg['e_max'] * icetray.I3Units.GeV,
        BreakEnergy=cfg['muongun_e_break'] * icetray.I3Units.GeV,
        GammaIndex=cfg['gamma'],
        ZenithRange=[cfg['zenith_min'] * icetray.I3Units.deg,
                     cfg['zenith_max'] * icetray.I3Units.deg])
    
    tray.AddSegment(
        segments.PropagateMuons,
        "PropagateMuons",
        RandomService=random_service_prop,
        **cfg['muon_propagation_config'])

    # Muon selection
    # tray.Add(start_timer, 'stimer', Streams=[icetray.I3Frame.DAQ])
    tray.Add(selection, 'Selector', percentage_energy_loss=cfg['percentage_energy_loss'], Streams=[icetray.I3Frame.DAQ])
    # tray.AddModule(SelectionModule, 'selectmodule', ELoss=0.01)
    # tray.Add(end_timer, 'etimer', Streams=[icetray.I3Frame.DAQ])

	# Add muon split module 
    # tray.AddModule(start_timer2, 'stimer2', Streams=[icetray.I3Frame.DAQ])
    tray.AddModule(SplitModule, 'SplitMuonTrack', NewPsi=cfg['new_psi'], RandomSeed=42)	
    # tray.AddModule(end_timer2, 'etimer2', Streams=[icetray.I3Frame.DAQ])


    if scratch:
        outfile = cfg['scratchfile_pattern'].format(**cfg)
    else:
        outfile = cfg['outfile_pattern'].format(**cfg)
    outfile = outfile.replace(' ', '0')
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
    click.echo('Scratch: {}'.format(scratch))
    
    tray.AddModule("TrashCan", "the can")
    tray.Execute()
    Endtime = time.time()
    print('time after execute: {}'.format(Endtime - starttime))

    tray.Finish()


if __name__ == '__main__':
    main()
