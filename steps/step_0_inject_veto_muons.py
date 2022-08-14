#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/icetray-start
#METAPROJECT /data/user/mhuennefeld/software/icecube/py3-v4.1.0/combo_V01-00-00/build
from __future__ import division
import click
import yaml
import numpy as np
import glob

from icecube.simprod import segments

from I3Tray import I3Tray, I3Units
from icecube import icetray, dataclasses

from utils import create_random_services, get_run_folder
from resources import geometry
from resources.import_events import ImportEvents
from resources.biased_simulation import BaseSimulationBias
from resources.veto_muon import InjectSingleVetoMuon, CombineMCTrees


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

    import_cfg['run_number'] = cfg['run_number']
    import_cfg['dataset_number'] = cfg['dataset_number']
    import_cfg['folder_num_pre_offset'] = cfg['run_number']//1000
    import_cfg['folder_num'] = (
        import_cfg['folder_offset'] + cfg['run_number']//1000
    )
    import_cfg['folder_pattern'] = import_cfg['folder_pattern'].format(
        **import_cfg)
    import_cfg['run_folder'] = import_cfg['folder_pattern'].format(
        **import_cfg)

    if isinstance(glob_files, str):
        # single string provided
        files = glob.glob(glob_files.format(**import_cfg))
    else:
        # list of file globs provided
        files = []
        for file_pattern in glob_files:
            files.extend(glob.glob(file_pattern.format(**import_cfg)))

    # sort files
    files = sorted(files)
    # ------------------------------

    click.echo('Run: {}'.format(run_number))
    click.echo('Outfile: {}'.format(outfile))
    click.echo('Keys to import: {}'.format(import_cfg['keys_to_import']))
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
        n_services=3,
        use_gslrng=cfg['random_service_use_gslrng'])

    # --------------------------------------
    # Build IceTray
    # --------------------------------------
    tray = I3Tray()

    # import events from another I3-file
    if 'num_events' not in import_cfg:
        import_cfg['num_events'] = None

    tray.AddModule(
        ImportEvents,
        'ImportEvents',
        files=files,
        num_events=import_cfg['num_events'],
        keys_to_import=import_cfg['keys_to_import'],
        rename_dict=import_cfg['rename_dict'],
        mctree_name=import_cfg['mctree_name'],
    )

    # inject coincident muon from some direction as neutrino
    if 'mctree_name' not in cfg['veto_muon_injection_config']:
        cfg['veto_muon_injection_config']['mctree_name'] = (
            import_cfg['mctree_name']
        )

    tray.AddModule(
        InjectSingleVetoMuon, 'InjectSingleVetoMuon',
        random_service=random_services[0],
        **cfg['veto_muon_injection_config']
    )

    t_name = import_cfg['mctree_name']

    # # rename I3MCTrees so that we can run PROPOSAL
    # def rename_keys(frame, rename_dict):
    #     for key, new_name in rename_dict.items():
    #         if key in frame:
    #             frame[new_name] = frame[key]
    #             del frame[key]
    #         print('key', key)

    # rename_dict = {
    #     t_name: '_' + t_name,
    #     t_name + '_preMuonProp': '_' + t_name + '_preMuonProp',
    #     t_name + '_preMuonProp_RNGState': (
    #         '_' + t_name + '_preMuonProp_RNGState',
    #     ),
    #     t_name + '_veto_muon': 'I3MCTree_preMuonProp',
    # }
    # tray.AddModule(rename_keys, 'TempRename', rename_dict=rename_dict)

    # temporarily save MMCTrackList
    tray.AddModule('Rename', keys=['MMCTrackList', '_MMCTrackList'])

    # propagate injected muon
    cfg['muon_propagation_config'].update({
        'InputMCTreeName': t_name + 'VetoMuon_preMuonProp',
        'OutputMCTreeName': t_name + 'VetoMuon',
    })
    tray.AddSegment(segments.PropagateMuons,
                    'propagate_muons',
                    RandomService=random_services[1],
                    **cfg['muon_propagation_config'])

    # undo renaming
    tray.AddModule('Rename', keys=['MMCTrackList', 'MMCTrackListVetoMuon'])
    tray.AddModule('Rename', keys=['_MMCTrackList', 'MMCTrackList'])

    # # now revert renaming
    # rename_dict = {
    #     'I3MCTree': t_name + '_veto_muon',
    #     'I3MCTree_preMuonProp': t_name + '_preMuonProp_veto_muon',
    #     '_' + t_name: t_name,
    #     '_' + t_name + '_preMuonProp': t_name + '_preMuonProp',
    #     '_' + t_name + '_preMuonProp_RNGState': (
    #         t_name + '_preMuonProp_RNGState',
    #     ),
    # }
    # tray.AddModule(rename_keys, 'UndoRenaming', rename_dict=rename_dict)

    # create combined I3MCTree for CLSIM
    tray.AddModule(
        CombineMCTrees, 'CombineMCTrees',
        tree1=t_name,
        tree2=t_name + 'VetoMuon',
        output_key='CombinedMuonVetoI3MCTree',
    )

    # Bias simulation if desired
    if 'ApplyBaseSimulationBias' in cfg and cfg['ApplyBaseSimulationBias']:
        tray.AddModule(
            BaseSimulationBias,
            'BaseSimulationBias',
            random_service=random_services[2],
            **cfg['BaseSimulationBiasSettings']
        )

    # --------------------------------------
    # Distance Splits
    # --------------------------------------
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
    # --------------------------------------

    click.echo('Scratch: {}'.format(scratch))
    tray.AddModule("TrashCan", "the can")
    tray.Execute()
    tray.Finish()


if __name__ == '__main__':
    main()
