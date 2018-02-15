#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT /home/mmeier/combo_stable/build
from __future__ import division
import click
import yaml
import copy

import numpy as np

from I3Tray import I3Tray, I3Units
from icecube import icetray, dataclasses
from icecube.simprod import segments

from utils import create_random_services, get_run_folder

import os
import sys

from icecube import VHESelfVeto


def scramble_azimuth(mctree, random_state):
    # For each particle in the MCTree, rotate the x-y coordinates
    # by an azimuth shift angle
    if not isinstance(random_state, np.random.RandomState):
        random_state = np.random.RandomState(random_state)

    azi_shift = np.deg2rad(random_state.uniform(0., 360.))

    rotation_matrix = np.array([[np.cos(azi_shift), -np.sin(azi_shift)],
                                [np.sin(azi_shift), np.cos(azi_shift)]])

    for i, particle in enumerate(mctree):
        pos = np.array([particle.pos.x,
                        particle.pos.y])
        rotated_pos = pos.dot(rotation_matrix)
        mctree[i].pos.x = rotated_pos[0]
        mctree[i].pos.y = rotated_pos[1]

        dir_ = np.array([particle.dir.x,
                         particle.dir.y])
        dir_z = particle.dir.z
        rotated_dir = dir_.dot(rotation_matrix)
        new_dir = dataclasses.I3Direction(
            rotated_dir[0], rotated_dir[1], dir_z)
        mctree[i].dir = new_dir

    return mctree


class MCTreeStripper(icetray.I3ConditionalModule):
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)

    def Configure(self):
        pass

    def DAQ(self, frame):
        mctree = frame['I3MCTree']
        new_tree = dataclasses.I3LinearizedMCTree()
        primaries = mctree.get_primaries()
        for primary in primaries:
            new_tree.add_primary(primary)
            daughters = mctree.get_daughters(primary)
            new_tree.append_children(primary, daughters)

        del frame['I3MCTree']
        frame['I3MCTree'] = new_tree

        frame['I3EventHeader_Original'] = copy.deepcopy(frame['I3EventHeader'])
        del frame['I3EventHeader']

        self.PushFrame(frame)

    def Physics(self, frame):
        pass


class QFactory(icetray.I3ConditionalModule):
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter('n_events_per_event', '', 10)
        self.AddParameter('random_state', '', 1337)

    def Configure(self):
        self.n_events = self.GetParameter('n_events_per_event')
        self.random_state = self.GetParameter('random_state')

    def DAQ(self, frame):
        for i in range(self.n_events):
            new_frame = copy.deepcopy(frame)
            new_frame['OversamplingIndex'] = icetray.I3Int(i)
            new_frame['CorsikaWeightMap']['OverSampling'] = \
                float(self.n_events)
            mctree = new_frame['I3MCTree']
            new_mctree = scramble_azimuth(mctree, self.random_state)
            del new_frame['I3MCTree']
            new_frame['I3MCTree'] = new_mctree
            self.PushFrame(new_frame)


@click.command()
@click.argument('cfg', click.Path(exists=True))
@click.argument('run_number', type=int)
@click.option('--scratch/--no-scratch', default=True)
def main(cfg, run_number, scratch):
    with open(cfg, 'r') as stream:
        cfg = yaml.load(stream)
    cfg['run_number'] = run_number
    cfg['run_folder'] = get_run_folder(run_number)
    if scratch:
        outfile = cfg['scratchfile_pattern'].format(**cfg)
    else:
        outfile = cfg['outfile_pattern'].format(**cfg)
    outfile = outfile.replace(' ', '0')

    print(cfg['input_pattern'])
    infile_cfg = {}
    infile_cfg['run_folder'] = cfg['run_folder']
    infile_cfg['run_number'] = str(cfg['run_number']).zfill(6)
    infile = cfg['input_pattern'].format(**infile_cfg)
    infile = infile.replace(' ', '0')

    if not os.path.isfile(infile):
        print(infile)
        print('File does not exist! Skipping it...')
        sys.exit(0)

    click.echo('Run: {}'.format(run_number))
    click.echo('Outfile: {}'.format(outfile))
    click.echo('n_events_per_event: {}'.format(cfg['n_events_per_event']))

    tray = I3Tray()
    random_services, _ = create_random_services(
        dataset_number=cfg['dataset_number'],
        run_number=cfg['run_number'],
        seed=cfg['seed'],
        n_services=2)

    # Read an input file
    tray.AddModule('I3Reader', 'read l2 corsika file',
                   Filenamelist=[cfg['gcd'], infile])

    # Apply a Qtot cut
    tray.AddModule('HomogenizedQTot', 'Qtot calculation',
                   Pulses=cfg['pulses'])

    tray.AddModule('VHESelfVeto', 'selfveto',
                   TimeWindow=3000,
                   VertexThreshold=250,
                   DustLayer=-160,
                   DustLayerWidth=60,
                   VetoThreshold=3,
                   Pulses=cfg['pulses'])

    tray.AddModule('HomogenizedQTot', 'Causal Qtot',
                   Pulses=cfg['pulses'],
                   Output='CausalQTot',
                   VertexTime='VHESelfVetoVertexTime')

    tray.AddModule(lambda frame: np.log10(frame['CausalQTot'].value) >= 3.)

    # Strip down the MCTree and inject n_events_per_event new Q Frames
    tray.AddModule(MCTreeStripper, 'Strip down the MCTree')

    keep_before_resim = [
        'I3EventHeader_Original',
        'I3MCTree',
        'CorsikaWeightMap'
    ]

    tray.AddModule("Keep", "keep_before_merge",
                   keys=keep_before_resim)

    tray.AddModule(QFactory, 'produce some q frames',
                   n_events_per_event=cfg['n_events_per_event'],
                   random_state=cfg['seed'] * cfg['run_number'])

    class convert_mctree(icetray.I3Module):
        def __init__(self, context):
            icetray.I3Module.__init__(self, context)
            self.AddOutBox('OutBox')

        def Configure(self):
            pass

        def DAQ(self, frame):
            mctree = frame['I3MCTree']
            new_mctree = dataclasses.I3LinearizedMCTree(mctree)
            del frame['I3MCTree']
            frame['I3MCTree'] = new_mctree
            self.PushFrame(frame)

    # Propagate Muons
    tray.AddSegment(segments.PropagateMuons,
                    'propagate_muons',
                    RandomService=random_services[1],
                    InputMCTreeName='I3MCTree')

    # Proposal likes to convert Linearized MCTrees to ordinary ones
    tray.AddModule(convert_mctree, 'convert')

    class KeyValueTester():
        def __init__(self, key, value):
            self.key = key
            self.value = value

        def __call__(self, frame):
            if frame.Stop == icetray.I3Frame.DAQ:
                if frame[self.key].value == self.value:
                    return True
                else:
                    return False
            else:
                return True

        def __str__(self):
            s = 'KeyValueTestObject'
            return s

    if cfg['write_multiple_files']:
        for i in range(cfg['n_events_per_event']):
            outfile_i = outfile.replace('.i3.bz2', '_{}.i3.bz2'.format(i))
            kv_tester_i = KeyValueTester('OversamplingIndex', i)
            tray.AddModule(
                'I3Writer', 'writer_{}'.format(i),
                Filename=outfile_i,
                Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Stream('M')],
                If=kv_tester_i)
    else:
        tray.AddModule(
            'I3Writer', 'write',
            Filename=outfile,
            Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Stream('M')])
    tray.AddModule('TrashCan', 'trash')
    tray.Execute()
    tray.Finish()
    del tray

if __name__ == '__main__':
    main()
