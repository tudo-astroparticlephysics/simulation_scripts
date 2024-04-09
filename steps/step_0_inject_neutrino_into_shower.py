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
from icecube.dataclasses import I3Particle

from utils import create_random_services, get_run_folder
from resources import geometry
from resources.oversampling import DAQFrameMultiplier
from resources.import_events import ImportEvents

class InjectNeutrinoIntoShower(icetray.I3ConditionalModule):
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("tree_key", "The key of the I3MCTree", "I3MCTree")

    def Configure(self):
        self._tree_key = self.GetParameter("tree_key")

    def DAQ(self, frame):
        """Inject neutrino into I3MCtree.

        Parameters
        ----------
        frame : icetray.I3Frame.DAQ
        tree_key : str
        """

        tree = frame[self._tree_key]

        assert len(tree.primaries) == 1
        primary = tree.primaries[0]

        muon_energy_sum = 0
        total_leaf_energy = 0
        
        for particle in tree.get_daughters(primary.id):
            if particle.pdg_encoding in [-14, 14]:
                muon_energy_sum += particle.energy
            total_leaf_energy += particle.energy


        if muon_energy_sum > 0:
            ## Create neutrino parameters
            particles = [
                dataclasses.I3Particle.NuMu,
                dataclasses.I3Particle.NuMuBar,
                dataclasses.I3Particle.NuE,
                dataclasses.I3Particle.NuEBar,
                dataclasses.I3Particle.NuTau,
                dataclasses.I3Particle.NuTauBar,
            ]
            type = particles[np.random.randint(0, len(particles))]
            e_min = 0
            e_max = primary.energy-total_leaf_energy
            energy = np.random.uniform(e_min, e_max)
            pos = primary.pos
            dir = primary.dir
            time = primary.time

            # Create neutrino
            neutrino = I3Particle()
            neutrino.type = type
            neutrino.energy = energy
            neutrino.pos = pos
            neutrino.dir = dir
            neutrino.time = time
            neutrino.shape = dataclasses.I3Particle.StartingTrack

            # Remove every neutrino from tree
            for particle in tree:
                if particle.is_neutrino:
                    assert len(tree.get_daughters(particle)) == 0
                    tree.erase(particle)

            # Add neutrino to tree
            tree.append_child(primary.id, neutrino)

            frame.Delete(self._tree_key)
            frame.Put(self._tree_key, tree)
            frame["NuGPrimary"] = neutrino
            self.PushFrame(frame)



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

    # inject neutrino into shower
    tray.AddModule(InjectNeutrinoIntoShower, "inject_neutrino")
    
    
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
