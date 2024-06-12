#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/icetray-start
#METAPROJECT icetray/v1.10.0

from __future__ import division

import sys

sys.path.append("/data/user/lbollmann/repositories/ic3-labels")  # >:(

import click
import yaml
import numpy as np
import glob

from icecube.simprod import segments
from icecube.sim_services.propagation import get_propagators
from icecube.icetray import I3Tray
from icecube import icetray, dataclasses, neutrino_generator, simclasses, dataio
from icecube.dataclasses import I3Particle

# print("ASDASDASD")
# import os
# print(os.environ["PYTHONPATH"])
# print(sys.path)

from ic3_labels.labels.utils import muon as mu_utils
from ic3_labels.labels.utils import high_level as hl
from ic3_labels.labels.utils.detector import icecube_hull

from utils import create_random_services, get_run_folder


# Transform uniform distribution to power law distribution
def p(x, x0, x1, n):
    return (x1**(n+1)-x0**(n+1)*x+x0**(n+1))**(1/(n+1))

class InjectNeutrinoIntoShower(icetray.I3ConditionalModule):
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("TreeKey", "The key of the I3MCTree", "I3MCTree")
        self.AddParameter(
            "MaxNeutrinoEnergy",
            "The maximum energy that injected neutrinos can have",
            10e6,
        )
        self.AddParameter(
            "MinNeutrinoEnergy",
            "The minimum energy that injected neutrinos can have",
            1e3,
        )
        self.AddParameter(
            "MinEnergyRatio",
            "The minimum energy ratio between neutrino and muon bundle",
            1.0,
        )
        self.AddParameter(
            "TypeRatios",
            """The ratios of all neutrino types. Must be an array with length 6 with values between 0 and 1.
        The order is [NuE, NuEBar, NuMu, NuMuBar, NuTau, NuTauBar]. Default is [1, 1, 1, 1, 1, 1].""",
            [1, 1, 1, 1, 1, 1],
        )
        self.AddParameter("RandomService", "The random service", None)
        self.AddParameter("MinMuonBundleEnergy", "The minimum energy of the muon bundle. If lower the frame is dropped.", 1e3)
        self.AddParameter("SpectralIndex", "The spectral index of the power law distribution", -1)


    def Configure(self):
        self._tree_key = self.GetParameter("TreeKey")
        self._max_neutrino_energy = self.GetParameter("MaxNeutrinoEnergy")
        self._min_neutrino_energy = self.GetParameter("MinNeutrinoEnergy")
        self._min_energy_ratio = self.GetParameter("MinEnergyRatio")
        self._ratios = self.GetParameter("TypeRatios")
        self._random_service = self.GetParameter("RandomService")
        self._min_muon_bundle_energy = self.GetParameter("MinMuonBundleEnergy")
        self._spectral_index = self.GetParameter("SpectralIndex")

    def DAQ(self, frame):
        """Inject neutrino into I3MCtree.

        Parameters
        ----------
        frame : icetray.I3Frame.DAQ
        tree_key : str
            The key of the I3MCTree
        max_neutrino_energy : float
            The maximum energy that injected neutrinos can have
        min_energy_ratio : float
            The minimum energy ratio between neutrino and muon bundle
        """

        track_cache, _ = mu_utils.get_muongun_track_cache(frame)

        labels = hl.get_muon_bundle_information(
            frame=frame,
            convex_hull=icecube_hull,
            track_cache=track_cache,
        )

        tree = frame[self._tree_key]

        assert len(tree.primaries) == 1
        primary = tree.primaries[0]

        bundle_energy_at_entry = labels["bundle_energy_at_entry"]

        # If the muon bundle energy is too high or zero then drop frame
        if (
            bundle_energy_at_entry <= self._min_muon_bundle_energy
            or bundle_energy_at_entry * self._min_energy_ratio
            >= self._max_neutrino_energy
        ):
            return

        ## Create neutrino parameters
        particles = [
            dataclasses.I3Particle.NuE,
            dataclasses.I3Particle.NuEBar,
            dataclasses.I3Particle.NuMu,
            dataclasses.I3Particle.NuMuBar,
            dataclasses.I3Particle.NuTau,
            dataclasses.I3Particle.NuTauBar,
        ]

        if self._random_service is None:
            raise ValueError(
                "random_service parameter of InjectNeutrinoIntoShower is not set"
            )

        # Sample type according to rations
        cdf = np.cumsum(self._ratios)
        num = self._random_service.uniform(float(cdf[-1]))

        nu_type = particles[np.searchsorted(cdf, num)]

        e_min = max(
            bundle_energy_at_entry * self._min_energy_ratio, self._min_neutrino_energy
        )
        e_max = self._max_neutrino_energy
        

        if self._spectral_index != -1:
            uniform = self._random_service.uniform(0, 1)
            energy = p(uniform, e_min, e_max, self._spectral_index)
        else:
            energy = self._random_service.uniform(e_min, e_max)

        pos = primary.pos
        dir = primary.dir
        time = primary.time

        # Create neutrino
        neutrino = I3Particle()
        neutrino.type = nu_type
        neutrino.energy = energy
        neutrino.pos = pos
        neutrino.dir = dir
        neutrino.time = time
        neutrino.shape = dataclasses.I3Particle.StartingTrack

        # Add neutrino to tree
        tree_pre_muon_prop = frame["I3MCTree_preMuonProp"]
        tree_pre_muon_prop.append_child(primary.id, neutrino)
        frame.Delete("I3MCTree_preMuonProp")
        frame.Put("I3MCTree_preMuonProp", tree_pre_muon_prop)

        frame["NuGPrimary"] = neutrino

        labels = dataclasses.I3MapStringDouble()

        labels["EnergyMin"] = e_min
        labels["EnergyMax"] = e_max
        labels["Energy"] = neutrino.energy
        labels["EnergyRatio"] = neutrino.energy / bundle_energy_at_entry
        labels["MinEnergyRatio"] = self._min_energy_ratio
        labels["MinNeutrinoEnergy"] = self._min_neutrino_energy
        labels["BundleEnergyAtEntry"] = bundle_energy_at_entry

        frame["InjectNeutrinoInfo"] = labels
        frame["InjectNeutrinoRatios"] = dataclasses.I3VectorDouble(self._ratios)
        self.PushFrame(frame)


@click.command()
@click.argument("cfg", type=click.Path(exists=True))
@click.argument("run_number", type=int)
@click.option("--scratch/--no-scratch", default=True)
def main(cfg, run_number, scratch):
    with open(cfg, "r") as stream:
        cfg = yaml.full_load(stream)
    cfg["run_number"] = run_number
    cfg["run_folder"] = get_run_folder(run_number)
    if scratch:
        outfile = cfg["scratchfile_pattern"].format(**cfg)
    else:
        outfile = cfg["outfile_pattern"].format(**cfg)
    outfile = outfile.replace(" ", "0")

    # ------------------------------
    # get list of files for this run
    # ------------------------------
    import_cfg = cfg["event_import_settings"]
    glob_files = import_cfg["input_file_glob_list"]
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

    click.echo("Run: {}".format(run_number))
    click.echo("Outfile: {}".format(outfile))
    click.echo("input Files:")
    for file in files:
        click.echo("\t{}".format(file))

    # crate random services
    if "random_service_use_gslrng" not in cfg:
        cfg["random_service_use_gslrng"] = False
    random_services, _ = create_random_services(
        dataset_number=cfg["dataset_number"],
        run_number=cfg["run_number"],
        seed=cfg["seed"],
        n_services=2,
        use_gslrng=cfg["random_service_use_gslrng"],
    )

    # --------------------------------------
    # Build IceTray
    # --------------------------------------
    tray = I3Tray()
    tray.context["I3RandomService"] = random_services[0]

    # import events from another I3-file
    tray.AddModule("I3Reader", "reader", FilenameList=files)

    if "muon_propagation_config" in cfg:
        tray.Add("Rename", Keys=["I3MCTree", "I3MCTree_preMuonProp"])
        tray.AddSegment(
            segments.PropagateMuons,
            "propagate_muons",
            RandomService=random_services[1],
            **cfg["muon_propagation_config"],
        )

    # inject neutrino into shower
    tray.AddModule(
        InjectNeutrinoIntoShower,
        "inject_neutrino",
        RandomService=random_services[0],
        **cfg["inject_neutrino_cfg"],
    )

    ## Delete propagated tree and MMCTracklist
    tray.Add("Delete", "delete_tree", Keys=["I3MCTree"])
    tray.Add("Rename", "rename_pre_prop_tree", Keys=["I3MCTree_preMuonProp", "I3MCTree"])
    tray.Add("Delete", "delete_mmctracklist", Keys=["MMCTrackList"])

    click.echo("Output: {}".format(outfile))
    tray.AddModule(
        "I3Writer",
        "writer",
        Filename=outfile,
        Streams=[
            icetray.I3Frame.DAQ,
            icetray.I3Frame.Physics,
            icetray.I3Frame.Stream("S"),
            icetray.I3Frame.Stream("M"),
        ],
    )
    # --------------------------------------

    click.echo("Scratch: {}".format(scratch))
    tray.AddModule("TrashCan", "the can")
    tray.Execute()
    tray.Finish()


if __name__ == "__main__":
    main()
