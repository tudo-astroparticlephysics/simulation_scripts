#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/icetray-start
#METAPROJECT /data/user/mhuennefeld/software/icecube/py3-v4.1.0/combo_V01-00-00/build
from __future__ import division
import os
import click
import yaml
import time

from I3Tray import I3Tray
from icecube import icetray, dataio, simclasses

from utils import create_random_services, get_run_folder


@click.command()
@click.argument('cfg', type=click.Path(exists=True))
@click.argument('run_number', type=int)
@click.option('--scratch/--no-scratch', default=True)
def main(cfg, run_number, scratch):

    start_time = time.time()

    with open(cfg, 'r') as stream:
        cfg = yaml.full_load(stream)
    cfg['run_number'] = run_number
    cfg['run_folder'] = get_run_folder(run_number)
    if scratch:
        outfile = cfg['scratchfile_pattern'].format(**cfg)
    else:
        outfile = cfg['outfile_pattern'].format(**cfg)
    outfile = outfile.replace(' ', '0')

    infile = cfg['infile_pattern'].format(**cfg)
    infile = infile.replace(' ', '0')

    click.echo('Run: {}'.format(run_number))
    click.echo('Outfile: {}'.format(outfile))

    # -----------------------------
    # Set PPC environment variables
    # -----------------------------
    ppc_config = cfg['ppc_config']

    # define default Environment variables
    ice_base = '$I3_BUILD/ice-models/resources/models/'
    default_ppc_environment_variables = {
        'PPCTABLESDIR': ice_base + 'spice_bfr-dv1_complete',
        'OGPU': '1',  # makes sure only GPUs are used (with OpenCL version)
    }
    ppc_environment_variables = dict(default_ppc_environment_variables)
    ppc_environment_variables.update(ppc_config['environment_variables'])

    ice_model = os.path.basename(ppc_environment_variables['ICEMODELDIR'])

    # ------------------------------------------------------------------------
    # Sanity checks to ensure that the environment variables are set correctly
    # ------------------------------------------------------------------------
    if ice_model in [
        "spice_3",
        "spice_3.1",
        "spice_3.2",
        "spice_3.2.1",
        "spice_3.2.2",
        "spice_3.2.2-for_clsim",
        "spice_3.2_emrmspice_ftp-v1",
        "spice_lea",
        "spice_mie",
        "spice_bfr-v1",
        "spice_bfr-v2",
        "spice_ftp-v0",
        "spice_ftp-v1",
        "spice_ftp-v2",
        "spice_ftp-v3",
    ]:
        if "BFRM" in ppc_environment_variables:
            if ppc_environment_variables["BFRM"] != "0":
                raise ValueError(
                    "The BFRM environment variable must be set to 0 for"
                    " the older ice models."
                )

    elif ice_model in [
        "spice_ftp-v3m",
    ]:
        if 'BFRM' not in ppc_environment_variables:
            ppc_environment_variables['BFRM'] = "2"
        else:
            if ppc_environment_variables['BFRM'] != "2":
                raise ValueError(
                    "The BFRM environment variable must be set to 2 for"
                    " the spice_ftp-v3m ice model."
                )
    else:
        if 'BFRM' not in ppc_environment_variables:
            # let's be safe and throw an error for newer ice models
            # We want to make sure that the BFRM and any other potentially
            # required environment variables are set correctly.
            raise ValueError(
                "The BFRM environment variable must be set!"
            )
        else:
            # peform some sanity checks
            if "ftp-v3m" in ice_model.lower():
                if not ppc_environment_variables["BFRM"] == "2":
                    raise ValueError(
                        "The BFRM environment variable must be set to 2 for"
                        " the spice_ftp-v3m ice model."
                    )
    # ------------------------------------------------------------------------

    # define default PPC arguments
    default_ppc_arguments = {
        'MCTree': 'I3MCTree',
    }
    if 'CUDA_VISIBLE_DEVICES' in os.environ:
        default_ppc_arguments['gpu'] = int(os.environ['CUDA_VISIBLE_DEVICES'])
    ppc_arguments = dict(default_ppc_arguments)
    ppc_arguments.update(ppc_config['arguments'])

    click.echo('PPC Settings:')
    for key, value in ppc_environment_variables.items():
        click.echo('\t{}: {}'.format(key, os.path.expandvars(value)))
        os.putenv(key, os.path.expandvars(value))

    click.echo('PPC Arguments:')
    for key, value in ppc_arguments.items():
        click.echo('\t{}: {}'.format(key, value))

    # importing ppc must be done *after* setting the environment variables
    from icecube import ppc
    # ------------------------------

    # get random service
    random_services, _ = create_random_services(
        dataset_number=cfg['dataset_number'],
        run_number=cfg['run_number'],
        seed=cfg['seed'],
        n_services=1,
        use_gslrng=cfg['random_service_use_gslrng'])
    random_service = random_services[0]

    # --------------------------------------
    # Build IceTray
    # --------------------------------------
    tray = I3Tray()
    tray.AddModule('I3Reader', 'i3 reader',
                   FilenameList=[cfg['gcd_pass2'], infile])

    # run PPC
    tray.context["I3RandomService"] = random_service
    tray.AddModule("i3ppc", 'ppc', **ppc_arguments)

    # add empty MCPESeriesMap
    def add_empty_pes(frame, MCPESeriesName="MCPESeriesMap"):
        if MCPESeriesName not in frame:
            frame[MCPESeriesName] = simclasses.I3MCPESeriesMap()
    tray.Add(add_empty_pes, streams=[icetray.I3Frame.DAQ])

    # sort the pulses in time
    # (This does not seem to be the case for PPC simulations)
    def sort_pulses(frame, MCPESeriesName="MCPESeriesMap"):
        new_map = simclasses.I3MCPESeriesMap()
        for omkey, pulses in frame[MCPESeriesName]:
            new_map[omkey] = simclasses.I3MCPESeries(
                sorted(pulses, key=lambda pulse: pulse.time)
            )
        del frame[MCPESeriesName]
        frame[MCPESeriesName] = new_map
    tray.Add(sort_pulses, streams=[icetray.I3Frame.DAQ])

    # rename MCPESeriesMap to I3MCPESeriesMap
    tray.Add("Rename", keys=["MCPESeriesMap", "I3MCPESeriesMap"])

    click.echo('Output: {}'.format(outfile))
    tray.AddModule("I3Writer", "writer",
                   Filename=outfile,
                   Streams=[icetray.I3Frame.TrayInfo,
                            icetray.I3Frame.Simulation,
                            icetray.I3Frame.Stream('M'),
                            icetray.I3Frame.Stream('S'),
                            icetray.I3Frame.DAQ,
                            icetray.I3Frame.Physics])
    # --------------------------------------

    click.echo('Scratch: {}'.format(scratch))
    tray.AddModule("TrashCan", "the can")
    tray.Execute()
    tray.Finish()

    end_time = time.time()
    print("That took "+str(end_time - start_time)+" seconds.")


if __name__ == '__main__':
    main()
