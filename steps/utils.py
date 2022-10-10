import numpy as np
import importlib

MAX_DATASET_NUMBER = 100000
MAX_RUN_NUMBER = 100000


def create_random_services_settings(
        dataset_number, run_number, seed, n_services=1, use_gslrng=False):

    if run_number < 0:
        raise RuntimeError("negative run numbers are not supported")
    elif run_number >= MAX_RUN_NUMBER:
        raise RuntimeError("run numbers > %u are not supported".format(
            MAX_RUN_NUMBER))

    if dataset_number < 0:
        raise RuntimeError("negative dataset numbers are not supported")

    max_run_num = MAX_RUN_NUMBER // 10

    int_run_number = dataset_number * max_run_num + run_number

    settings_list = []
    for i in range(n_services):
        streamnum = run_number + (MAX_RUN_NUMBER * i)

        if use_gslrng:
            if MAX_RUN_NUMBER * n_services > (2 ** 31 - 1 - streamnum) / seed:
                raise ValueError(
                    'Integer 32bit overflow encountered in seed generation. '
                    'This can cause issues with CORSIKA simulations. '
                    + 'Combined GSLRNG seed: {} | '.format(
                        seed*MAX_RUN_NUMBER*n_services + streamnum)
                    + 'MAX_RUN_NUMBER: {} | n_services: {} | seed: {}'.format(
                        MAX_RUN_NUMBER, n_services, seed)
                )
            settings_list.append(dict(
                seed=seed*MAX_RUN_NUMBER*n_services + streamnum))
        else:
            settings_list.append(dict(
                seed=seed,
                nstreams=MAX_RUN_NUMBER * n_services,
                streamnum=streamnum))
    return settings_list, int_run_number


def create_random_services(dataset_number, run_number, seed, n_services=1,
                           use_gslrng=False):
    from icecube import phys_services, icetray, dataclasses

    # get seed settings
    settings_list, int_run_number = create_random_services_settings(
        dataset_number=dataset_number,
        run_number=run_number,
        seed=seed,
        n_services=n_services,
        use_gslrng=use_gslrng,
    )

    random_services = []
    for i in range(n_services):

        if use_gslrng:
            random_services.append(phys_services.I3GSLRandomService(
                **settings_list[i]))
        else:
            random_services.append(phys_services.I3SPRNGRandomService(
                **settings_list[i]))

    return random_services, int_run_number


def get_run_folder(run_number, runs_per_folder=1000):
    fill = int(np.log10(MAX_RUN_NUMBER) + 0.5)
    start = (run_number // runs_per_folder) * runs_per_folder
    stop = start + runs_per_folder - 1
    return '{}-{}'.format(str(start).zfill(fill), str(stop).zfill(fill))


def load_class(full_class_string):
    """
    dynamically load a class from a string

    Parameters
    ----------
    full_class_string : str
        The full class string to the given python clas.
        Example:
            my_project.my_module.my_class

    Returns
    -------
    python class
        PYthon class defined by the 'full_class_string'
    """

    class_data = full_class_string.split(".")
    module_path = ".".join(class_data[:-1])
    class_str = class_data[-1]

    module = importlib.import_module(module_path)
    # Finally, we retrieve the Class
    return getattr(module, class_str)


muongun_keys = ['MCOversizeStreamDefault',
                'MCOversizeStream0',
                'MCOversizeStream1',
                'MCOversizeStream2',
                'MCOversizeStream3',
                'MCMuon',
                'MuonEffectiveArea',
                'MCOversizing',
                'GenerateCosmicRayMuons',
                'MCDomThresholds',
                'MCDistanceCuts']
