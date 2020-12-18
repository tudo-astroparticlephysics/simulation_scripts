import os
import json


def write_config(file_path, setting_updates, default_config=os.path.expandvars(
        "$I3_BUILD/PROPOSAL/resources/config_icesim.json")):
    """Write PROPOSAL config file

    Parameters
    ----------
    file_path : str
        The path to which the config will be written to.
    setting_updates : dict
        Updates to the default PROPOSAL config as specicied in
        `default_config`. The updates must be provided as dictionary of
        key: value pairs, where value is the updated value and key is the
        name of the setting to update. Nested settings may be accessed by
        adding '.' in the key.
    default_config : str, optional
        The file path to the default PROPOSAL config.
    """
    with open(default_config) as json_file:
        data = json.load(json_file)

    # remove writable tables path and disable readonly:
    # allows PROPOSAL to build on-the-fly interpolations in memory
    data['global']['interpolation']['path_to_tables'] = []
    data['global']['interpolation']['just_use_readonly_path'] = False

    # now go through settings and update values
    for key, value in setting_updates.items():
        layers = key.split('.')

        dict_element = data
        for layer in layers[:-1]:
            dict_element = data[layer]

        dict_element[layers[-1]] = value

    with open(file_path, "w") as text_file:
        json.dump(data, text_file, indent=4, sort_keys=True)
