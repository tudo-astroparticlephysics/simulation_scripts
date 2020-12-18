
DEFAULT_SETTINGS = {
    'global.scattering': 'NoScattering',
    'global.stopping_decay': True,
    'global.brems_multiplier': 1,
    'global.photo_multiplier': 1,
    'global.ioniz_multiplier': 1,
    'global.epair_multiplier': 1,
    'global.epair': 'EpairKelnerKokoulinPetrukhin',
    'global.brems': 'BremsKelnerKokoulinPetrukhin',
    'global.photo': 'PhotoAbramowiczLevinLevyMaor97',
    'global.lpm': True,
}

CONFIG_TEMPLATE_STRING = """{
        "global":
        {
                "seed" : 1,
                "continous_loss_output" : false,
                "only_loss_inside_detector" : true,

                "interpolation":
                {
                        "do_interpolation" : true,
                        "path_to_tables" : [
                        ],
                        "path_to_tables_readonly" : [
                                "$I3_TESTDATA/PROPOSAL/resources/tables",
                                "$I3_BUILD/PROPOSAL/resources/tables"
                        ],
                        "do_binary_tables" : false,
                        "just_use_readonly_path" : false
                },

                "exact_time" : true,
                "scattering" : {global.scattering},
                "stopping_decay" : {global.stopping_decay},

                "brems_multiplier" : {global.brems_multiplier},
                "photo_multiplier" : {global.photo_multiplier},
                "ioniz_multiplier" : {global.ioniz_multiplier},
                "epair_multiplier" : {global.epair_multiplier},
                "epair" : {global.epair},
                "brems" : {global.brems},
                "photo" : {global.photo},
                "lpm" : {global.lpm},

                "cuts_infront":
                {
                        "e_cut": -1,
                        "v_cut": 0.05,
                        "cont_rand": true
                },
                "cuts_inside":
                {
                        "e_cut": 500,
                        "v_cut": -1,
                        "cont_rand": false
                },
                "cuts_behind":
                {
                        "e_cut": -1,
                        "v_cut": -1,
                        "cont_rand": false
                }
        },

        "sectors": [
                {
                        "hierarchy": 1,
                        "medium": "air",
                        "density_correction": 0.673,

                        "geometry":
                        {
                                "shape": "sphere",
                                "origin": [0, 0, -6372186],
                                "outer_radius": 1000000000,
                                "inner_radius": 6374134
                        },
                        "cuts_inside":
                        {
                                "e_cut": 500,
                                "v_cut": -1,
                                "cont_rand": false
                        },
                        "cuts_infront":
                        {
                                "e_cut": -1,
                                "v_cut": 0.05,
                                "cont_rand": true
                        },
                        "cuts_behind":
                        {
                                "e_cut": -1,
                                "v_cut": -1,
                                "cont_rand": false
                        }
                },
                {
                        "hierarchy": 1,
                        "medium": "ice",
                        "density_correction": 0.832,

                        "geometry":
                        {
                                "shape": "sphere",
                                "origin": [0, 0, -6372186],
                                "outer_radius": 6374134,
                                "inner_radius": 6373934
                        },
                        "cuts_inside":
                        {
                                "e_cut": 500,
                                "v_cut": -1,
                                "cont_rand": false
                        },
                        "cuts_infront":
                        {
                                "e_cut": -1,
                                "v_cut": 0.05,
                                "cont_rand": true
                        },
                        "cuts_behind":
                        {
                                "e_cut": -1,
                                "v_cut": -1,
                                "cont_rand": false
                        }
                },
                {
                        "hierarchy": 1,
                        "medium": "ice",
                        "density_correction": 1.005,

                        "geometry":
                        {
                                "shape": "sphere",
                                "origin": [0, 0, -6372186],
                                "outer_radius": 6373934,
                                "inner_radius": 6371324
                        },
                        "cuts_inside":
                        {
                                "e_cut": 500,
                                "v_cut": -1,
                                "cont_rand": false
                        },
                        "cuts_infront":
                        {
                                "e_cut": -1,
                                "v_cut": 0.05,
                                "cont_rand": true
                        },
                        "cuts_behind":
                        {
                                "e_cut": -1,
                                "v_cut": -1,
                                "cont_rand": false
                        }
                },
                {
                        "hierarchy": 1,
                        "medium": "standardrock",
                        "density_correction": 1.0,

                        "geometry":
                        {
                                "shape": "sphere",
                                "origin": [0, 0, -6372186],
                                "outer_radius": 6371324,
                                "inner_radius": 0
                        },
                        "cuts_inside":
                        {
                                "e_cut": 500,
                                "v_cut": -1,
                                "cont_rand": false
                        },
                        "cuts_infront":
                        {
                                "e_cut": -1,
                                "v_cut": 0.05,
                                "cont_rand": true
                        },
                        "cuts_behind":
                        {
                                "e_cut": -1,
                                "v_cut": -1,
                                "cont_rand": false
                        }
                }
        ],

        "detector":
        {
                "shape": "cylinder",
                "origin" : [0, 0, 0],
                "outer_radius": 800,
                "inner_radius": 0,
                "height": 1600
        }
}
"""


def write_config(file_path, kwargs):
    """Write PROPOSAL config file

    Parameters
    ----------
    file_path : str
        The path to which the config will be written to.
    **kwargs
        Keyword arguments that will be used to replace placeholders
        in the CONFIG_TEMPLATE_STRING.
    """
    settings = dict(DEFAULT_SETTINGS)
    settings.update(kwargs)

    config_str = str(CONFIG_TEMPLATE_STRING).replace('        ', '\t')
    for key, value in settings.items():

        # sanity and safety check to let user know if a key does not exist
        if '{' + key + '}' not in config_str:
            raise KeyError('Key {} does not exist!'.format(key))

        if isinstance(value, str):
            config_str = config_str.replace('{' + key + '}', '"' + value + '"')
        elif isinstance(value, bool):
            if value:
                config_str = config_str.replace('{' + key + '}', 'true')
            else:
                config_str = config_str.replace('{' + key + '}', 'false')
        else:
            config_str = config_str.replace('{' + key + '}', str(value))

    with open(file_path, "w") as text_file:
        text_file.write(config_str)
