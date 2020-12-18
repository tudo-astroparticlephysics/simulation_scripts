from icecube import icetray, dataclasses

from .config import write_config


class PROPOSALStorm:

    def __init__(
            self, config_file_path, random_service,
            uniform_ranges={}, discrete_options={},
            ):
        """Set up PROPOSAL Storm

        Parameters
        ----------
        config_file_path : str
            Path to which the config file will be written to.
        random_service : obj
            A random number service.
        uniform_ranges : dict, optional
            A dictionary of settings which will be sampled uniformly in the
            provided range. Format:
                {'name': [min, manx]}
        discrete_options : dict, optional
            A dictionary of settings which will be sampled uniformly from
            the provided options. Format:
                {'name': [option1, option2]}
            [Note: this option is not yet supported!]
        """
        self.config_file_path = config_file_path
        self._rnd = random_service
        self.uniform_ranges = uniform_ranges
        self.discrete_options = discrete_options

        if self.discrete_options != {}:
            raise NotImplementedError('Discrete options not yet supported!')

        # sample setting values
        self.sampled_settings = self.sample_settings()

        # write config file with these settings
        write_config(self.config_file_path, self.sampled_settings)

    def sample_settings(self):
        """Sample PROPOSAL Settings
        """
        settings_dict = {}

        # Uniform ranges
        for key, value_range in self.uniform_ranges.items():
            settings_dict[key] = float(self._rnd.uniform(*value_range))

        # Uniform options
        for key, options in self.discrete_options.items():
            raise NotImplementedError('Discrete options not yet supported!')

        return settings_dict


class PROPOSALStormWriter(icetray.I3Module):
    """
    Add a "m" frame for the PROPOSAL settings.
    """
    def __init__(self, ctx):
        super(PROPOSALStormWriter, self).__init__(ctx)
        self.AddParameter("PROPOSALStormObject", "PROPOSALStorm object.")
        self.AddParameter(
            "OutputKey", "Output key for the PROPOSAL-Storm data.",
            'PROPOSALStorm')

    def Configure(self):
        self.proposal_storm = self.GetParameter("PROPOSALStormObject")
        self._output_key = self.GetParameter("OutputKey")
        self._frame_has_been_pushed = False

    def Process(self):

        # get next frame
        frame = self.PopFrame()

        if not self._frame_has_been_pushed:

            # create settings frame and push it
            settings_frame = icetray.I3Frame('m')

            # add meta data on ranges
            settings_data = {}
            for key, value_range in self.proposal_storm.uniform_ranges.items():
                settings_data[key+'RangeMin'] = value_range[0]
                settings_data[key+'RangeMax'] = value_range[1]

            # write to frame
            settings_frame[self._output_key+'UniformRanges'] = (
                dataclasses.I3MapStringDouble(settings_data)
            )
            settings_frame[self._output_key] = dataclasses.I3MapStringDouble(
                self.proposal_storm.sampled_settings)

            self.PushFrame(settings_frame)

            self._frame_has_been_pushed = True

        self.PushFrame(frame)
