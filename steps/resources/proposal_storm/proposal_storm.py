from icecube import icetray, dataclasses

from .config import write_config


class PROPOSALStorm(icetray.I3Module):
    """
    Add a "m" frame for the PROPOSAL settings.
    """
    def __init__(self, ctx):
        super(PROPOSALStorm, self).__init__(ctx)
        self.AddParameter("ConfigFilePath", "Path to output config file.")
        self.AddParameter("RandomService", "RandomService to use.")
        self.AddParameter(
            "UniformRanges",
            "Dictionary of settings which will be sampled uniformly "
            "in the given range. Format: {'name': [min, max]}", {})
        self.AddParameter(
            "DiscreteOptions",
            "Dictionary of settings which will be sampled uniformly "
            "from the given options. Format: {'name': [option1, option2]}", {})
        self.AddParameter(
            "OutputKey", "Output key of the meta data.", 'PROPOSALStorm')

    def Configure(self):
        self.config_file_path = self.GetParameter("ConfigFilePath")
        self._rnd = self.GetParameter("RandomService")
        self._ranges_dict = self.GetParameter("UniformRanges")
        self._options_dict = self.GetParameter("DiscreteOptions")
        self._output_key = self.GetParameter("OutputKey")

        if self._options_dict != {}:
            raise NotImplementedError('Discrete options not yet supported!')

        self._frame_has_been_pushed = False

        # sample setting values
        self.sampled_settings = self.sample_settings()

        # write config file with these settings
        write_config(self.config_file_path, self.sampled_settings)

    def Process(self):

        # get next frame
        frame = self.PopFrame()

        if not self._frame_has_been_pushed:

            # create settings frame and push it
            settings_frame = icetray.I3Frame('m')

            settings_data = {}

            # add meta data on ranges
            for key, value_range in self._ranges_dict.items():
                settings_data[key+'RangeMin'] = value_range[0]
                settings_data[key+'RangeMax'] = value_range[1]

            settings_data.update(self.sampled_settings)
            settings_frame[self._output_key] = dataclasses.I3MapStringDouble(
                settings_data)

            self.PushFrame(settings_frame)

            self._frame_has_been_pushed = True

        self.PushFrame(frame)

    def sample_settings(self):
        """Sample PROPOSAL Settings
        """
        settings_dict = {}

        # Uniform ranges
        for key, value_range in self._ranges_dict.items():
            settings_dict[key] = float(self._rnd.uniform(*value_range))

        # Uniform options
        for key, options in self._options_dict.items():
            raise NotImplementedError('Discrete options not yet supported!')

        return settings_dict
