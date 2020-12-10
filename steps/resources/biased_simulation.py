from __future__ import division
import numpy as np

from I3Tray import I3Tray
from icecube import icetray, dataclasses

from .bias_utils import bias_functions


class BaseSimulationBias(icetray.I3ConditionalModule):

    """Class to bias simulation based on their I3MCTree.

    Depending on the specified bias function, this may be intended to run
    after PROPOSAL.

    Attributes
    ----------
    mctree_name : str
        The name of the I3MCTree key.
    """

    def __init__(self, context):
        """Class to import events from another I3-File

        Parameters
        ----------
        context : TYPE
            Description
        """
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddOutBox('OutBox')
        self.AddParameter(
            'bias_function',
            'The bias function class. This function computes and returns the '
            'probability that the given event is kept plus a dictionary which '
            'holds any addtional infor that will be written to the bias info. '
            'This probability should always be in (0, 1].'
            'The input to this function is a dictionary containing '
            'the data for the bias function. If the `bias_function` '
            'is provided as a str, the bias function class will be loaded '
            'from .bias_utils.bias_functions')
        self.AddParameter(
            'bias_function_kwargs',
            'This is only relevant if `bias_function` is provided as a str. '
            'In this case, the `bias_function_kwargs` will be provided as '
            'keyword arguments to the function class initializer. '
            '`bias_function_kwargs` must be a dictionary.',
            {})
        self.AddParameter(
            'lower_probability_bound',
            'A lower bound of this value is applied to the computed keep '
            'probability.',
            1e-6)
        self.AddParameter(
            'keep_all_events',
            'If True, all events are kept and the bias results '
            'are only written to the frame',
            False)
        self.AddParameter(
            'verbose_output',
            'If True, additional bias info is written to the '
            'output key.',
            True)
        self.AddParameter(
            'output_key',
            'The output base key to which bias weights will be saved.',
            'BiasedSimulationWeight')
        self.AddParameter(
            'random_service',
            'The random service or seed to use. If this is an '
            'integer, a numpy random state will be created with '
            'the seed set to `random_service`',
            42)
        self.AddParameter(
            'verbose', 'If True, more verbose output is provided.', False)

    def Configure(self):
        """Configures MuonLossProfileFilter.
        """
        self.mctree_name = self.GetParameter('mctree_name')
        self.bias_function = self.GetParameter('bias_function')
        self.bias_function_kwargs = self.GetParameter('bias_function_kwargs')
        self.lower_probability_bound = self.GetParameter(
            'lower_probability_bound')
        self.keep_all_events = self.GetParameter('keep_all_events')
        self.verbose_output = self.GetParameter('verbose_output')
        self.output_key = self.GetParameter('output_key')
        self.random_service = self.GetParameter('random_service')
        self.verbose = self.GetParameter('verbose')

        if isinstance(self.random_service, int):
            self.random_service = np.random.RandomState(self.random_service)

        # set up bias function class
        if isinstance(self.bias_function, str):
            self.bias_function = getattr(bias_functions, self.bias_function)(
                **self.bias_function_kwargs)

    def DAQ(self, frame):
        """Bias events based on the specified bias function.
        """

        # start timer
        t_0 = timeit.default_timer()

        # compute keep probability
        bias_data = {
            'frame': frame,
        }
        keep_prob, additional_bias_info = self.bias_function(bias_data)

        keep_prob = float(keep_prob)
        assert keep_prob > 0. and keep_prob <= 1., keep_prob
        keep_prob = np.clip(keep_prob, self.lower_probability_bound, 1.0)

        passed = self.random_service.uniform() <= keep_prob

        bias_weights = dataclasses.I3MapStringDouble({
            'keep_prob': keep_prob,
            'weight_multiplier': 1. / keep_prob,
            'passed': float(passed),
        })

        # stop timer
        t_1 = timeit.default_timer()

        # add verbose output if desired
        if self.verbose_output:
            bias_weights['runtime_total'] = t_1 - t_0
            for key, value in additional_bias_info.items():
                bias_weights[key] = float(value)

        if self.verbose:
            print('Biasing took: {:3.3f}ms'.format((t_1 - t_0) * 1000))

        frame[self.output_key] = bias_weights

        # push frame to next modules
        if self.keep_all_events:
            self.PushFrame(frame)
        else:
            if passed:
                self.PushFrame(frame)
