from icecube import icetray, dataclasses, simclasses
from I3Tray import I3Units
from icecube.filterscripts import filter_globals
from icecube.filterscripts.baseproc import BaseProcessing
from icecube.icetray.i3logging import log_error, log_fatal
from icecube.STTools.seededRT.configuration_services import \
    I3DOMLinkSeededRTConfigurationService
from icecube import filter_tools
import numpy as np


@icetray.traysegment
def GetPulses(tray, name,
              simulation=False,
              decode=False,
              sdstarchive=False,
              slop_split_enabled=True,
              needs_wavedeform_spe_corr=False,
              If=lambda f: True,
              ):
    '''
    Relevant part of OnlineL2 tray segment that creates pulses from
    InIceRawData. Taken from:
        https://code.icecube.wisc.edu/projects/icecube/browser/IceCube/
        projects/filterscripts/trunk/python/all_filters.py
    '''

    # Create a SeededRT configuration object with the standard RT settings.
    # This object will be used by all the different SeededRT modules, i.e. the
    # modules use the same causial space and time conditions, but can use
    # different seed algorithms.
    seededRTConfig = I3DOMLinkSeededRTConfigurationService(
                         ic_ic_RTRadius              = 150.0*I3Units.m,
                         ic_ic_RTTime                = 1000.0*I3Units.ns,
                         treat_string_36_as_deepcore = False,
                         useDustlayerCorrection      = False,
                         allowSelfCoincidence        = True
                     )
    # base processing requires:  GCD and frames being fed in by reader or Inlet
    # base processing include:
    #     decoding, TriggerCheck, Bad Dom cleaning, calibrators,
    #     Feature extraction, pulse cleaning (seeded RT, and Time Window),
    #     PoleMuonLineit, PoleMuonLlh, Cuts module and Mue on PoleMuonLlh
    if sdstarchive:
        tray.AddSegment(BaseProcessing, "BaseProc",
                        pulses=filter_globals.CleanedMuonPulses,
                        decode=decode,
                        simulation=False,
                        needs_calibration=False, needs_superdst=False,
                        do_slop=slop_split_enabled,
                        needs_trimmer=False, seededRTConfig=seededRTConfig
                        )
    else:
        tray.AddSegment(BaseProcessing, "BaseProc",
                        pulses=filter_globals.CleanedMuonPulses,
                        decode=decode,
                        simulation=simulation,
                        do_slop=slop_split_enabled,
                        seededRTConfig=seededRTConfig,
                        needs_wavedeform_spe_corr=needs_wavedeform_spe_corr
                        )


class MergeOversampledEvents(icetray.I3ConditionalModule):

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter('OversamplingFactor', 'Oversampling factor.', None)
        self.AddParameter('PulseKey', 'The pulse key over which to aggregate.',
                          'InIceDSTPulses')
        self.AddParameter('MinPulseTimeSeparation',
                          'The minimum time in ns a merged pulse may be '
                          'separated from a previous one. If it is less, then '
                          'the two pulses will be merged into a single one',
                          1.)

    def Configure(self):
        self.oversampling_factor = self.GetParameter('OversamplingFactor')
        self.pulse_key = self.GetParameter('PulseKey')
        self.min_separation = self.GetParameter('MinPulseTimeSeparation')
        self.current_time_shift = None
        self.current_event_counter = None
        self.current_aggregation_frame = None
        self.current_daq_frame = None
        self.oversampling_counter = None
        self.pushed_frame_already = False

    def push_aggregated_frame(self):

        # adjust charges of pulses
        for om_key in self.merged_pulse_series.keys():
            for pulse in self.merged_pulse_series[om_key]:
                pulse.charge /= self.oversampling_counter

        self.current_aggregation_frame['AggregatedPulses'] = \
            self.merged_pulse_series

        # update oversampling dictionary
        dic = dict(self.current_aggregation_frame['oversampling'])
        del self.current_aggregation_frame['oversampling']
        dic['num_aggregated_pulses'] = self.oversampling_counter
        dic['time_shift'] = self.current_time_shift
        self.current_aggregation_frame['oversampling'] = \
            dataclasses.I3MapStringDouble(dic)

        self.PushFrame(self.current_daq_frame)
        self.PushFrame(self.current_aggregation_frame)
        self.current_daq_frame = None
        self.pushed_frame_already = True

    def merge_pulse_series(self, pulse_series, new_pulses, time_shift):
        """Merge two pulse series.

        Assumes that new_pulses are to be merged into existing pulse_series,
        e.g. new_pulses are smaller than pulse_series

        Parameters
        ----------
        pulse_series : dataclasses.I3RecoPulseSeriesMap
            Pulse series map to which the new pulses will be added
        new_pulses : dataclasses.I3RecoPulseSeriesMap
            New pulse series that will be merged into the existing pulse
            series.
        time_shift : float
            The time shift of the new pulses.

        Returns
        -------
        dataclasses.I3RecoPulseSeriesMap
            Merged pulse series map.
        """
        pulse_series = dataclasses.I3RecoPulseSeriesMap(pulse_series)

        # calculate relative time difference to first oversampling frame
        delta_t = time_shift - self.current_time_shift

        # Go through every key in new pulses:
        for key, new_hits in new_pulses:
            if key not in pulse_series:
                # DOM has not been previously hit: just merge all hits
                pulse_series[key] = new_hits

                # correct times:
                for new_hit in pulse_series[key]:
                    new_hit.time += delta_t
            else:
                # DOM already has hits:
                #   now need to merge new pulses in existing series
                # Loop through existing pulses and sort them in
                merged_hits = list(pulse_series[key])
                len_merged_hits = len(merged_hits)
                index = 0
                for new_hit in new_hits:
                    pulse_is_merged = False
                    combine_pulses = False

                    # correct for relative time shift difference
                    new_hit.time += delta_t

                    # sort the pulse into existing list
                    while not pulse_is_merged:
                        if (index >= len(merged_hits) or
                                new_hit.time < merged_hits[index].time):

                            time_diff = abs(
                                new_hit.time - merged_hits[index - 1].time)
                            combine_pulses = time_diff < self.min_separation
                            if combine_pulses:
                                # the pulses are close in time: merge
                                merged_hits[index - 1].charge += new_hit.charge
                            else:
                                # insert pulse
                                merged_hits.insert(index, new_hit)
                                len_merged_hits += 1

                            pulse_is_merged = True

                        # only increase index if we did not combine 2 pulses
                        if not combine_pulses:
                            index += 1

                # overwrite old pulse series
                pulse_series[key] = merged_hits

                # # sanity check
                # assert len(pulse_series[key]) == len_merged_hits

            # # sanity checks
            # t_previous = pulse_series[key][0].time
            # for p in pulse_series[key][1:]:
            #     assert p.time >= t_previous
            #     t_previous = p.time

        return pulse_series

    def _get_pulses(self, frame):
        """Get the I3RecoPulseSeriesMap from the frame.

        Parameters
        ----------
        frame : I3Frame
            The current I3Frame.
        """
        pulses = dataclasses.I3RecoPulseSeriesMap.from_frame(
            frame, self.pulse_key,
        )
        return dataclasses.I3RecoPulseSeriesMap(pulses)

    def DAQ(self, frame):
        if self.current_daq_frame is None:
            self.current_daq_frame = frame

    def Physics(self, frame):
        if 'oversampling' in frame:
            oversampling = frame['oversampling']

            # Find out if a new event started
            if (oversampling['event_num_in_run']
                    != self.current_event_counter):
                # new event started:
                # push aggregated frame if it hasn't been yet
                if (self.current_aggregation_frame is not None and
                        self.pushed_frame_already is False):
                    self.push_aggregated_frame()

                # reset values for new event
                self.current_time_shift = frame['TimeShift'].value
                self.current_aggregation_frame = frame
                self.current_event_counter = oversampling['event_num_in_run']
                self.merged_pulse_series = self._get_pulses(frame)
                self.oversampling_counter = 1
                self.pushed_frame_already = False

            else:
                # same event, keep aggregating pulses
                new_pulses = self._get_pulses(frame)
                self.merged_pulse_series = self.merge_pulse_series(
                                        self.merged_pulse_series,
                                        new_pulses,
                                        frame['TimeShift'].value)
                self.oversampling_counter += 1

            # Find out if event ended
            if (self.oversampling_factor
                    == 1 + oversampling['oversampling_num']):

                if self.current_aggregation_frame is not None:
                    self.push_aggregated_frame()

        else:
            self.PushFrame(frame)


class GetMCPulses(icetray.I3ConditionalModule):

    """Creates I3RecoPulseSeriesMap from I3MCPESeriesMap and optionally
    creates and inserts new Physics-frames.
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter('I3MCPEorPulseSeriesMap', 'MC PE or Pulse series to use.',
                          'I3MCPESeriesMapWithoutNoise')
        self.AddParameter('OutputKey',
                          'Output key to which the MC hits will be stored.',
                          'MCPulses')
        self.AddParameter('CreatePFrames', 'Create P frames from q frames?.',
                          True)
        self.AddParameter('WriteToQFrame', 'Add MC pulses to Q-Frame instead of P-frame.', False)
        self.AddParameter('PulseWidth', 'The width of each MCPE in ns. ', 0.2)

    def Configure(self):
        """Configure the module.
        """
        self._mcpe_series = self.GetParameter('I3MCPEorPulseSeriesMap')
        self._output_key = self.GetParameter('OutputKey')
        self._create_p_frames = self.GetParameter('CreatePFrames')
        self._write_to_q_frame = self.GetParameter('WriteToQFrame')
        self._pulse_width = self.GetParameter('PulseWidth')

        assert isinstance(self._create_p_frames, bool), \
            'Expected CreatePFrames to be a boolean, but got {!r}'.format(
                self._create_p_frames)

    def DAQ(self, frame):
        """Create P-frames and add MC pulses.

        Parameters
        ----------
        frame : I3Frame
            The current I3Frame.
        """
        if self._write_to_q_frame:
            # add MC reco pulses from I3MCPESeriesMap
            self._add_mc_pulses(
                frame, frame[self._mcpe_series], self._pulse_width
            )

            # Detector simulation creates trigger and shifts times relative
            # to this trigger. If detector simulation is skipped,
            # we must manually add the TimeShift key to the frame.
            if 'TimeShift' not in frame:
                frame['TimeShift'] = dataclasses.I3Double(0.)

        self.PushFrame(frame)

        if self._create_p_frames:
            p_frame = icetray.I3Frame(icetray.I3Frame.Physics)

            if not self._write_to_q_frame:
                # add MC reco pulses from I3MCPESeriesMap
                self._add_mc_pulses(
                    p_frame, frame[self._mcpe_series], self._pulse_width
                )

                # Detector simulation creates trigger and shifts times relative
                # to this trigger. If detector simulation is skipped,
                # we must manually add the TimeShift key to the frame.
                if 'TimeShift' not in frame:
                    p_frame['TimeShift'] = dataclasses.I3Double(0.)
            self.PushFrame(p_frame)

    def Physics(self, frame):
        """Add MC pulses to P-frame

        Parameters
        ----------
        frame : I3Frame
            The current I3Frame.
        """
        if not self._create_p_frames and not self._write_to_q_frame:

            # add MC reco pulses from I3MCPESeriesMap
            self._add_mc_pulses(
                frame, frame[self._mcpe_series], self._pulse_width
            )

        # push frame on to next modules
        self.PushFrame(frame)

    def _add_mc_pulses(self, frame, mcpe_series_map, pulse_width):
        '''Create MC reco pulses from I3MCPESeriesMap

        This is a dirty hack, so that other modules can be used without
        changing them. However, this will use up unecessary space, because
        I3RecoPulses have more data fields, which are not required by an
        MC hit (width, ATWD, ...) .

        Parameters
        ----------
        frame : I3Frame
            The I3Frame to which the MC Pulses will be added to.
        mcpe_series_map : I3MCPESeriesMap
            The I3MCPESeriesMap which will be converted.
        pulse_width : float
            The width of each MCPE in ns. The created
            I3RecoPulses will have this width.
        '''
        mc_pulse_map = dataclasses.I3RecoPulseSeriesMap()
        for omkey, mcpe_series in mcpe_series_map.items():

            mc_pulses = []
            for idx, mcpe in enumerate(mcpe_series):

                if isinstance(mcpe, simclasses.I3MCPE):
                    # create I3RecoPulse with corresponding time and 'charge'
                    # The charge is set to the number of photo electrons (npe)
                    mc_pulse = dataclasses.I3RecoPulse()
                    mc_pulse.time = mcpe.time
                    mc_pulse.charge = mcpe.npe
                    mc_pulse.width = pulse_width

                    # mis-use flags to store if the pulse is from a
                    # primary particle or noise
                    if mcpe.ID == dataclasses.I3ParticleID(0, 0):
                        mc_pulse.flags = 0
                    else:
                        mc_pulse.flags = 1

                elif isinstance(mcpe, simclasses.I3MCPulse):
                    # create I3RecoPulse with corresponding time and 'charge'
                    # The charge is set to the number of photo electrons (npe)
                    mc_pulse = dataclasses.I3RecoPulse()
                    mc_pulse.time = mcpe.time
                    mc_pulse.charge = mcpe.charge
                    mc_pulse.width = pulse_width

                    # check if we can find the origin particle ID
                    id_map = frame[self._mcpe_series + "ParticleIDMap"][omkey]
                    origin_id = None
                    for particle_id, idx_list in id_map:
                        if idx in idx_list:
                            origin_id = particle_id
                            break

                    # mis-use flags to store if the pulse is from a
                    # primary particle or noise
                    if origin_id is None:
                        mc_pulse.flags = 0
                    else:
                        mc_pulse.flags = 1

                else:
                    msg = 'Expected I3MCPulse or I3MCPE, but got {!r}'
                    raise TypeError(msg.format(type(mcpe)))

                # append pulse
                mc_pulses.append(mc_pulse)

            mc_pulse_map[omkey] = dataclasses.vector_I3RecoPulse(mc_pulses)

        # write to frame
        frame[self._output_key] = mc_pulse_map


class CompressPulses(icetray.I3ConditionalModule):
    """Compress charge and time from pulses.

    Note: compression of pulses will modify the pulse flags
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter('PulseKeys', 'The pulses to compress.',)
        self.AddParameter(
            'OutputKeys', 'Output keys for compressed pulses.',  None,
        )
        self.AddParameter(
            'RunOnQFrames',
            'If True, run on Q-frames, otherwise P-frames',
            True,
        )

    def Configure(self):
        """Configure the module.
        """
        self._pulse_keys = self.GetParameter('PulseKeys')
        self._output_keys = self.GetParameter('OutputKeys')
        self._run_on_q_frames = self.GetParameter('RunOnQFrames')

        if self._output_keys is None:
            self._output_keys = self._pulse_keys
        else:
            assert len(self._output_keys) == len(self._pulse_keys)

    def DAQ(self, frame):
        """Merge pulses on DAQ frames.

        Parameters
        ----------
        frame : I3Frame
            The current I3Frame.
        """
        if self._run_on_q_frames:
            self.compress_pulses(frame)
        self.PushFrame(frame)

    def Physics(self, frame):
        """Merge pulses on Physics frames.

        Parameters
        ----------
        frame : I3Frame
            The current I3Frame.
        """
        if not self._run_on_q_frames:
            self.compress_pulses(frame)
        self.PushFrame(frame)

    def discard_early_pulses(
        self, pulses, min_time=-512, max_removed=100, max_fraction=0.45,
    ):
        """Discard pulses that are too early.

        I3SuperDST expects pulses to have times later than -512ns. This
        function removes pulses that are too early.

        Parameters
        ----------
        pulses : I3RecoPulseSeriesMap
            The pulses to filter.
        min_time : float, optional
            The minimum time in ns a pulse must have to be kept.
        max_removed : int, optional
            The maximum number of pulses that are allowed to be removed.
            If more pulses are removed, an error is raised.
        max_fraction : float, optional
            The maximum fraction of pulses that are allowed to be removed.
            If more pulses are removed, an error is raised.

        Returns
        -------
        Same as input, but with pulses that are too early removed.
        """
        n_removed = 0
        n_total = 0
        new_pulses = type(pulses)()
        for omkey, pulse_series in pulses.items():
            new_pulses[omkey] = type(pulse_series)([
                pulse for pulse in pulse_series if pulse.time > min_time
            ])
            n_removed += len(pulse_series) - len(new_pulses[omkey])
            n_total += len(pulse_series)

        if n_removed > 0:
            log_error('Removed {} pulses with time < {}ns.'.format(
                n_removed, min_time))

        fraction_removed = n_removed / n_total
        if n_removed > max_removed or fraction_removed > max_fraction:
            print(pulses)
            log_fatal(
                'Removed {} pulses from {} [{:6.3f} %], but only {}'
                '[{:6.3f} %] pulses are allowed to be removed.'.format(
                    n_removed,
                    n_total,
                    fraction_removed*100,
                    max_removed,
                    max_fraction*100,
            ))

        return new_pulses

    def compress_pulses(self, frame):
        for pulse_key, output_key in zip(
            self._pulse_keys, self._output_keys,
        ):
            pulses = frame[pulse_key]
            pulses_trimmed = self.discard_early_pulses(pulses)
            if output_key in frame:
                del frame[output_key]
            frame[output_key] = dataclasses.I3SuperDST(pulses_trimmed)

class MergePulsesNearbyInTime(icetray.I3ConditionalModule):

    """Merge pulses nearby in time.

    Note: merging of pulses will disregard the pulse flags and widths!
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter('PulseKeys', 'The pulses to merge.',)
        self.AddParameter(
            'TimeThreshold',
            'Pulses within this time difference are merged togehter. '
            'Note that this will dicsard the pulse flags and widths!',
        )
        self.AddParameter(
            'OutputKeys',
            'Output keys for merged pulses. If None, input is overwritten.',
            None,
        )
        self.AddParameter(
            'RunOnQFrames',
            'If True, run on Q-frames, otherwise P-frames',
            True,
        )

    def Configure(self):
        """Configure the module.
        """
        self._pulse_keys = self.GetParameter('PulseKeys')
        self._time_threshold = self.GetParameter('TimeThreshold')
        self._output_keys = self.GetParameter('OutputKeys')
        self._run_on_q_frames = self.GetParameter('RunOnQFrames')

        if self._output_keys is None:
            self._output_keys = self._pulse_keys
        else:
            assert len(self._output_keys) == len(self._pulse_keys)

    def DAQ(self, frame):
        """Merge pulses on DAQ frames.

        Parameters
        ----------
        frame : I3Frame
            The current I3Frame.
        """
        if self._run_on_q_frames:
            # add MC reco pulses from I3MCPESeriesMap
            merge_pulses_in_time(
                frame=frame,
                pulse_keys=self._pulse_keys,
                time_threshold=self._time_threshold,
                output_keys=self._output_keys,
            )

        self.PushFrame(frame)

    def Physics(self, frame):
        """Merge pulses on Physics frames.

        Parameters
        ----------
        frame : I3Frame
            The current I3Frame.
        """
        if not self._run_on_q_frames:
            # add MC reco pulses from I3MCPESeriesMap
            merge_pulses_in_time(
                frame=frame,
                pulse_keys=self._pulse_keys,
                time_threshold=self._time_threshold,
                output_keys=self._output_keys,
            )

        self.PushFrame(frame)


def merge_pulses_in_time(frame, pulse_keys, time_threshold, output_keys=None):
    """Merge pulses nearby in time.

    Note: merging of pulses will disregard the pulse flags and widths!

    Parameters
    ----------
    frame : I3Frame
        The current I3Frame.
    pulse_keys : list of str
        The pulse keys to merge.
    time_threshold : float
        The time threshold in ns. Pulses that fall within this time window
        will be merged together.
    output_keys : list of str, optional
        The output keys to store the merged pulses. If None, the original
        pulse key will be overwritten.
    """
    from ic3_data.ext_boost import merge_pulses

    if output_keys is None:
        output_keys = pulse_keys
    else:
        assert len(output_keys) == len(pulse_keys)


    for pulse_key, output_key in zip(pulse_keys, output_keys):
        merged_pulses = merge_pulses(frame, pulse_key, time_threshold)

        if output_key in frame:
            del frame[output_key]
        frame[output_key] = merged_pulses
