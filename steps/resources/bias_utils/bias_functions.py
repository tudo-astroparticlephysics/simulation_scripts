import numpy as np

from icecube.icetray.i3logging import log_error, log_warn

from ic3_labels.labels.utils import geometry
from ic3_labels.labels.utils import detector
from ic3_labels.labels.utils import muon as mu_utils


class BaseBiasFunction:

    """Base bias function class
    """

    def __init__(self, **settings):
        """Summary

        Parameters
        ----------
        settings : dict
            The settings for this class.
        """
        pass

    def __call__(self, bias_data):
        """Apply Bias Function

        Parameters
        ----------
        bias_data : dict
            Dictionary of bias input data.
            Contents may include:
            {
                'frame': the current I3Frame,
            }

        Returns
        -------
        float
            Keep probability: probability with which this event should be kept.
        dict
            A dictionary with additional bias information.
        """
        return 1., {}

    def sigmoid(self, x, s=1., b=0.):
        """Compute Sigmoid Function

        Parameters
        ----------
        x : array_like
            The input data.
        s : float, optional
            The scale parameter of the sigmoid.
        b : float, optional
            The offset parameter of the sigmoid.

        Returns
        -------
        array_like
            Sigmoid results
        """
        xt = (x - b) / s
        return 1./(1 + np.exp(-xt))


class UpgoingMuonStochasticity(BaseBiasFunction):

    """Biases simulation towards upgoing muons with large relative
    energy losses.
    """

    def __init__(
            self,
            cos_zenith_sigmoid_scale=0.03,
            cos_zenith_sigmoid_bias=-0.2,
            track_length_sigmoid_scale=15,
            track_length_sigmoid_bias=120,
            muon_loss_sigmoid_scale=.03,
            muon_loss_sigmoid_bias=0.25,
            mctree_name='I3MCTree',
            ):
        """Summary

        Parameters
        ----------
        mctree_name : str, optional
            The name of the I3MCTree. The I3MCTree may be before CLSIM, but
            PROPOSAL needs to have run.
        """
        self.mctree_name = mctree_name
        self.cos_zenith_sigmoid_scale = cos_zenith_sigmoid_scale
        self.cos_zenith_sigmoid_bias = cos_zenith_sigmoid_bias
        self.track_length_sigmoid_scale = track_length_sigmoid_scale
        self.track_length_sigmoid_bias = track_length_sigmoid_bias
        self.muon_loss_sigmoid_scale = muon_loss_sigmoid_scale
        self.muon_loss_sigmoid_bias = muon_loss_sigmoid_bias

    def __call__(self, bias_data):
        """Apply Bias Function

        Parameters
        ----------
        bias_data : dict
            Dictionary of bias input data.
            Contents may include:
            {
                'frame': the current I3Frame,
            }

        Returns
        -------
        float
            Keep probability: probability with which this event should be kept.
        """

        frame = bias_data['frame']

        # get primary
        mc_tree = frame[self.mctree_name]
        primaries = mc_tree.get_primaries()
        assert len(primaries) == 1, 'Expected only 1 Primary!'

        # get muon
        muon = mu_utils.get_muon(
            frame, primaries[0], detector.icecube_hull,
            mctree_name=self.mctree_name,
        )

        if muon is None:

            # if muon did not hit the convex hull, or if no muon exists,
            # it will be None. In this case we set default values
            found_muon = False
            cos_zen = np.cos(primaries[0].dir.zenith)
            track_length = 0.
            max_rel_loss = 0.

        else:
            found_muon = True
            cos_zen = np.cos(muon.dir.zenith)
            track_length = mu_utils.get_muon_track_length_inside(
                muon, detector.icecube_hull)

            # get muon energy losses
            losses = [
                loss for loss in mc_tree.get_daughters(muon) if
                not mu_utils.is_muon(loss) and
                geometry.is_in_detector_bounds(loss.pos, extend_boundary=60)
            ]

            # compute relative energy losses
            rel_losses = []
            loss_energies = []
            for loss in losses:

                # get energy of muon prior to energy loss
                distance = (muon.pos - loss.pos).magnitude
                energy = mu_utils.get_muon_energy_at_distance(
                    frame, muon, np.clip(distance - 1, 0., float('inf')))

                # If the loss is at the muon decay point, the returned energy
                # might be NaN, assert this and set default value of 1 GeV
                if not np.isfinite(energy):
                    assert np.abs(distance - muon.length) < 1, (energy, muon)
                    energy = 1

                rel_loss = loss.energy / energy
                if rel_loss > 1. or rel_loss < 0.:
                    msg = 'Found out of bounds rel_loss: {:3.3f}. '.format(
                        rel_loss)
                    msg += 'Clipping value to [0, 1]'
                    log_warn(msg)
                    rel_loss = np.clip(rel_loss, 0., 1.)

                loss_energies.append(loss.energy)
                rel_losses.append(rel_loss)
            if rel_losses:
                max_rel_loss = rel_losses[np.argmax(loss_energies)]
            else:
                max_rel_loss = 0.

        # bias based on zenith
        if self.cos_zenith_sigmoid_scale is None:
            zenith_keep_prob = 1.0
        else:
            zenith_keep_prob = self.sigmoid(
                -cos_zen,
                s=self.cos_zenith_sigmoid_scale,
                b=self.cos_zenith_sigmoid_bias,
            )

        # bias based on in detector track length
        if self.track_length_sigmoid_scale is None:
            track_length_prob = 1.0
        else:
            track_length_prob = self.sigmoid(
                track_length,
                s=self.track_length_sigmoid_scale,
                b=self.track_length_sigmoid_bias,
            )

        # bias based on maximum relative energy loss
        if self.muon_loss_sigmoid_scale is None:
            max_rel_loss_prob = 1.
        else:
            max_rel_loss_prob = self.sigmoid(
                max_rel_loss,
                s=self.muon_loss_sigmoid_scale,
                b=self.muon_loss_sigmoid_bias,
            )

        bias_info = {
            'found_muon': found_muon,
            'cos_zenith': cos_zen,
            'track_length_in_detector': track_length,
            'max_relative_energy_loss': max_rel_loss,
        }

        keep_prob = zenith_keep_prob * track_length_prob * max_rel_loss_prob
        return keep_prob, bias_info


class DowngoingFirstPrimary(BaseBiasFunction):

    """Biases simulation towards downgoing events

    This will be based on the first primary in the provide I3MCTree.
    """

    def __init__(
            self,
            cos_zenith_sigmoid_scale=0.03,
            cos_zenith_sigmoid_bias=0.,
            mctree_name='I3MCTree',
            ):
        """Summary

        Parameters
        ----------
        mctree_name : str, optional
            The name of the I3MCTree. The I3MCTree may be before CLSIM, but
            PROPOSAL needs to have run.
        """
        self.mctree_name = mctree_name
        self.cos_zenith_sigmoid_scale = cos_zenith_sigmoid_scale
        self.cos_zenith_sigmoid_bias = cos_zenith_sigmoid_bias

    def __call__(self, bias_data):
        """Apply Bias Function

        Parameters
        ----------
        bias_data : dict
            Dictionary of bias input data.
            Contents may include:
            {
                'frame': the current I3Frame,
            }

        Returns
        -------
        float
            Keep probability: probability with which this event should be kept.
        """

        frame = bias_data['frame']

        # get primary
        mc_tree = frame[self.mctree_name]
        primaries = mc_tree.get_primaries()
        cos_zen = np.cos(primaries[0].dir.zenith)

        # bias based on zenith
        if self.cos_zenith_sigmoid_scale is None:
            zenith_keep_prob = 1.0
        else:
            zenith_keep_prob = self.sigmoid(
                cos_zen,
                s=self.cos_zenith_sigmoid_scale,
                b=self.cos_zenith_sigmoid_bias,
            )

        bias_info = {
            'cos_zenith': cos_zen,
        }

        return zenith_keep_prob, bias_info
