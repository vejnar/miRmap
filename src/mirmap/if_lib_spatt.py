#
# Copyright © 2024 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Ctypes interface classe with the [Spatt](http://www.mi.parisdescartes.fr/~nuel/spatt) C library."""

from ctypes import POINTER, c_bool, c_char, c_char_p, c_double, c_long, c_short, c_ulong, cast, cdll

from . import utils


class Spatt:
    """Interface class for the Spatt library."""

    def __init__(self, path_library="libspatt2.so"):
        """Create new interface to the Spatt library."""
        self._library = cdll.LoadLibrary(path_library)
        # Functions arguments and result types
        self._library.spatt_exact.argtypes = [
            c_char_p,
            c_char_p,
            POINTER(c_double),
            c_short,
            c_bool,
            c_long,
            c_ulong,
            c_char,
        ]
        self._library.spatt_exact.restype = c_double

    def transitions_l2c(self, transitions):
        """Take the transition matrix a list with element in C-order."""
        ctransitions = (c_double * len(transitions))()
        i = 0
        for v in transitions:
            ctransitions[i] = v
            i += 1
        return ctransitions

    def get_exact_prob(self, motif, nobs, length_seq, alphabet, transitions, markov_order, direction):
        """Compute exact probability."""
        ctransitions = self.transitions_l2c(utils.flatten(transitions))
        return self._library.spatt_exact(
            "".join(alphabet).encode("ascii"),
            motif.encode("ascii"),
            cast(ctransitions, POINTER(c_double)),
            markov_order,
            False,
            nobs,
            length_seq,
            direction.encode("ascii"),
        )
