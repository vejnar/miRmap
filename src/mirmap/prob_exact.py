#
# Copyright © 2024 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Probability based on exact distribution feature."""

from . import if_lib_spatt, prob


def calc_prob_exact(target, libspatt=None, markov_order=1, alphabet=None, transitions=None):
    """Compute the *P.over binomial* score.

    Args:
        target: Target
        libspatt: Link to the Spatt library
        markov_order: Markov Chain order
        alphabet: List of nucleotides to consider in the sequences
        transitions: Transition matrix of the Markov Chain model
    """
    # Parameters
    if libspatt is not None:
        if_spatt = if_lib_spatt.Spatt()
    if alphabet is None:
        alphabet = list(set(target.host_seq))
    if transitions is None:
        transitions = prob.get_transitions(target.host_seq, alphabet, markov_order)
    # Target seed binding sequence
    target_seed_seq = target.host_seq[target.seed.start : target.seed.end]
    return if_spatt.get_exact_prob(
        target_seed_seq,
        target.host_seq.count(target_seed_seq),
        len(target.host_seq),
        alphabet,
        transitions,
        markov_order,
        "o",
    )
