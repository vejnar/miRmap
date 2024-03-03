#
# Copyright © 2024 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Probability based on binomial distribution feature."""

import math

from . import prob


def binomial_cdf(x, n, p):
    """Compute binomial CDF exxact solution using log to avoid float overflow errors.
    (https://stackoverflow.com/a/45869209)
    """
    cdf = 0
    b = 0
    for k in range(x + 1):
        if k > 0:
            b += math.log(n - k + 1) - math.log(k)
        log_pmf_k = b + k * math.log(p) + (n - k) * math.log(1 - p)
        cdf += math.exp(log_pmf_k)
    return cdf


def calc_prob_binomial(target, markov_order=1, alphabet=None, transitions=None):
    """Compute the *P.over binomial* score.

    Args:
        target: Target
        markov_order: Markov Chain order
        alphabet: List of nucleotides to consider in the sequences
        transitions: Transition matrix of the Markov Chain model
    """
    # Parameters
    if alphabet is None:
        alphabet = list(set(target.host_seq))
    if transitions is None:
        transitions = prob.get_transitions(target.host_seq, alphabet, markov_order)
    # Target seed binding sequence
    target_seed_seq = target.host_seq[target.seed.start : target.seed.end]
    return 1.0 - binomial_cdf(
        target.host_seq.count(target_seed_seq),
        len(target.host_seq) - len(target_seed_seq) + 1,
        prob.prob_motif(target_seed_seq, alphabet, markov_order, transitions),
    )
