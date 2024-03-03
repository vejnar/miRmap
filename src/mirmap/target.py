#
# Copyright © 2024 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

import re

from . import utils
from .interval import Interval


class Target:
    """miRNA target class."""

    def __init__(self, host_seq, mirna_interval, mirna_seq, seed_interval, seed_seq, seed_start_on_mirna):
        """Create new Target object."""
        self.host_seq = host_seq
        self.mirna = mirna_interval
        self.mirna_length = len(mirna_interval)
        self.mirna_seq = mirna_seq
        self.seed = seed_interval
        self.seed_length = len(seed_interval)
        self.seed_seq = seed_seq
        self.seed_start_on_mirna = seed_start_on_mirna

    def report(self, extension=10):
        """Report miRNA and host sequence pairing."""
        window_start = self.mirna.start - extension
        window_end = self.mirna.end + extension
        lines = [
            " " * extension
            + str(self.mirna.start)
            + " " * (self.mirna_length - len(str(self.mirna.start)) - 1)
            + str(self.mirna.end - 1),
            " " * extension + "|" + " " * (self.mirna_length - 2) + "|",
            self.host_seq[window_start:window_end],
            " " * (self.seed.start - window_start) + "|" * self.seed_length,
            " " * extension + self.mirna_seq[::-1],
        ]
        return "\n".join(lines)


def find_targets_with_seed(host_seq, mirna_seq, seed_start_on_mirna=1, seed_lengths=None):
    """Find target sites in host sequence searching for miRNA seeds.

    Args:
        host_seq: Host sequence
        mirna_seq: miRNA sequence
        seed_start_on_mirna: Start position of the seed in the miRNA (from the 5')
        seed_lengths: List of seed length(s)
    """
    # Parameter
    if seed_lengths is None:
        seed_lengths = [6, 7]

    host_seq_rc = utils.reverse_complement(host_seq)
    targets = []
    target_ends = set()
    for seed_length in sorted(seed_lengths, reverse=True):
        seed_seq = mirna_seq[seed_start_on_mirna : seed_start_on_mirna + seed_length]
        for m in re.finditer(rf"(?=({seed_seq}))", host_seq_rc):
            seed_interval = Interval(len(host_seq) - m.end(1), len(host_seq) - m.start(1))
            # Check miRNA is fully within the host sequence
            if (
                seed_interval.start < (len(mirna_seq) - len(seed_seq) - seed_start_on_mirna)
                or seed_interval.end > len(host_seq) - seed_start_on_mirna
            ):
                continue
            if seed_interval.end not in target_ends:
                mirna_interval = Interval(
                    seed_interval.end - len(mirna_seq) + seed_start_on_mirna, seed_interval.end + seed_start_on_mirna
                )
                # Check interval, length and sequence fit together
                assert len(seed_interval) == seed_length, f"{len(seed_interval)} must be == {seed_length}"
                assert len(mirna_interval) == len(mirna_seq), f"{len(mirna_interval)} must be == {len(mirna_seq)}"
                targets.append(
                    Target(host_seq, mirna_interval, mirna_seq, seed_interval, seed_seq, seed_start_on_mirna)
                )
                target_ends.add(seed_interval.end)

    # Sort targets by their seed end position
    targets.sort(key=lambda x: x.seed.end, reverse=True)

    return targets
