#
# Copyright © 2024 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Interval class."""

from dataclasses import dataclass


@dataclass(frozen=True, order=True)
class Interval:
    """Coordinates interval (0-based with inclusive start and exclusive end).

    Args:
        start: Inclusive 0-based start coordinate
        end: Exclusive 0-based end coordinate
    """

    start: int
    end: int

    def __post_init__(self):
        """Check coordinates are start then end."""
        assert self.start <= self.end, f"{self.start} must be <= {self.end}"

    def __len__(self):
        """Length of interval or number of nucleotides in the interval."""
        return self.end - self.start
