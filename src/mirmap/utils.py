#
# Copyright © 2024 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Utility functions."""

import contextlib
import gzip
import itertools
import os
import pathlib
import tempfile


_ttable = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq):
    """Return the reverse complement sequence.

    Args:
        seq: Sequence
    """
    return seq[::-1].translate(_ttable)


def flatten(l1d):
    """Flatten list by dimension."""
    return list(itertools.chain.from_iterable(l1d))


def grouper(n, iterable, fillvalue=None):
    """https://docs.python.org/3/library/itertools.html#itertools-recipes"""
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)


def clean_seq(seq, alphabet, replace_string):
    """Clean sequence by replacing string not in declared alphabet."""
    for nt in tuple(set(seq)):
        if nt not in alphabet:
            seq = seq.replace(nt, replace_string)
    return seq


def load_fasta(fasta, as_string=False, upper=False):
    """Parse FASTA."""
    if as_string:
        ff = fasta.split("\n")
    else:
        if pathlib.Path(fasta).suffix == ".gz":
            ff = gzip.open(fasta, "rb")
        else:
            ff = open(fasta)
    seqs = {}
    name_seq = ""
    for line in ff:
        line = line.strip()
        if line.startswith(">"):
            name_seq = line[1:].strip()
            seqs[name_seq] = ""
        else:
            if upper:
                seqs[name_seq] += line.upper()
            else:
                seqs[name_seq] += line
    if as_string is False:
        ff.close()
    return seqs


@contextlib.contextmanager
def closeable_temp_file(*args, **kwds):
    """Yield a closeable temporary file."""
    ftmp = tempfile.NamedTemporaryFile(delete=False, mode="wt", suffix=kwds["suffix"])
    try:
        yield ftmp
    finally:
        os.remove(ftmp.name)
