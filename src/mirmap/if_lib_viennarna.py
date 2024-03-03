#
# Copyright © 2024 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

from dataclasses import dataclass

from ViennaRNA import RNA


@dataclass(frozen=True)
class RnacofoldResult:
    """RNAcofold results."""

    mfe_structure: str
    mfe: float


@dataclass(frozen=True)
class RnacofoldPartfuncResult:
    """RNAcofold with partition function results."""

    mfe_structure: str
    mfe: float
    efe_structure: str
    efe: float
    efe_binding: float


def get_model_details(min_loop_size=None, temperature=None):
    """New folding model with user parameters."""
    md = RNA.md()
    if min_loop_size is not None:
        md.min_loop_size = min_loop_size
    if temperature is not None:
        md.temperature = temperature
    return md


def RNAfold(seqs, constraint=None, partfunc=False, md=None):
    """Function RNAfold."""
    # Init
    fc = RNA.fold_compound("&".join(seqs), md)
    # Add constraint
    if constraint is not None:
        fc.constraints_add(constraint, RNA.CONSTRAINT_DB_DEFAULT)
    # MFE
    mfe_structure, mfe = fc.mfe()
    # Partition function
    if partfunc:
        efe_structure, fa, fb, fcab, efe = fc.pf_dimer()
        efe_binding = fcab - fa - fb
        return RnacofoldPartfuncResult(mfe_structure, mfe, efe_structure, efe, efe_binding)
    else:
        return RnacofoldResult(mfe_structure, mfe)
