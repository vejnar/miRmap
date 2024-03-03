#
# Copyright © 2024 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Interface classes with the [PHAST](http://compgen.bscb.cornell.edu/phast) executable programs."""

import contextlib
import re
import subprocess

from .utils import closeable_temp_file


def phylofit(
    subst_model=None,
    aln_fname=None,
    aln=None,
    aln_format="FASTA",
    tree=None,
    use_em=None,
    path_exe="phyloFit",
):
    """Run phyloFit."""
    assert aln_fname is not None or aln is not None, "Input alignment is required"
    # Cmd
    cmd = [path_exe, "--precision", "HIGH", "--out-root", "-", "--msa-format", aln_format]
    # Options
    if subst_model is not None:
        cmd.append("--subst-mod")
        cmd.append(subst_model)
    if use_em is True:
        cmd.append("--EM")
    # Tree
    with (
        closeable_temp_file(suffix=".nh") if tree is not None else contextlib.nullcontext() as ftree,
        closeable_temp_file(suffix=".fa") if aln is not None else contextlib.nullcontext() as faln,
    ):
        if tree is not None:
            ftree.write(tree)
            ftree.close()
            cmd.append("--tree")
            cmd.append(ftree.name)
        if aln is not None:
            faln.write(aln)
            faln.close()
            cmd.append(faln.name)
        else:
            cmd.append(aln_fname)
        # Run
        p = subprocess.run(cmd, capture_output=True, text=True, check=True)
    # Parsing results
    decoded = re.match(
        (
            r"ALPHABET: (?P<alphabet>[^\n]+)\nORDER: (?P<order>\S+)\nSUBST_MOD: (?P<subst_mod>\S+)\n"
            r"TRAINING_LNL: (?P<training_lnl>\S+)\nBACKGROUND: (?P<background>[^\n]+)\nRATE_MAT:\n"
            "(?P<rate_mat>.+)\nTREE: (?P<tree>.+;)"
        ),
        p.stdout,
        re.DOTALL,
    )
    result = decoded.groupdict()
    result["mod_raw"] = p.stdout
    result["training_lnl"] = float(result["training_lnl"])
    return result


def phylop(
    method,
    mode,
    mod_fname,
    aln_fname=None,
    aln=None,
    aln_format="FASTA",
    gff_fname=None,
    branch=None,
    prune=None,
    path_exe="phyloP",
):
    """Run phyloP."""
    assert aln_fname is not None or aln is not None, "Input alignment is required"
    # Cmd
    cmd = [str(path_exe), "--method", method, "--mode", mode, "--msa-format", aln_format]
    # Options
    if gff_fname is not None:
        cmd.append("--features")
        cmd.append(gff_fname)
    if branch is not None:
        cmd.append("--branch")
        cmd.append(branch)
    if prune is False:
        cmd.append("--no-prune")
    # Model
    cmd.append(mod_fname)
    # MSA
    with closeable_temp_file(suffix=".fa") if aln is not None else contextlib.nullcontext() as faln:
        if aln is not None:
            faln.write(aln)
            faln.close()
            cmd.append(faln.name)
        else:
            cmd.append(aln_fname)
        # Run
        p = subprocess.run(cmd, capture_output=True, text=True, check=True)
    # Parse results
    decoded = re.search(r"p-value of conservation: (?P<prob>\S+)", p.stdout)
    return float(decoded.groupdict()["prob"])
