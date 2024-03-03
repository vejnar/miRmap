#
# Copyright © 2024 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Evolutionary features."""

import copy
from dataclasses import dataclass

import dendropy

from . import if_exe_phast, utils


@dataclass(frozen=True)
class TargetAln:
    """Multi-species alignment of target."""

    ref_species: str
    seqs: dict
    species: set


def extract_tree_from_mod(mod=None, mod_fname=None):
    """Extract tree from model or model file."""
    if mod_fname is not None:
        with open(mod_fname, "rt") as f:
            mod = f.read()
    start = mod.find("TREE: ") + 6
    end = mod.find(";", start) + 1
    return mod[start:end]


def get_coord_vec(seq, alphabet):
    """Get coordinates of positions with letter from alphabet."""
    coord_vec = []
    for i, s in enumerate(seq):
        if s in alphabet:
            coord_vec.append(i)
    return coord_vec


def remove_gap_column(aln):
    """Remove only-gap position(s) in multi-sequence alignment."""
    clean_aln = copy.copy(aln)
    keep_cols = []
    alns = list(aln.values())
    for i in range(len(alns[0])):
        only_gap = True
        for seq in alns:
            if seq[i] != "-":
                only_gap = False
                continue
        if only_gap is False:
            keep_cols.append(i)
    for seq_name in aln.keys():
        clean_aln[seq_name] = "".join([aln[seq_name][i] for i in keep_cols])
    return clean_aln


def get_target_alns(target, aln_fname=None, aln=None, aln_alphabet=("A", "C", "G", "T", "N")):
    """Extract target alignments from alignment."""
    assert aln_fname is not None or aln is not None, "Input alignment is required"
    # Load alignment and remove gaps
    if aln_fname is not None:
        seqs = utils.load_fasta(aln_fname, as_string=False, upper=True)
    else:
        seqs = utils.load_fasta(aln, as_string=True, upper=True)
    # First sequence is reference
    ref_species = list(seqs.keys())[0]
    ref_seq_coords = get_coord_vec(seqs[ref_species], aln_alphabet)
    # Target seed binding sequence
    target_seed_seq = target.host_seq[target.seed.start : target.seed.end]
    # Extract alignment
    start_seed_in_aln = ref_seq_coords[target.seed.start]
    end_seed_in_aln = ref_seq_coords[target.seed.end]
    partial_seqs = {}
    with_motifs = set()
    for seq_name, seq in seqs.items():
        aln_seq = seq[start_seed_in_aln:end_seed_in_aln]
        if len(aln_seq) != aln_seq.count("-"):
            partial_seqs[seq_name] = aln_seq
            if target_seed_seq == utils.clean_seq(aln_seq, ["A", "C", "G", "T"], ""):
                with_motifs.add(seq_name)
    return TargetAln(ref_species, remove_gap_column(partial_seqs), with_motifs)


def calc_cons_bls(
    tree,
    fitting_tree=True,
    target_alns=None,
    target=None,
    aln_fname=None,
    aln=None,
    aln_alphabet=("A", "C", "G", "T", "N"),
    subst_model="REV",
    use_em=True,
    path_phyfit="phyloFit",
):
    """Compute the Branch Length Score (*BLS*).

    Args:
        tree: Tree in the Newick format
        fitting_tree: Fitting or not the tree on the alignment
        target_alns: Target alignments
        target: Target
        aln_fname: Alignment filename
        aln: Alignment
        aln_alphabet: List of nucleotides to consider in the aligned sequences (others get filtered)
        subst_model: PhyloFit substitution model (e.g. REV)
        use_em: Fitting or not the tree with Expectation-Maximization algorithm
        path_phyfit: Path to phyloFit executable
    """
    # Get alignments
    if target_alns is None:
        target_alns = get_target_alns(target, aln_fname, aln, aln_alphabet)
    # BLS
    if target_alns.ref_species in target_alns.species and len(target_alns.species) > 1:
        # Fitting tree if necessary
        if fitting_tree:
            if aln_fname is not None:
                fitted_tree = if_exe_phast.phylofit(
                    subst_model=subst_model, aln_fname=aln_fname, tree=tree, use_em=use_em, path_exe=path_phyfit
                )["tree"]
            elif aln is not None:
                fitted_tree = if_exe_phast.phylofit(
                    subst_model=subst_model, aln=aln, tree=tree, use_em=use_em, path_exe=path_phyfit
                )["tree"]
            else:
                raise ValueError("Missing alignment")
        else:
            fitted_tree = tree
        # Compute BLS
        dtree = dendropy.Tree.get_from_string(fitted_tree, schema="newick", preserve_underscores=True)
        dtree.retain_taxa_with_labels(target_alns.species)
        return sum([edge.length for edge in dtree.postorder_edge_iter()][:-1])
    else:
        return 0.0


def calc_selec_phylop(
    mod_fname,
    target_alns=None,
    target=None,
    aln_fname=None,
    aln=None,
    aln_alphabet=("A", "C", "G", "T", "N"),
    method="SPH",
    mode="CONACC",
    path_phylop="phyloP",
):
    """Compute the *PhyloP* score.

    Args:
        mod_fname: Model filename
        target_alns: Target alignments
        target: Target
        aln_fname: Alignment filename
        aln: Alignment
        aln_alphabet: List of nucleotides to consider in the aligned sequences (others get filtered)
        method: Test name performed by PhyloP (e.g. SPH)
        mode: Testing for conservation (CON), acceleration (ACC) or both (CONACC)
        path_phylop: Path to phyloP executable
    """
    # Get alignments
    if target_alns is None:
        target_alns = get_target_alns(target, aln_fname, aln, aln_alphabet)
    # Run phyloP
    if target_alns.ref_species in target_alns.species and len(target_alns.species) > 1:
        # Extract alignment
        aln = "\n".join([f">{seq_name}\n{seq}" for seq_name, seq in target_alns.seqs.items()])
        # Compute p-value
        return if_exe_phast.phylop(
            method=method, mode=mode, mod_fname=mod_fname, aln=aln, aln_format="FASTA", path_exe=path_phylop
        )
    else:
        return 1.0
