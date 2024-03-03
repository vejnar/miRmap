#
# Copyright © 2024 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

from . import evolution, model, prob_binomial, prob_exact, targetscan, thermo


score_labels = {
    "dg_duplex": "ΔG duplex (kcal/mol)",
    "dg_binding": "ΔG binding (kcal/mol)",
    "dg_open": "ΔG open (kcal/mol)",
    "dg_total": "ΔG total (kcal/mol)",
    "tgs_au": "AU content",
    "tgs_position": "UTR position",
    "tgs_pairing3p": "3' pairing",
    "tgs_score": "TargetScan score",
    "prob_exact": "Probability (Exact)",
    "prob_binomial": "Probability (Binomial)",
    "cons_bls": "Conservation (BLS)",
    "selec_phylop": "Conservation (PhyloP)",
    "mirmap_score": "miRmap score",
}


score_agg_funcs = {
    "tgs_au": max,
    "tgs_position": min,
    "tgs_pairing3p": max,
    "tgs_score": sum,
    "dg_duplex": min,
    "dg_binding": min,
    "dg_duplex_seed": min,
    "dg_binding_seed": min,
    "dg_open": min,
    "dg_total": min,
    "prob_exact": min,
    "prob_binomial": min,
    "cons_bls": max,
    "selec_phylop": min,
    "mirmap_score": sum,
}


def calc_scores(
    target,
    rna_md=None,
    path_aln=None,
    path_mod=None,
    tree=None,
    if_spatt=None,
    path_phylofit=None,
    path_phylop=None,
):
    """Compute all scores including *miRmap* score.

    Args:
        target: Target
        rna_md: Folding model
        path_aln: Path to multiple sequence alignment
        path_mod: Path to evolutionary model
        tree: Newick species tree
        if_spatt: Interface object to the Spatt library
        path_phylofit: Path to the phyloFit executable
        path_phylop: Path to the phyloP executable
    """
    scores = {}

    # TargetScan features
    scores["tgs_au"] = targetscan.calc_tgs_au(target)
    scores["tgs_position"] = targetscan.calc_tgs_position(target)
    scores["tgs_pairing3p"] = targetscan.calc_tgs_pairing3p(target)
    scores["tgs_score"] = targetscan.calc_tgs_score(
        target,
        score_tgs_au=scores["tgs_au"],
        score_tgs_position=scores["tgs_position"],
        score_tgs_pairing3p=scores["tgs_pairing3p"],
    )

    # Thermodynamics features
    scores |= thermo.calc_dg_duplex(target, md=rna_md)
    scores |= thermo.calc_dg_duplex_seed(target, md=rna_md)
    scores |= thermo.calc_dg_open(target, md=rna_md)
    scores |= thermo.calc_dg_total(
        target,
        dg_duplex=scores["dg_duplex"],
        dg_open=scores["dg_open"],
        md=rna_md,
    )

    # Probabilistic features
    scores["prob_exact"] = prob_exact.calc_prob_exact(target, if_spatt)
    scores["prob_binomial"] = prob_binomial.calc_prob_binomial(target)

    # Evolutionary features
    scores["cons_bls"] = 0.0
    scores["selec_phylop"] = 1.0
    if path_aln is not None:
        target_alns = evolution.get_target_alns(target, aln_fname=path_aln)

        # BLS
        if tree is not None:
            # Fitting the species tree
            scores["cons_bls"] = evolution.calc_cons_bls(
                tree=tree,
                fitting_tree=True,
                aln_fname=path_aln,
                target_alns=target_alns,
                path_phylofit=path_phylofit,
            )
        elif path_mod is not None:
            # Using fitted tree
            scores["cons_bls"] = evolution.calc_cons_bls(
                tree=evolution.extract_tree_from_mod(mod_fname=path_mod),
                fitting_tree=False,
                target_alns=target_alns,
            )

        # PhyloP
        if path_mod is not None:
            scores["selec_phylop"] = evolution.calc_selec_phylop(
                mod_fname=path_mod,
                target_alns=target_alns,
                path_phylop=path_phylop,
            )

    # miRmap score
    scores["mirmap_score"] = model.calc_mirmap(target, model.full_mirmap_models, scores)

    return scores


def agg_scores(targets_scores):
    """Aggregate score from multiple targets.

    Args:
        targets_scores: Scores for all targets on one transcript
    """
    agg_scores = {}
    for name, fn in score_agg_funcs.items():
        agg_scores[name] = fn([target[name] for target in targets_scores])
    return agg_scores


def report_scores(scores):
    """Pretty target scores report.

    Args:
        scores: Target scores
    """
    lines = []
    for name, label in score_labels.items():
        lines.append(f" {label:<25}{scores[name]:.4}")
    return "\n".join(lines)
