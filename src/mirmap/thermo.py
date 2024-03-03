#
# Copyright © 2024 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Thermodynamics features."""

from .if_lib_viennarna import RNAfold, get_model_details


def _get_default_md():
    return get_model_details(min_loop_size=2)


def calc_dg_duplex(target, md=None):
    """Compute the *ΔG duplex*, *ΔG binding* scores.

    Args:
        target: Target
        md: Folding model
    """
    # Default model
    if md is None:
        md = _get_default_md()
    # Target site sequence
    target_site_seq = target.host_seq[target.mirna.start : target.mirna.end]
    # Constraint sequence
    constraint_up = ["."] * target.mirna_length
    constraint_up[(target.seed_length + target.seed_start_on_mirna) * -1 : target.seed_start_on_mirna * -1] = [
        "("
    ] * target.seed_length
    constraint_dn = ["."] * target.mirna_length
    constraint_dn[(target.seed_length + target.seed_start_on_mirna) * -1 : target.seed_start_on_mirna * -1] = [
        ")"
    ] * target.seed_length
    # Fold
    result = RNAfold(
        [target_site_seq, target.mirna_seq],
        constraint="".join(constraint_up + constraint_dn[::-1]),
        partfunc=True,
        md=md,
    )
    return {
        "dg_duplex": result.mfe,
        "dg_binding": result.efe_binding,
    }


def calc_dg_duplex_seed(target, md=None):
    """Compute the *ΔG seed duplex* and *ΔG seed binding* scores.

    Args:
        target: Target
        md: Folding model
    """
    # Default model
    if md is None:
        md = _get_default_md()
    # Target seed binding sequence
    target_seed_seq = target.host_seq[target.seed.start : target.seed.end]
    # Fold
    result = RNAfold(
        [target_seed_seq, target.seed_seq],
        partfunc=True,
        md=md,
    )
    return {
        "dg_duplex_seed": result.mfe,
        "dg_binding_seed": result.efe_binding,
    }


def calc_dg_open(target, upstream_rest=20, downstream_rest=20, dg_binding_area=70, md=None):
    """Compute the *ΔG open* score.

    Args:
        target: Target
        upstream_rest: Upstream unfolding length
        downstream_rest: Downstream unfolding length
        dg_binding_area: Supplementary sequence length to fold (applied upstream and downstream)
        md: Folding model
    """
    # Default model
    if md is None:
        md = _get_default_md()
    start_dgopen_tmp = target.mirna.start - upstream_rest - dg_binding_area
    end_dgopen_tmp = target.mirna.end + downstream_rest + dg_binding_area
    if start_dgopen_tmp < 0:
        start_dgopen = 0
        len_polya_upstream = abs(start_dgopen_tmp)
    else:
        start_dgopen = start_dgopen_tmp
        len_polya_upstream = 0
    if end_dgopen_tmp > len(target.host_seq):
        end_dgopen = len(target.host_seq)
        len_polya_downstream = end_dgopen_tmp - len(target.host_seq)
    else:
        end_dgopen = end_dgopen_tmp
        len_polya_downstream = 0
    seq_for_dg_open = len_polya_upstream * "A" + target.host_seq[start_dgopen:end_dgopen] + len_polya_downstream * "A"
    # Constraint sequence
    constraint = (
        "." * dg_binding_area + "x" * (upstream_rest + target.mirna_length + downstream_rest) + "." * dg_binding_area
    )
    # Folding
    # dg0
    result_dg0 = RNAfold(
        [seq_for_dg_open],
        partfunc=True,
        md=md,
    )
    # dg1
    result_dg1 = RNAfold(
        [seq_for_dg_open],
        constraint=constraint,
        partfunc=True,
        md=md,
    )
    # dg_open
    return {"dg_open": result_dg1.efe - result_dg0.efe}


def calc_dg_total(target, dg_duplex, dg_open, md=None):
    """Compute the *ΔG total* score combining *ΔG duplex* and *ΔG open* scores.

    Args:
        target: Target
        dg_duplex: ΔG duplex score
        dg_open: ΔG open score
        md: Folding model
    """
    # Default model
    if md is None:
        md = _get_default_md()
    # Calculate ΔG duplex and ΔG open
    if dg_duplex is None:
        dg_duplex = calc_dg_duplex(target, md=md)
    if dg_open is None:
        dg_open = calc_dg_open(target, md=md)
    return {"dg_total": dg_duplex + dg_open}
