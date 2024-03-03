#
# Copyright © 2024 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""TargetScan features."""

import operator

from .targetscan_model import ts_types


def get_targetscan_ts_type(seed_length, nt1):
    """Return the target site type according to TargetScan terminology.

    Args:
        seed_length: Seed length
        nt1: Last nucleotide of the target on the host sequence
    """
    if seed_length >= 7:
        if nt1 == "A":
            return "8mer"
        else:
            return "7mer-m8"
    elif seed_length == 6:
        if nt1 == "A":
            return "7mer-A1"
        else:
            return "6mer"
    return None


def _binarize(c):
    if c == "T" or c == "A":
        return 1.0
    else:
        return 0.0


def calc_tgs_au(target, ts_type=None, ts_types=ts_types, ca_window_length=30, with_correction=False):
    """Calculate the *AU content* score.

    Args:
        target: Target
        ts_type: Target site type (6mer, 7mer-A1, 7mer-m8 or 8mer)
        ts_types: TargetScan model parameters
        ca_window_length: Sequence window length used to compute the TargetScan score
        with_correction: Apply the linear regression correction or not
    """
    if ts_type is None:
        ts_type = get_targetscan_ts_type(target.seed_length, target.host_seq[target.seed.end])
    if ts_type:
        tts = ts_types[ts_type]
        seq_up = target.host_seq[
            max(0, target.seed.end + 1 + tts.up_shift - ca_window_length) : target.seed.end + 1 + tts.up_shift
        ]
        seq_down = target.host_seq[
            target.seed.end + tts.down_shift : min(
                len(target.host_seq), target.seed.end + ca_window_length + tts.down_shift
            )
        ]
        wup = tts.ca_weights_up[len(tts.ca_weights_up) - len(seq_up) :]
        wdn = tts.ca_weights_down[: len(seq_down)]
        content = sum(map(operator.truediv, map(_binarize, seq_up), wup)) + sum(
            map(operator.truediv, map(_binarize, seq_down), wdn)
        )
        content = content / (
            sum(map(operator.truediv, [1.0] * len(wup), wup)) + sum(map(operator.truediv, [1.0] * len(wdn), wdn))
        )
        if with_correction:
            return content * tts.ca_fc_slope + tts.ca_fc_intercept - tts.fc_mean
        else:
            return content
    else:
        return None


def calc_tgs_position(target, ts_type=None, ts_types=ts_types, with_correction=False):
    """Calculate the *UTR position* score.

    Args:
        target: Target
        ts_type: Target site type (6mer, 7mer-A1, 7mer-m8 or 8mer)
        ts_types: TargetScan model parameters
        with_correction: Apply the linear regression correction or not
    """
    if ts_type is None:
        ts_type = get_targetscan_ts_type(target.seed_length, target.host_seq[target.seed.end])
    if ts_type:
        tts = ts_types[ts_type]
        closest_term = min(
            target.seed.end + 1 + tts.up_shift, len(target.host_seq) - (target.seed.end + 1) + tts.down_shift
        )
        if closest_term > 1500:
            closest_term = 1500
        if with_correction:
            return float(closest_term) * tts.po_fc_slope + tts.po_fc_intercept - tts.fc_mean
        else:
            return float(closest_term)
    else:
        return None


def _align(utr_3p_seq, mir_3p_seq, mir_offset, utr_offset, overhang):
    score = 0
    tempscore = 0
    prevmatch = 0
    i = 0
    offset = max(mir_offset, utr_offset)
    while (i < len(mir_3p_seq) - mir_offset) and (i < len(utr_3p_seq) - utr_offset):
        if (
            (utr_3p_seq[i + utr_offset] == "A" and mir_3p_seq[i + mir_offset] == "T")
            or (utr_3p_seq[i + utr_offset] == "T" and mir_3p_seq[i + mir_offset] == "A")
            or (utr_3p_seq[i + utr_offset] == "G" and mir_3p_seq[i + mir_offset] == "C")
            or (utr_3p_seq[i + utr_offset] == "C" and mir_3p_seq[i + mir_offset] == "G")
        ):
            if (i + mir_offset - overhang >= 4) and (i + mir_offset - overhang <= 7):
                if prevmatch == 0:
                    tempscore = 0
                tempscore += 1
            else:
                if prevmatch == 0:
                    tempscore = 0
                tempscore += 0.5
            prevmatch += 1
        elif prevmatch >= 2:
            if tempscore > score:
                score = tempscore
            tempscore = 0
            prevmatch = 0
        else:
            tempscore = 0
            prevmatch = 0
        i += 1
    if prevmatch >= 2:
        if tempscore > score:
            score = tempscore
        tempscore = 0
        prevmatch = 0
    score = score - max(0, ((offset - 2) / 2.0))
    return score


def calc_tgs_pairing3p(target, ts_type=None, ts_types=ts_types, with_correction=False):
    """Calculate the *3' pairing* score.

    Args:
        target: Target
        ts_type: Target site type (6mer, 7mer-A1, 7mer-m8 or 8mer)
        ts_types: TargetScan model parameters
        with_correction: Apply the linear regression correction or not
    """
    if ts_type is None:
        ts_type = get_targetscan_ts_type(target.seed_length, target.host_seq[target.seed.end])
    if ts_type:
        tts = ts_types[ts_type]
        utr_3p_seq = target.host_seq[
            max(0, (target.seed.end + 1) - tts.pa_mirna_seed_start - 15) : (target.seed.end + 1)
            - tts.pa_mirna_seed_start
        ][::-1]
        mir_3p_seq = target.mirna_seq[tts.pa_mirna_seed_start :]
        maxscore = max(len(utr_3p_seq), len(mir_3p_seq))
        scores_mir = []
        scores_utr = []
        for offset in range(maxscore):
            scores_mir.append(_align(utr_3p_seq, mir_3p_seq, offset, 0, tts.pa_mirna_seed_overhang))
            scores_utr.append(_align(utr_3p_seq, mir_3p_seq, 0, offset, tts.pa_mirna_seed_overhang))
        if with_correction:
            return float(max(scores_mir + scores_utr)) * tts.pa_fc_slope + tts.pa_fc_intercept - tts.fc_mean
        else:
            return float(max(scores_mir + scores_utr))
    else:
        return None


def calc_tgs_score(
    target,
    ts_types=ts_types,
    with_correction=False,
    score_tgs_au=None,
    score_tgs_position=None,
    score_tgs_pairing3p=None,
):
    """Compute the *TargetScan* score combining *AU content*, *UTR position* and *3' pairing* scores.

    Args:
        target: Target
        ts_types: TargetScan model parameters
        with_correction: Apply the linear regression correction or not
        score_tgs_au: AU content score
        score_tgs_position: UTR position score
        score_tgs_pairing3p: 3' pairing score
    """
    ts_type = get_targetscan_ts_type(target.seed_length, target.host_seq[target.seed.end])
    if ts_type:
        tts = ts_types[ts_type]
        if score_tgs_au is None:
            score_tgs_au = calc_tgs_au(target, ts_type=ts_type, ts_types=ts_types, with_correction=with_correction)
        if score_tgs_position is None:
            score_tgs_position = calc_tgs_position(
                target, ts_type=ts_type, ts_types=ts_types, with_correction=with_correction
            )
        if score_tgs_pairing3p is None:
            score_tgs_pairing3p = calc_tgs_pairing3p(
                target, ts_type=ts_type, ts_types=ts_types, with_correction=with_correction
            )
        try:
            score = score_tgs_au + score_tgs_position + score_tgs_pairing3p
            if with_correction:
                return score + tts.fc_mean
            else:
                return score
        except TypeError:
            return None
