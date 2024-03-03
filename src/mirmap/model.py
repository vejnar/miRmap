#
# Copyright © 2024 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""miRmap models and scores."""

from dataclasses import dataclass, fields


@dataclass(frozen=True, order=True)
class MirmapModel:
    """miRmap model base class."""

    def apply_on_target(self, scores):
        """Compute miRmap score using individual feature scores.

        Args:
            scores: Features scores
        """
        score = self.intercept
        for field in fields(self):
            try:
                score += scores[field.name] * getattr(self, field.name)
            except KeyError:
                if field.name != "intercept":
                    raise
        return score


@dataclass(frozen=True, order=True)
class FullMirmapModel(MirmapModel):
    """Full miRmap model including all features."""

    tgs_au: float
    tgs_position: float
    tgs_pairing3p: float
    dg_duplex: float
    dg_binding: float
    dg_duplex_seed: float
    dg_binding_seed: float
    dg_open: float
    prob_exact: float
    prob_binomial: float
    cons_bls: float
    selec_phylop: float
    intercept: float


@dataclass(frozen=True, order=True)
class PythonOnlyMirmapModel(MirmapModel):
    """miRmap model limited to Python-implemented features."""

    tgs_au: float
    tgs_position: float
    tgs_pairing3p: float
    prob_binomial: float
    cons_bls: float
    intercept: float


full_mirmap_models = {
    6: FullMirmapModel(
        tgs_au=-0.275016235769136,
        tgs_position=5.43367028065211e-06,
        tgs_pairing3p=-0.00233278119760994,
        dg_duplex=0.00772658898496047,
        dg_binding=-0.00303683833660696,
        dg_duplex_seed=0.0496909801533612,
        dg_binding_seed=-0.048931930580652,
        dg_open=0.000674676164622922,
        prob_exact=0.16111635592018,
        prob_binomial=-0.0388333740708671,
        cons_bls=-0.00426314077593848,
        selec_phylop=-0.0112455248228072,
        intercept=0.148300586692704,
    ),
    7: FullMirmapModel(
        tgs_au=-0.402470212080983,
        tgs_position=6.89249707831041e-05,
        tgs_pairing3p=-0.0129891251446967,
        dg_duplex=0.0141332997802509,
        dg_binding=-0.0132159175462755,
        dg_duplex_seed=-0.0814445085121904,
        dg_binding_seed=0.115558118311931,
        dg_open=0.00331507347139685,
        prob_exact=0.792962156550929,
        prob_binomial=-0.22119499646323,
        cons_bls=-0.0355840335642203,
        selec_phylop=-0.0127531995991629,
        intercept=0.349448109979275,
    ),
}


python_only_mirmap_models = {
    6: PythonOnlyMirmapModel(
        tgs_au=-0.275594504153219,
        tgs_position=9.44582844229299e-06,
        tgs_pairing3p=-0.0111209267382849,
        prob_binomial=0.0701619992923641,
        cons_bls=-0.00646548621345819,
        intercept=0.121104869645859,
    ),
    7: PythonOnlyMirmapModel(
        tgs_au=-0.443606032336791,
        tgs_position=6.34603935320321e-05,
        tgs_pairing3p=-0.0207672870210752,
        prob_binomial=0.378665477250754,
        cons_bls=-0.0552713344740971,
        intercept=0.150015113841088,
    ),
}


def calc_mirmap(target, models=None, scores=None):
    """Compute the *miRmap* score.

    Args:
        target: Target
        models: miRmap model
        scores: Target scores
    """
    return models[target.seed_length].apply_on_target(scores)
