import pytest

import mirmap.target
import mirmap.utils
from mirmap.interval import Interval
from mirmap.target import Target


@pytest.mark.unit()
@pytest.mark.parametrize(
    "fname_mirna,fname_transcript,expected_targets",
    [
        (
            "hsa-miR-30a-3p.fa",
            "NM_024573.fa",
            [
                Target(
                    host_seq="",
                    mirna_interval=Interval(start=909, end=931),
                    mirna_seq="CTTTCAGTCGGATGTTTGCAGC",
                    seed_interval=Interval(start=923, end=930),
                    seed_seq="TTTCAGT",
                    seed_start_on_mirna=1,
                )
            ],
        ),
        (
            "hsa-miR-3124-5p.fa",
            "NM_024573.fa",
            [],
        ),
        (
            "hsa-miR-3124-5p.fa",
            "ENST00000597389.fa",
            [
                Target(
                    host_seq="",
                    mirna_interval=Interval(start=3279, end=3300),
                    mirna_seq="TTCGCGGGCGAAGGCAAAGTC",
                    seed_interval=Interval(start=3293, end=3299),
                    seed_seq="TCGCGG",
                    seed_start_on_mirna=1,
                ),
                Target(
                    host_seq="",
                    mirna_interval=Interval(start=2794, end=2815),
                    mirna_seq="TTCGCGGGCGAAGGCAAAGTC",
                    seed_interval=Interval(start=2807, end=2814),
                    seed_seq="TCGCGGG",
                    seed_start_on_mirna=1,
                ),
                Target(
                    host_seq="",
                    mirna_interval=Interval(start=1994, end=2015),
                    mirna_seq="TTCGCGGGCGAAGGCAAAGTC",
                    seed_interval=Interval(start=2008, end=2014),
                    seed_seq="TCGCGG",
                    seed_start_on_mirna=1,
                ),
                Target(
                    host_seq="",
                    mirna_interval=Interval(start=530, end=551),
                    mirna_seq="TTCGCGGGCGAAGGCAAAGTC",
                    seed_interval=Interval(start=544, end=550),
                    seed_seq="TCGCGG",
                    seed_start_on_mirna=1,
                ),
            ],
        ),
        (
            "hsa-miR-3124-5p.fa",
            "ENST00000355526.fa",
            [
                Target(
                    host_seq="",
                    mirna_interval=Interval(start=0, end=21),
                    mirna_seq="TTCGCGGGCGAAGGCAAAGTC",
                    seed_interval=Interval(start=14, end=20),
                    seed_seq="TCGCGG",
                    seed_start_on_mirna=1,
                )
            ],
        ),
    ],
)
def test_target(path_root_test, fname_mirna, fname_transcript, expected_targets):
    _mirs = mirmap.utils.load_fasta(path_root_test.joinpath("data", fname_mirna))
    _mrnas = mirmap.utils.load_fasta(path_root_test.joinpath("data", fname_transcript))

    targets = mirmap.target.find_targets_with_seed(
        _mrnas[fname_transcript[:-3]], _mirs[fname_mirna[:-3]].upper().replace("U", "T")
    )
    assert len(expected_targets) == len(targets)

    for itarget, target in enumerate(targets):
        for att in ("mirna", "mirna_length", "mirna_seq", "seed", "seed_length", "seed_seq", "seed_start_on_mirna"):
            assert getattr(expected_targets[itarget], att) == getattr(target, att)
