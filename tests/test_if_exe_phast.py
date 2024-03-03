import platform

import pytest

import mirmap.if_exe_phast


@pytest.mark.unit()
def test_get_exact_prob(path_root_test):
    path_phylop = path_root_test.joinpath("..", "bin", f"{platform.system().lower()}_{platform.machine()}", "phyloP")
    result = mirmap.if_exe_phast.phylop(
        method="SPH",
        mode="CONACC",
        mod_fname=path_root_test.joinpath("data", "NM_024573.mod"),
        aln_fname=path_root_test.joinpath("data", "NM_024573_ts1.fa"),
        aln_format="FASTA",
        prune=True,
        path_exe=path_phylop,
    )
    assert 0.5303034 == pytest.approx(result)
