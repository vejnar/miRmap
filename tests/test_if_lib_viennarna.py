import pytest

import mirmap.if_lib_viennarna


@pytest.mark.unit()
def test_fold():
    result = mirmap.if_lib_viennarna.RNAfold(
        ["AUCGAUGCGAUCGAGGGGCGCCCUUAAAGCUCUGAGGCGGCCCCCCA"],
        partfunc=True,
    )
    assert -20.20 == pytest.approx(result.mfe)
