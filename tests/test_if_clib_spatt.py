import platform

import pytest

import mirmap.if_lib_spatt


@pytest.mark.unit()
def test_get_exact_prob(path_root_test):
    path_libspatt = path_root_test.joinpath(
        "..", "bin", f"{platform.system().lower()}_{platform.machine()}", "libspatt2.so"
    )
    if_spatt = mirmap.if_lib_spatt.Spatt(path_libspatt)
    result = if_spatt.get_exact_prob(
        "AUUAAAA",
        1,
        1000,
        ["A", "U"],
        [[0.1780821917808219, 0.821917808219178], [0.6593406593406593, 0.34065934065934067]],
        1,
        "o",
    )
    assert 0.3704134091486 == pytest.approx(result)
