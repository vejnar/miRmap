import pytest


@pytest.fixture(scope="module")
def path_root_test(request):
    """Return the directory of the currently running test script."""
    return request.path.parent
