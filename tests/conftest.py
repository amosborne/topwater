import pytest

from topwater import HullMesh
from .reference_hull import input_data as reference_input_data


@pytest.fixture
def input_data():
    return reference_input_data()


@pytest.fixture
def mesh(input_data):
    return HullMesh(**input_data["body"]["mesh"])
