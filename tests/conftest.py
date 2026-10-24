import pytest

from topwater import HullMesh
import json
from .reference_case import reference_inputs


@pytest.fixture
def input_data():
    return json.loads(reference_inputs().to_json())


@pytest.fixture
def mesh(input_data):
    return HullMesh(**input_data["body"]["mesh"])
