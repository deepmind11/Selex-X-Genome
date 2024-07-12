import pytest
from ENCODE.experiment import TFChipSeq
from ENCODE.search import EncodeSearch, EncodeSearchError
from motifs.parse_motifcentral_json import MOTIFCENTRAL


@pytest.fixture
def my_valid_search():
    return EncodeSearch.search_using_MOTIFCENTRAL_json(MOTIFCENTRAL[222], "10")


@pytest.fixture
def my_invalid_search():
    return EncodeSearch.search_using_MOTIFCENTRAL_json(MOTIFCENTRAL[0])


def test_valid_search(my_valid_search):
    assert my_valid_search.tf == "CTCF"
    assert my_valid_search.organism == "Homo+sapiens"
    assert my_valid_search.limit == "10"
    assert my_valid_search.search_result is not None
    assert len(my_valid_search.search_result) == 10


def test_get_experiment_valid_search(my_valid_search):
    my_Experiment = my_valid_search.get_experiments()
    assert len(my_Experiment) == 10
    assert isinstance(my_Experiment[0], TFChipSeq)
    assert not my_Experiment[0].control


def test_invalid_search(my_invalid_search):
    assert my_invalid_search.tf == "TFE3"
    assert my_invalid_search.organism == "Mus+musculus"
    assert my_invalid_search.limit == "all"
    assert my_invalid_search.search_result is None


def test_get_experiment_invalid_search(my_invalid_search):
    with pytest.raises(EncodeSearchError):
        my_invalid_search.get_experiments()
