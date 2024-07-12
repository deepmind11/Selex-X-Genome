# Experiment:ENCSR068HEE  Control:ENCSR608IVH
import json
import sqlite3

import pytest
from ENCODE.base import RunType
from ENCODE.experiment import Control, TFChipSeq
from ENCODE.library import Library


@pytest.fixture
def my_TFChipSeq():
    return TFChipSeq("ENCSR068HEE")


def test_TFChipSeq_Initialization(my_TFChipSeq):
    assert my_TFChipSeq.accession == "ENCSR068HEE"
    assert my_TFChipSeq.expr_data is not None
    assert not my_TFChipSeq.control


def test_get_url(my_TFChipSeq):
    assert my_TFChipSeq.get_url() == "https://www.encodeproject.org/ENCSR068HEE"


def test_fetchData(my_TFChipSeq):
    my_TFChipSeq.fetchData()
    assert my_TFChipSeq.expr_data is not None


def test_get_run_type(my_TFChipSeq):
    my_TFChipSeq.fetchData()
    TFChipSeq.get_run_type(my_TFChipSeq.expr_data) == RunType.SE


def test_get_controls(my_TFChipSeq):
    my_Control = my_TFChipSeq.get_controls()
    assert my_Control.accession == "ENCSR608IVH"
    assert my_Control.expr_data is not None
    assert my_Control.control
    assert isinstance(my_Control, Control)


def test_get_libraries(my_TFChipSeq):
    my_Libraries = my_TFChipSeq.get_libraries()
    assert len(my_Libraries) == 4

    expected_libraries = {
        "ENCLB386VDZ",
        "ENCLB191KSM",
        "ENCLB534KAZ",
        "ENCLB951UFA",
    }

    actual_libraries = {library.accession for library in my_Libraries}

    assert expected_libraries == actual_libraries

    for library in my_Libraries:
        assert isinstance(library, Library)


def test_update_db(my_TFChipSeq):
    my_TFChipSeq.update_database()

    with sqlite3.connect(
        "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/database/Selex_X_Genome.db"
    ) as conn:

        cursor = conn.cursor()

        cursor.execute(
            """
            SELECT "control" FROM "experiments" WHERE "accession" =?
            """,
            (my_TFChipSeq.accession,),
        )
        control = cursor.fetchone()[0]

    assert control == "n"
