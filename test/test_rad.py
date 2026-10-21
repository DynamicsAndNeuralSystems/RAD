import csv
import math
import pathlib
import sys
import warnings

import numpy as np
import pytest

HERE = pathlib.Path(__file__).parent
sys.path.insert(0, str(HERE.parent))
from RAD import RAD  # noqa: E402

SERIES = {row[0]: np.array(row[1:], dtype=float)
          for row in csv.reader(open(HERE / "series.csv"))}
EXPECTED = list(csv.DictReader(open(HERE / "expected.csv")))


@pytest.mark.parametrize("row", EXPECTED, ids=lambda r: f"{r['series']}-tau{r['tau']}-c{r['doCentre']}")
def test_matches_reference(row):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        f = RAD(SERIES[row["series"]], int(row["tau"]), bool(int(row["doCentre"])))
    expected = float(row["RAD"])
    if math.isnan(expected):
        assert math.isnan(f)
    else:
        assert f == pytest.approx(expected, rel=1e-12)


def test_defaults_are_tau_1_centred():
    x = SERIES["ar1"]
    assert RAD(x) == RAD(x, 1, True)


def test_degenerate_input_warns():
    with pytest.warns(RuntimeWarning):
        assert math.isnan(RAD(SERIES["constant"]))


def test_row_and_column_vectors():
    x = SERIES["ar1"]
    assert RAD(x.reshape(1, -1)) == RAD(x) == RAD(x.reshape(-1, 1))


def test_old_argument_order_is_rejected():
    with pytest.raises(TypeError):
        RAD(SERIES["ar1"], True, 2)
    with pytest.raises(TypeError):
        RAD(SERIES["ar1"], 1, 2)
    with pytest.raises(ValueError):
        RAD(SERIES["ar1"], 0)
