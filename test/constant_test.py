import pytest
import mega

CATALAN_VALUES: list[int] = [
    1,
    1,
    2,
    5,
    14,
    42,
    132,
    429,
    1430,
    4862,
    16796,
]

LUCAS_VALUES: dict = {
    0: 2,
    1: 1,
    2: 3,
    3: 4,
    4: 7,
    5: 11,
    6: 18,
    7: 29,
    8: 47,
    9: 76,
    10: 123,
    20: 15127,
    30: 1860498,
    40: 228826127,
}

GOLDEN_RATIO_VALUE: dict = {
    1: 2.0,
    2: 1.5,
    3: 1.6666666666666667,
    4: 1.6,
    5: 1.625,
    10: 1.6179775280898876,
    20: 1.6180339886754843,
    30: 1.6180339887498948,
}

TRUE_PHI: float = (1 + (5**0.5)) / 2


def test_scalar_input_catalan_number() -> None:
    for i, expected in enumerate(CATALAN_VALUES):
        assert mega.catalan_number(i) == expected


def test_negative_index_catalan_number() -> None:
    with pytest.raises(ValueError):
        mega.catalan_number(-1)


def test_known_values_lucas_number() -> None:
    for n, expected in LUCAS_VALUES.items():
        assert mega.lucas_number(n) == expected, f"lucas_number({n}) must be {expected}"


def test_zero_one_lucas_number() -> None:
    assert mega.lucas_number(0) == 2
    assert mega.lucas_number(1) == 1


def test_large_value_lucas_number() -> None:
    assert mega.lucas_number(20) == 15127
    assert mega.lucas_number(30) == 1860498
    assert mega.lucas_number(40) == 228826127


def test_value_golden_ratio() -> None:
    for iterations, expected in GOLDEN_RATIO_VALUE.items():
        result = mega.golden_ratio(iterations)
        assert (
            abs(result - expected) < 1e-3
        ), f"golden_ratio({iterations}) = {result}, expected {expected}"


def test_coverage_to_phi_golden_ratio() -> None:
    result = mega.golden_ratio()
    assert abs(result - TRUE_PHI) < 1e-3


def test_return_type_golden_ratio() -> None:
    result = mega.golden_ratio()
    assert isinstance(result, float)
