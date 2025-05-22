import pytest
import mega.utils.constant as constant

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

def test_scalar_input_catalan_number() -> None:
    for i, expected in enumerate(CATALAN_VALUES):
        assert constant.catalan_number(i) == expected

def test_negative_index_catalan_number() -> None:
    with pytest.raises(ValueError):
        constant.catalan_number(-1)

def test_known_values_lucas_number() -> None:
    for n, expected in LUCAS_VALUES.items():
        assert constant.lucas_number(n) == expected, f"lucas_number({n}) must be {expected}"

def test_zero_one_lucas_number() -> None:
    assert constant.lucas_number(0) == 2
    assert constant.lucas_number(1) == 1

def test_large_value_lucas_number() -> None:
    assert constant.lucas_number(20) == 15127
    assert constant.lucas_number(30) == 1860498
    assert constant.lucas_number(40) == 228826127
