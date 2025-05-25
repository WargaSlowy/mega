import math
import mega
import pytest

HAVERSINE_VALUE: dict = {
    0.0: 0.0,
    math.pi / 2: (1 - 0) / 2 == 0.5,  
    math.pi: (1 - (-1)) / 2 == 1.0,  
    3 * math.pi / 2: (1 - 0) / 2 == 0.5,
    2 * math.pi: (1 - 1) / 2 == 0.0,
}

GAMMA_VALUE: list = [
    (0.5, math.sqrt(math.pi)),
    (1.0, 1.0),
    (2.0, 1.0),
    (3.0, 2.0),
    (4.0, 6.0), 
    (5.0, 24.0),
]

JORDAN_TOTIEN_VALUE: list = [
    (6, 0, 0),
    (12, 0, 0),
    (1, 1, 1),
    (2, 1, 1),
    (3, 1, 2),
    (5, 1, 4),
    (6, 1, 2),
    (7, 1, 6),
    (10, 1, 4),
    (1, 2, 1),
    (2, 2, 3),
    (3, 2, 8),
    (5, 2, 24),
    (6, 2, 24),
    (1, 3, 1),
    (2, 3, 7),
    (3, 3, 26),
    (5, 3, 124),
    (6, 3, 182),
]

def test_small_theta_haversine() -> None:
    theta: float = 0.001
    expected: float = (1.0 - math.cos(theta)) / 2.0
    result = mega.Haversine(theta).compute()
    assert abs(result - expected) < 1e-3

def test_negative_angle_haversine() -> None:
    theta: float = -math.pi / 2
    result_neg: float = mega.Haversine(theta).compute()
    result_pos: float = mega.Haversine(-theta).compute()
    assert abs(result_neg - result_pos) < 1e-12

def test_haversine_with_earth_distance_formula() -> None:
    def earth_distance(lat1, lon1, lat2, lon2):
        R: float = 6371.0
        dlat = math.radians(lat2 - lat1)
        dlon = math.radians(lon2 - lon1)
        a = mega.Haversine(dlat).compute() + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * mega.Haversine(dlon).compute()
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
        return R * c

    dist = earth_distance(40.7128, -74.0060, 40.7128, -74.0060)
    assert dist < 1e-9, f"distance should be zero, got {dist} km"

def test_value_gamma() -> None:
    for point, expected in GAMMA_VALUE:
        result = mega.Gamma(point).compute()
        assert abs(result - expected) < 1e-5, (
            f"gamma({point}) = {result}, expected {expected}"
        )

def test_integer_input_gamma() -> None:
    assert mega.Gamma(6).compute() == pytest.approx(120.0, rel=1e-3)
    assert mega.Gamma(7).compute() == pytest.approx(720.0, rel=1e-10)

def test_reflection_gamma() -> None:
    z: float = 1.0 / 3.0
    expected = math.pi / math.sin(math.pi * z)
    result = mega.Gamma(z).compute() * mega.Gamma(1.0 - z).compute()
    assert abs(result - expected) < 1e-10

def test_value_jordan_totient() -> None:
    for n, k, expected in JORDAN_TOTIEN_VALUE:
        result = mega.JordanTotient(n, k).compute()
        assert result == expected, f"JordanTotient({n}, {k}) = {result}, expected = {expected}"

def test_k_zero_jordan_totient() -> None:
    for n in range(1, 21):
        assert mega.JordanTotient(n, 0).compute() == 0

def test_large_input_jordan_totient() -> None:
    assert mega.JordanTotient(100, 2).compute() == 5184
    assert mega.JordanTotient(12, 2).compute() == 72
