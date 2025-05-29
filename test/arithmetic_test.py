import mega
import pytest
import math

SIGMA_VALUES: list = [
    (1, 0, 1),
    (6, 0, 4),
    (12, 0, 6),
    (28, 0, 6),
    (1, 1, 1),
    (6, 1, 1 + 2 + 3 + 6),
    (12, 1, 1 + 2 + 3 + 4 + 6 + 12),
    (28, 1, 1 + 2 + 4 + 7 + 14 + 28),
    (1, 2, 1),
    (6, 2, 1**2 + 2**2 + 3**2 + 6**2),
    (10, 2, 1**2 + 2**2 + 5**2 + 10**2),
    (100, 0, 9),
    (100, 1, 217),
]

SIGMA_COMPLEX_VALUES: list = [
    (
        6,
        0.5 + 1j,
        1 ** (0.5 + 1j) + 2 ** (0.5 + 1j) + 3 ** (0.5 + 1j) + 6 ** (0.5 + 1j),
    ),
    (6, 1 + 1j, 1 ** (1 + 1j) + 2 ** (1 + 1j) + 3 ** (1 + 1j) + 6 ** (1 + 1j)),
    (4, 0.5 + 1j, 1 ** (0.5 + 1j) + 2 ** (0.5 + 1j) + 4 ** (0.5 + 1j)),
]

EULER_PHI_VALUE: list = [
    (1, 1),
    (2, 1),
    (3, 2),
    (4, 2),
    (5, 4),
    (6, 2),
    (7, 6),
    (8, 4),
    (9, 6),
    (10, 4),
    (12, 4),
    (25, 20),
    (48, 16),
    (100, 40),
    (105, 48),
]

CHEBYSHEV_VALUE: list = [
    (2.0, math.log(2)),
    (3.0, math.log(2) + math.log(3)),
    (5.0, math.log(2) + math.log(3) + math.log(5)),
    (7.0, sum(math.log(p) for p in [2, 3, 5, 7])),
    (10.0, sum(math.log(p) for p in [2, 3, 5, 7])),
    (11.0, sum(math.log(p) for p in [2, 3, 5, 7, 11])),
    (20.0, sum(math.log(p) for p in [2, 3, 5, 7, 11, 13, 17, 19])),
]


def test_value_sigma() -> None:
    for n, z, expected in SIGMA_VALUES:
        sigma = mega.SigmaZ(n, z)
        result = sigma.compute()
        assert (
            result == expected
        ), f"SigmaZ({n}, {z}) -> {result}, expected = {expected}"


def test_value_complex_sigma() -> None:
    import numpy as np

    for n, z, expected in SIGMA_COMPLEX_VALUES:
        sigma = mega.SigmaZ(n, z)
        result = sigma.compute()

        assert isinstance(
            result, complex
        ), f"Expected complex output for SigmaZ({n}, {z}), got {type(result)}"

        assert np.isclose(result, expected, atol=1e-10), (
            f"SigmaZ({n}, {z}) → {result}, " f"expected {expected}"
        )


def test_perfect_square_sigma() -> None:
    sig = mega.SigmaZ(4, 1)
    assert sig.compute() == 1 + 2 + 4


def test_large_input_sigma() -> None:
    sig = mega.SigmaZ(1000, 0)
    assert sig.compute() == 16

    sig = mega.SigmaZ(1000, 1)
    assert sig.compute() == 2340


def test_value_euler_phi() -> None:
    for n, expected in EULER_PHI_VALUE:
        phi = mega.EulerPhi(n)
        result = phi.compute()
        assert (
            result == expected
        ), f"EulerPhi({n}).compute() -> {result}, expected {expected}"


def test_invalid_euler_phi() -> None:
    with pytest.raises(ValueError):
        mega.EulerPhi(0)


def test_case_prime_euler_phi() -> None:
    primes: list[int] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    for p in primes:
        phi = mega.EulerPhi(p)
        assert phi.compute() == p - 1, f"phi({p}) should -> {p - 1}"


def test_perfect_power_euler_phi() -> None:
    p: int = 2
    for k in range(1, 6):
        power = p**k
        expected = power - (p ** (k - 1))
        phi = mega.EulerPhi(power)
        assert phi.compute() == expected, f"phi({power}) should be {expected}"


def test_large_input_euler_phi() -> None:
    phi = mega.EulerPhi(1_000_000)
    assert phi.compute() == 400_000


def test_value_chebyshev() -> None:
    for x, expected in CHEBYSHEV_VALUE:
        ch = mega.Chebyshev(x)
        result = ch.compute()
        assert result == pytest.approx(expected, abs=1e-5), f"Chebyshev({x}).compute()"


def test_edge_cases_chebyshev() -> None:
    res = mega.Chebyshev(10.9)
    assert res.compute() == pytest.approx(
        sum(math.log(p) for p in [2, 3, 5, 7]), abs=1e-10
    )


def test_setitem_manual_cache_chebyshev() -> None:
    res = mega.Chebyshev(10.0)
    res[10] = 5.3471
    assert res[10] == 5.3471


def test_large_input_chebyshev() -> None:
    res = mega.Chebyshev(100.0)
    prime = [
        2,
        3,
        5,
        7,
        11,
        13,
        17,
        19,
        23,
        29,
        31,
        37,
        41,
        43,
        47,
        53,
        59,
        61,
        67,
        71,
        73,
        79,
        83,
        89,
        97,
    ]
    expected = sum(math.log(p) for p in prime)
    res_ch = res.compute()
    assert res_ch == pytest.approx(expected, abs=1e-5)
