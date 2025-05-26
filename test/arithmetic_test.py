import mega
import pytest

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

def test_value_sigma() -> None:
    for n, z, expected in SIGMA_VALUES:
        sigma = mega.SigmaZ(n, z)
        result = sigma.compute()
        assert result == expected, f"SigmaZ({n}, {z}) -> {result}, expected = {expected}"

def test_invalid_sigma() -> None:
    with pytest.raises(ValueError):
        mega.SigmaZ(0, 1)
    with pytest.raises(ValueError):
        mega.SigmaZ(-5, 2)

def test_invalid_z_sigma() -> None:
    with pytest.raises(ValueError):
        mega.SigmaZ(5, -1)

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
        assert result == expected, f"EulerPhi({n}).compute() -> {result}, expected {expected}"

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
        power = p ** k
        expected = power - (p ** (k - 1))
        phi = mega.EulerPhi(power)
        assert phi.compute() == expected, f"phi({power}) should be {expected}"

def test_large_input_euler_phi() -> None:
    phi = mega.EulerPhi(1_000_000)
    assert phi.compute() == 400_000
