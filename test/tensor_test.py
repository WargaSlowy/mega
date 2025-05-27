import mega
import pytest

def test_dtype_support_tensor() -> None:
    t_int = mega.Tensor((2, 2), dtype="int")
    t_long = mega.Tensor((2, 2), dtype="long")
    t_float = mega.Tensor((2, 2), dtype="float")
    t_double = mega.Tensor((2, 2), dtype="double")

    assert t_int.dtype() == "int"
    assert t_long.dtype() == "long"
    assert t_float.dtype() == "float"
    assert t_double.dtype() == "double"

def test_invalid_dtype_value_error_tensor() -> None:
    with pytest.raises(ValueError):
        mega.Tensor((2, 2), dtype="str")

def test_invalid_dimension_size_tensor() -> None:
    with pytest.raises(ValueError):
        mega.Tensor((0, 2))
    with pytest.raises(ValueError):
        mega.Tensor((-1, 3))

def test_initial_data_int_tensor() -> None:
    data: list[int] = [1, 2, 3, 4]
    tensor = mega.Tensor((4,), data=data, dtype="int")
    assert tensor[0] == 1
    assert tensor[3] == 4

def test_initial_data_long_tensor() -> None:
    data: list[int] = [10, 20, 30]
    tensor = mega.Tensor((3,), data=data, dtype="long")
    assert tensor[1] == 20

def test_indexing_1d_tensor() -> None:
    tensor = mega.Tensor((5,), dtype="long")
    tensor[0] = 100
    tensor[4] = 200
    assert tensor[0] == 100
    assert tensor[4] == 200

def test_indexing_2d_tensor() -> None:
    tensor = mega.Tensor((2, 3), dtype="long")
    tensor[0, 0] = 1
    tensor[0, 0] = 1
    tensor[0, 1] = 2
    tensor[1, 2] = 6

    assert tensor[0, 0] == 1
    assert tensor[0, 1] == 2
    assert tensor[1, 2] == 6

def test_indexing_3d_tensor() -> None:
    tensor = mega.Tensor((2, 2, 2), dtype="long")
    tensor[0, 0, 0] = 1
    tensor[1, 1, 1] = 8

    assert tensor[0, 0, 0] == 1
    assert tensor[1, 1, 1] == 8

def test_out_of_bound_indexing_tensor() -> None:
    tensor = mega.Tensor((3,), dtype="long")
    with pytest.raises(IndexError):
        _ = tensor[3]

def test_repr_tensor() -> None:
    tensor = mega.Tensor((5,), dtype="long")
    for i in range(5):
        tensor[i] = i * 10
    repr_string = repr(tensor)
    assert "shape=(5)" in repr_string
    assert "dtype=long" in repr_string

    tensor_small_repr = mega.Tensor((2,), data=[1, 2], dtype="int")
    small_repr = repr(tensor_small_repr)
    assert "[1, 2]" in small_repr
