![GitHub License](https://img.shields.io/github/license/WargaSlowy/mega?style=for-the-badge)
[![](https://dcbadge.limes.pink/api/server/https://discord.gg/xkvjwsDrnx)](https://discord.gg/xkvjwsDrnx)

Python package for perform curated arithmetic calculation with __Cython__. this project was inspired with math (theory) and implementing on python library. this package using cython for speed up calculation and more optiomal and faster code.

## Build Package

> [!NOTE]  
> Make sure to install uv for build this package
> for information you can check [here](https://docs.astral.sh/uv/getting-started/installation/) for initialization

```sh
git clone https://github.com/mega
cd mega
```

build package
```sh
uv run setup.py build_ext --inplace;
uv pip install -e .
```

## Usage

```python
import mega.op.function as function

number: float = 5.2
print(function.gamma(number))
```
