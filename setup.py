from setuptools import Extension, find_packages, setup
from Cython.Build import cythonize
import numpy as np
import os

ext = [
    Extension(
        "mega.utils.constant",
        ["mega/utils/constant.pyx"],
    ),
    Extension(
        "mega.op.function",
        ["mega/op/function.pyx"],
        include_dirs=[np.get_include()],
    ),
    Extension(
        "mega.op.arithmethic",
        ["mega/op/arithmethic.pyx"],
        include_dirs=[np.get_include()],
    ),
]

os.makedirs("mega/utils", exist_ok=True)
os.makedirs("mega/op", exist_ok=True)

setup(
    name="mega number utils",
    version="1.0",
    packages=find_packages(),
    ext_modules=cythonize(ext, language_level=3),
    include_dirs=[np.get_include()],
)
