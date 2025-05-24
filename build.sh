#usr/bin/bash

make_libs () {
  echo "make project"
  uv run setup.py build_ext --inplace;
}

make_as_package() {
  make_libs
  uv pip install -e .
}

clean_project() {
  echo "clean project"
  find . -type d -name "__pycache__" -exec rm -rf {} +
  find . -type f -name "*.so" -exec rm -rf {} +
  find . -type f -name "*.c" -exec rm -rf {} +
  find . -type f -name "*.h" -exec rm -rf {} +
  echo "cleaning build folders egg info"
  rm -rf mega.egg-info build
}

testing_package() {
  make_as_package
  uv pip install .
  uv run pytest --verbose
}

if [ "$1" == "make" ]; then
  make_libs
elif [ "$1" == "clean" ]; then
  clean_project
elif [ "$1" == "package" ]; then
  make_as_package
elif [ "$1" == "testing" ]; then
  testing_package
fi

