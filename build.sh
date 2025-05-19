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
  echo "cleaning build folders"
  rm -rf build
}

if [ "$1" == "make" ]; then
  make_libs
elif [ "$1" == "clean" ]; then
  clean_project
elif [ "$1" == "package" ]; then
  make_as_package
fi

