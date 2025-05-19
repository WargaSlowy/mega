#usr/bin/bash

make_libs () {
  uv run setup.py build_ext --inplace;
}

clean_project() {
  find . -type f $ -name "*.so" -o -name "*.c" -o -name "*.pyc" $ -exec rm -f {} +
find . -type d -name "__pycache__" -exec rm -rf {} +
}

if [ "$1" == "make" ]; then
  make_libs
elif [ "$1" == "clean" ]; then
  clean_project
fi

