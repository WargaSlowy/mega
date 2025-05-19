#usr/bin/bash

make_libs () {
  echo "make project"
  uv run setup.py build_ext --inplace;
}

clean_project() {
  echo "clean project"
  find . -type d -name "__pycache__" -exec rm -rf {} +
}

if [ "$1" == "make" ]; then
  make_libs
elif [ "$1" == "clean" ]; then
  clean_project
fi

