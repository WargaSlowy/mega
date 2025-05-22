import pytest
import sys

def run_test() -> None:
    print("running testing")

    if len(sys.argv) > 1:
        test_files = sys.argv[1:]
    else:
        test_files = ["constant_test.py", "function_test.py"]

    return_code = pytest.main(test_files + ["-v"])

    if return_code == pytest.ExitCode.OK:
        print("\ntesting pass")
    else:
        print("\nfailed")

if __name__ == "__main__":
    run_test()
