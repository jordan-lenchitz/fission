import fission_core

def test_init():
    try:
        result = fission_core.sum_as_string(5, 7)
        print(f"Fission Core loaded! Result of sum_as_string(5, 7): {result}")
        if result == "12":
            print("QA Checkpoint Passed: Rust and Python are talking!")
        else:
            print(f"QA Checkpoint Failed: Unexpected result {result}")
    except Exception as e:
        print(f"QA Checkpoint Failed: Could not call fission_core function. Error: {e}")

if __name__ == "__main__":
    test_init()
