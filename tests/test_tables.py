import runpy, pathlib
T = pathlib.Path(__file__).parent
def test_s1_table(): runpy.run_path(str(T / "test_s1.py"))
def test_s3_table(): runpy.run_path(str(T / "test_s3.py"))
