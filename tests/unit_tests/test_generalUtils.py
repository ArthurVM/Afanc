import pytest
from pathlib import Path
from io import StringIO
from Afanc.utilities.generalUtils import isFile, isDir, checkDate
from Afanc.utilities.generalUtils import vprint
from Afanc.utilities.exceptions import FileNotFoundErrorAfanc, DirectoryNotFoundErrorAfanc, InvalidFileFormatError


def test_isFile_exists(tmp_path):
    """Test isFile when the file exists."""
    d = tmp_path / "sub"
    d.mkdir()
    p = d / "hello.txt"
    p.write_text("content")
    assert isFile(str(p)) == p.resolve()


def test_isFile_not_exists(tmp_path):
    """Test isFile when the file does not exist."""
    non_existent_file = tmp_path / "non_existent.txt"
    with pytest.raises(FileNotFoundErrorAfanc):
        isFile(str(non_existent_file))


def test_isDir_exists(tmp_path):
    """Test isDir when the directory exists."""
    d = tmp_path / "test_dir"
    d.mkdir()
    assert isDir(str(d)) == d.resolve()


def test_isDir_not_exists(tmp_path):
    """Test isDir when the directory does not exist."""
    non_existent_dir = tmp_path / "non_existent_dir"
    with pytest.raises(DirectoryNotFoundErrorAfanc):
        isDir(str(non_existent_dir))


@pytest.mark.parametrize(
    "date_str, expected_valid",
    [
        ("2022-05-01", True),
        ("2023-05-31", True),
        ("202-05-01", False),  ## invalid year
        ("2022-06-01", False), ## invalid month (must be 05)
        ("2022-05-1", False),  ## invalid day format
        ("2022/05/01", False), ## invalid separator
        ("invalid-date", False),
    ],
)
def test_checkDate(date_str, expected_valid):
    """Test checkDate with various date string formats."""
    if expected_valid:
        assert checkDate(date_str) == date_str
    else:
        with pytest.raises(InvalidFileFormatError):
            checkDate(date_str)


def test_vprint_tees_log_stream_to_console(capsys):
    log_stream = StringIO()

    vprint("TEST", "hello", "prYellow", log_stream)

    assert "TEST" in log_stream.getvalue()
    assert ":: hello" in log_stream.getvalue()
    captured = capsys.readouterr()
    assert "TEST" in captured.out
    assert ":: hello" in captured.out


def test_vprint_tees_error_log_stream_to_stderr(capsys):
    log_stream = StringIO()

    vprint("ERROR", "bad", "prRed", log_stream)

    captured = capsys.readouterr()
    assert "ERROR" in log_stream.getvalue()
    assert ":: bad" in log_stream.getvalue()
    assert "ERROR" in captured.err
    assert ":: bad" in captured.err
