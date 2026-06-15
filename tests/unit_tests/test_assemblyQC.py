import pytest
from unittest.mock import patch, mock_open, MagicMock
import numpy as np
import pandas as pd
from pathlib import Path

from Afanc.autodatabase.assemblyQC import mash, buildMatrix, fastaMove

@pytest.fixture
def mock_args():
    """Fixture to create mock arguments."""
    args = MagicMock()
    args.threads = 4
    args.stdout = MagicMock() # Mock stdout file object
    args.stderr = MagicMock() # Mock stderr file object
    args.fasta_WDir = "/fake/fasta_WDir"
    args.cleanFasta_WDir = "/fake/cleanFasta_WDir"
    return args

@patch("Afanc.autodatabase.assemblyQC.command")
def test_mash_command_execution(mock_command_class, mock_args, tmp_path, monkeypatch):
    """Test that the mash function calls the 'mash' command correctly."""
    monkeypatch.chdir(tmp_path)
    mock_cmd_instance = MagicMock()
    mock_command_class.return_value = mock_cmd_instance

    taxon_id = "12345"
    fastas = [str(tmp_path / "f1.fasta"), str(tmp_path / "f2.fasta")]
    
    # Create dummy fasta files if their content or existence is checked by mash sketch
    for f_path in fastas:
        Path(f_path).touch()

    expected_mashdist_out = str(tmp_path / f"{taxon_id}_mashdist.txt")

    result_path = mash(mock_args, taxon_id, fastas)

    assert result_path == expected_mashdist_out
    
    expected_sketch_line = f"mash sketch -o ref {fastas[0]} {fastas[1]}"
    expected_dist_line = f"mash dist ref.msh ref.msh > {expected_mashdist_out}"

    assert mock_command_class.call_count == 2
    mock_command_class.assert_any_call(expected_sketch_line, "MASH")
    mock_command_class.assert_any_call(expected_dist_line, "MASH")
    
    assert mock_cmd_instance.run_comm_quiet.call_count == 2

    mock_cmd_instance.run_comm_quiet.assert_any_call(0, mock_args.stdout, mock_args.stderr)


@patch("Afanc.autodatabase.assemblyQC.move") # Patch shutil.move where it's used in assemblyQC
@patch("Afanc.autodatabase.assemblyQC.np.savetxt") # Still called to save _mash.txt
def test_buildMatrix_fewer_than_3_assemblies(mock_numpy_savetxt, mock_shutil_move, mock_args, tmp_path, monkeypatch):
    """Test buildMatrix when there are fewer than 3 assemblies."""
    monkeypatch.chdir(tmp_path)
    # Setup for 2 assemblies
    fasta1_path = tmp_path / "asm1.fa"
    fasta1_path.touch()
    fasta2_path = tmp_path / "asm2.fa"
    fasta2_path.touch()

    # Mash output content for 2 assemblies
    mash_output_content = (
        f"{fasta1_path}\t{fasta1_path}\t0.0\t0\t1000/1000\n"
        f"{fasta1_path}\t{fasta2_path}\t0.1\t0\t900/1000\n"
        f"{fasta2_path}\t{fasta1_path}\t0.1\t0\t900/1000\n"
        f"{fasta2_path}\t{fasta2_path}\t0.0\t0\t1000/1000\n"
    )
    inMash_path = tmp_path / "123_mashdist.txt"
    inMash_path.write_text(mash_output_content)

    # The function writes output files to the current working directory.
    # For cleanup, we'll track these files.
    expected_warning_file_path = tmp_path / "123_warning.txt"
    expected_mash_list_path = tmp_path / "123_mash.txt"

    try:
        calcArray, tax, modeVal = buildMatrix(mock_args, str(inMash_path))

        assert calcArray is None
        assert tax is None
        assert modeVal is None

        mock_numpy_savetxt.assert_called_once() # _mash.txt is still created

        # Assert shutil.move was called correctly
        assert mock_shutil_move.call_count == 2 # For fasta1_path and fasta2_path
        mock_shutil_move.assert_any_call(str(fasta1_path), mock_args.cleanFasta_WDir)
        mock_shutil_move.assert_any_call(str(fasta2_path), mock_args.cleanFasta_WDir)

        # Check warning file was created and has expected content
        assert expected_warning_file_path.exists()
        assert "Warning: Taxon 123 has fewer than 3 samples." in expected_warning_file_path.read_text()

    finally:
        # Clean up files created in CWD
        if expected_warning_file_path.exists():
            expected_warning_file_path.unlink()
        if expected_mash_list_path.exists():
            expected_mash_list_path.unlink()


@patch("Afanc.autodatabase.assemblyQC.np.savetxt")
def test_buildMatrix_sufficient_assemblies(mock_numpy_savetxt, mock_args, tmp_path, monkeypatch):
    """Test buildMatrix with 3 or more assemblies, exercising the main QC logic."""
    monkeypatch.chdir(tmp_path)
    # Setup for 3 assemblies
    fasta1_path = tmp_path / "asmA.fa"
    fasta1_path.touch()
    fasta2_path = tmp_path / "asmB.fa"
    fasta2_path.touch()
    fasta3_path = tmp_path / "asmC.fa"
    fasta3_path.touch()

    mash_output_content = (
        f"{fasta1_path}\t{fasta1_path}\t0.0\t0\t1000/1000\n"
        f"{fasta1_path}\t{fasta2_path}\t0.1\t0\t900/1000\n"  # asmA vs asmB
        f"{fasta1_path}\t{fasta3_path}\t0.2\t0\t800/1000\n"  # asmA vs asmC
        f"{fasta2_path}\t{fasta1_path}\t0.1\t0\t900/1000\n"  # asmB vs asmA
        f"{fasta2_path}\t{fasta2_path}\t0.0\t0\t1000/1000\n"
        f"{fasta2_path}\t{fasta3_path}\t0.05\t0\t950/1000\n" # asmB vs asmC
        f"{fasta3_path}\t{fasta1_path}\t0.2\t0\t800/1000\n"  # asmC vs asmA
        f"{fasta3_path}\t{fasta2_path}\t0.05\t0\t950/1000\n" # asmC vs asmB
        f"{fasta3_path}\t{fasta3_path}\t0.0\t0\t1000/1000\n"
    )
    inMash_path = tmp_path / "456_mashdist.txt"
    inMash_path.write_text(mash_output_content)
    
    expected_mash_list_path = tmp_path / "456_mash.txt"

    try:
        calcArray, tax, mode_val = buildMatrix(mock_args, str(inMash_path))

        assert tax == ['456']
        # avArray for asmA: (0.0+0.1+0.2)/3 = 0.1
        # avArray for asmB: (0.1+0.0+0.05)/3 = 0.05
        # avArray for asmC: (0.2+0.05+0.0)/3 = 0.08333...
        # rdArray (rounded): [0.1, 0.05, 0.0833] (approx)
        # mode(rdArray)[0] (if unique, smallest): np.array([0.05])
        assert np.isclose(np.atleast_1d(mode_val)[0], 0.05)

        # calcArray is fileavArray after np.unique on average distances, sorted by unique avg_dist
        # Expected order: asmB (0.05), asmC (0.0833), asmA (0.1)
        assert calcArray.shape == (3, 2)
        assert calcArray[0, 0] == str(fasta2_path) # asmB
        assert np.isclose(float(calcArray[0, 1]), 0.05)
        assert calcArray[1, 0] == str(fasta3_path) # asmC
        assert np.isclose(float(calcArray[1, 1]), 0.0833333333)
        assert calcArray[2, 0] == str(fasta1_path) # asmA
        assert np.isclose(float(calcArray[2, 1]), 0.1)

        mock_numpy_savetxt.assert_called_once()
    finally:
        if expected_mash_list_path.exists():
            expected_mash_list_path.unlink()
