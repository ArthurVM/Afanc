import pytest

from Afanc.utilities.exceptions import CommandExecutionError
from Afanc.utilities.runCommands import command


def test_command_accepts_argv_input():
    stdout, stderr = command(["/bin/echo", "hello"], "TEST").run_comm_quiet(1)

    assert stdout.decode() == "hello\n"
    assert stderr.decode() == ""


def test_command_pipefail_raises_on_failed_pipeline():
    with pytest.raises(CommandExecutionError):
        command("false | true", "TEST", pipefail=True).run_comm_quiet(0)


def test_command_timeout_raises_cleanly():
    with pytest.raises(CommandExecutionError) as exc:
        command("sleep 2", "TEST", timeout=1).run_comm_quiet(0)

    assert "timed out" in str(exc.value).lower()
