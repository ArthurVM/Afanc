import sys
import shlex
from os import kill, path
from signal import alarm, signal, SIGALRM, SIGKILL
from subprocess import PIPE, Popen

from .exceptions import CommandExecutionError
from .generalUtils import vprint


class Error(Exception):
    """ Base class for exceptions in this module.
    """
    pass


class InputError(Error):
    """ Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message

        sys.stderr.write(f"INPUT ERROR : {self.expression}\n{self.message}")
        sys.exit(1)


class command():
    """A class for handling shell executed commands."""
    def __init__(
        self,
        command,
        subprocessID,
        shell=None,
        shell_executable=None,
        pipefail=False,
        cwd=None,
        env=None,
        timeout=360000,
    ):
        self.command = command
        self.subprocessID = str(subprocessID)
        self.shell_executable = shell_executable
        self.pipefail = pipefail
        self.cwd = cwd
        self.env = env
        self.timeout = timeout

        ## infer shell use from the command type
        if shell is None:
            self.shell = not isinstance(command, (list, tuple))
        else:
            self.shell = shell

        ## pipefail requires bash
        if self.pipefail:
            self.shell = True
            if self.shell_executable == None:
                self.shell_executable = "/bin/bash"

        self.command_display = self._format_command_for_display(command)

    def run(self, timeout = -1):

        kill_tree = True

        class Alarm(Exception):
            pass

        def alarm_handler(signum, frame):
            raise Alarm

        popen_command = self._build_popen_command()
        popen_kwargs = {
            "shell": self.shell,
            "stdout": PIPE,
            "stderr": PIPE,
            "env": self.env,
            "cwd": self.cwd,
        }

        if self.shell and self.shell_executable != None:
            popen_kwargs["executable"] = self.shell_executable

        p=Popen(popen_command, **popen_kwargs)

        if timeout != -1:
            signal(SIGALRM, alarm_handler)
            alarm(timeout)

        try:
            stdout, stderr = p.communicate()
            if timeout != -1:
                alarm(0)

        except Alarm as a:
            pids = [p.pid]
            if kill_tree:
                pids.extend(self.get_process_children(p.pid))
            for pid in pids:
                ## process may exit before cleanup
                try:
                    kill(pid, SIGKILL)
                except OSError:
                    pass
            stderr = f"Command timed out after {timeout} seconds.".encode()
            return -1, b'', stderr

        return p.returncode, stdout, stderr

    def run_comm(self, if_out_return, so=None, se=None, raise_on_error=True):
        """ Run the command with or without an exit code
        """
        return self._run_comm(
            if_out_return,
            so=so,
            se=se,
            raise_on_error=raise_on_error,
            quiet=False,
        )

    def run_comm_quiet(self, if_out_return, so=None, se=None, raise_on_error=True):
        """ Run the command quietly with or without an exit code
        """
        return self._run_comm(
            if_out_return,
            so=so,
            se=se,
            raise_on_error=raise_on_error,
            quiet=True,
        )

    def _run_comm(self, if_out_return, so=None, se=None, raise_on_error=True, quiet=False):
        """Shared implementation for visible and quiet command execution."""
        if not quiet:
            print(f"COMMAND={self.command_display}")

        self._write_command_headers(so=so, se=se)

        returncode, stdout, stderr = self.run(self.timeout)

        stdout_decoded = stdout.decode(errors='replace')
        stderr_decoded = stderr.decode(errors='replace')

        self._write_command_output(stdout_decoded, stderr_decoded, so=so, se=se)

        if returncode != 0:
            self._handle_error(returncode, stdout_decoded, stderr_decoded, raise_on_error, quiet)

        if if_out_return:
            return stdout, stderr

    def _write_command_headers(self, so=None, se=None):
        if so != None:
            print(f"{self.subprocessID}\nCOMMAND={self.command_display}\n", file=so, flush=True)
        if se != None:
            print(f"{self.subprocessID}\nCOMMAND={self.command_display}\n", file=se, flush=True)

    def _write_command_output(self, stdout_decoded, stderr_decoded, so=None, se=None):
        if so != None:
            print(f"{stdout_decoded}\n", file=so, flush=True)
        if se != None:
            print(f"{stderr_decoded}\n", file=se, flush=True)

    def _handle_error(self, returncode, stdout_decoded, stderr_decoded, raise_on_error, quiet):
        if raise_on_error:
            raise CommandExecutionError(self.command_display, returncode, stdout_decoded, stderr_decoded)

        if quiet:
            print(f"Quiet command '{self.command_display}' failed. STDOUT:\n{stdout_decoded}\nSTDERR:\n{stderr_decoded}\nContinuing as raise_on_error is False.")
            return

        vprint(self.subprocessID, f"Command '{self.command_display}' failed. Continuing as raise_on_error is False.", "prRed")
        print(f"STDOUT :\n{stdout_decoded}\nSTDERR :\n{stderr_decoded}\n")

    def _build_popen_command(self):
        if self.pipefail:
            command_str = self._format_command_for_shell(self.command)
            return f"set -o pipefail; {command_str}"

        if self.shell:
            return self._format_command_for_shell(self.command)

        return self.command

    def _format_command_for_shell(self, command):
        if isinstance(command, (list, tuple)):
            return shlex.join([str(c) for c in command])
        return str(command)

    def _format_command_for_display(self, command):
        if isinstance(command, (list, tuple)):
            return shlex.join([str(c) for c in command])
        return str(command)

    def get_process_children(self, pid):
        p = Popen(f"ps --no-headers -o pid --ppid {pid}", shell = True,
	         stdout = PIPE, stderr = PIPE)
        stdout, stderr = p.communicate()
        return [int(p) for p in stdout.split()]
