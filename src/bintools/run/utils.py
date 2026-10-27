import os
import re
import shlex
import subprocess
import sys
from textwrap import dedent
from typing import Any
from typing import Dict
from typing import Optional
from typing import Tuple


def joint_kwargs(**kwargs) -> str:
    """
    builds "key value" pairs; values are converted to str and shell-quoted
    (so each value is a single argument, e.g. a path containing spaces).
    An empty value yields the bare key (flag without argument).
    """
    return " ".join(
        [
            (k + " " + shlex.quote(str(v))) if str(v) != "" else k
            for k, v in kwargs.items()
        ]
    )


def sub(
    cmd: str, cwd: Optional[str] = None, stdout: bool = False, stder: bool = False
) -> Tuple[Any, Any]:
    if cwd is not None:
        print(cmd + " " + cwd)
    else:
        print(cmd)
    # shell=False: redirections ("> f", "1> f", "2> f", ">> f") are not
    # interpreted by a shell, so they are extracted here and applied to the
    # child's stdout/stderr (they take precedence over the stdout/stder flags).
    cmd_ = []
    redirs = {1: None, 2: None}
    merge_stderr = False
    tokens = shlex.split(cmd)
    i = 0
    while i < len(tokens):
        m = re.fullmatch(r"([12]?)(>>?)(.*)", tokens[i])
        if tokens[i] == "2>&1":
            merge_stderr = True
        elif m is not None and (m.group(3) != "" or i + 1 < len(tokens)):
            if m.group(3) != "":
                path = m.group(3)
            else:
                i += 1
                path = tokens[i]
            redirs[int(m.group(1) or 1)] = (path, "a" if m.group(2) == ">>" else "w")
        else:
            cmd_.append(tokens[i])
        i += 1

    files = []
    try:
        streams: Dict[int, Any] = {
            1: subprocess.PIPE if stdout else None,
            2: subprocess.PIPE if stder else None,
        }
        for fd, spec in redirs.items():
            if spec is not None:
                fh = open(spec[0], spec[1])
                files.append(fh)
                streams[fd] = fh
        if merge_stderr:
            streams[2] = subprocess.STDOUT
        subp = subprocess.Popen(
            cmd_,
            stdout=streams[1],
            stderr=streams[2],
            shell=False,
            cwd=cwd,
        )
        return subp.communicate()
    finally:
        for fh in files:
            fh.close()


def run_shell_command(cmd, raise_errors=False, extra_env=None, cwd=None):
    """
    Run the given command string via Bash with error checking.
    Returns True if the command exits normally. Returns False if the command
    exits with failure and "raise_errors" is False (the default).  When
    "raise_errors" is True, exceptions are rethrown.
    If an *extra_env* mapping is passed, the provided keys and values are
    overlayed onto the default subprocess environment.
    """
    return ShellCommandRunner(
        cmd, raise_errors=raise_errors, extra_env=extra_env, cwd=cwd
    ).run()


class ShellCommandRunner:
    """
    Run the given command string via Bash with error checking.
    TODO move to method docstrings. Returns True if the command exits normally.  Returns False if the command
    exits with failure and "raise_errors" is False (the default).  When
    "raise_errors" is True, exceptions are rethrown.
    If an *extra_env* mapping is passed, the provided keys and values are
    overlayed onto the default subprocess environment.
    """

    def __init__(self, cmd, *, raise_errors=False, extra_env=None, cwd=None):
        self.cmd = cmd
        self.raise_errors = raise_errors
        self.extra_env = extra_env
        self.cwd = cwd

    def run(self):
        try:
            self.invoke_command()
        except Exception as error:
            self.print_error_message(error)

            if self.raise_errors:
                raise error

            return False

        return True

    def invoke_command(self):
        return subprocess.check_output(
            self.shell_executable + self.shell_args,
            shell=False,
            stderr=subprocess.STDOUT,
            env=self.modified_env,
            cwd=self.cwd,
        )

    @property
    def shell_executable(self):
        if os.name == "posix":
            return ["/bin/bash"]
        else:
            # We try best effort on other systems. For now that means nt/java.
            return ["env", "bash"]

    @property
    def shell_args(self):
        return ["-c", "set -euo pipefail; " + self.cmd]

    @property
    def modified_env(self):
        env = os.environ.copy()

        if self.extra_env:
            env.update(self.extra_env)

        return env

    def print_error_message(self, error):
        if isinstance(error, subprocess.CalledProcessError):
            message = f"{error.output}\nshell exited {error.returncode} when running: {self.cmd}"

            if error.returncode == 127:
                message += "\nAre you sure this program is installed?"
        elif isinstance(error, FileNotFoundError):
            shell = " and ".join(self.shell_executable)

            message = f"""
                Unable to run shell commands using {shell}!
                Augur requires {shell} to be installed.  Please open an issue on GitHub
                <https://github.com/nextstrain/augur/issues/new> if you need assistance.
                """
        else:
            message = str(error)

        self.print_error(message)

    @staticmethod
    def print_error(message):
        """Prints message to STDERR formatted with textwrap.dedent"""
        print("\nERROR: " + dedent(message).lstrip("\n") + "\n", file=sys.stderr)
