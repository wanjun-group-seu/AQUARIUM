# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#
import shlex
import subprocess
import tempfile
import sys

import pysrc.body.logger
from pysrc.body import utilities

__doc__ = ''' wrapper of subprocess, got a job , and run it
'''
__author__ = 'zerodel'

_logger = pysrc.body.logger.default_logger("external_command")


class ErrCmd(Exception):
    def __int__(self, process, description):
        self.process = process
        self.description = description

    def __str__(self):
        return "{when} @ {where}".format(when=self.description,
                                         where=str(self.process))


class Cmd(object):
    def __init__(self, cmd, pipe_in=None, env=None, use_shell=False, target_file=None, error_log=None):
        self.raw_cmd = cmd
        self.stdin = subprocess.PIPE
        self.stdout = None
        self.stderr = None
        self.use_shell = use_shell
        self.os_error = None
        self.env = env
        self.return_value = None
        self.use_tmp_file = False  # no , I do not want use tmp file right now.
        self.process = None
        self.pipe_in = pipe_in

        if isinstance(cmd, str):
            self.cmd_list = shlex.split(cmd)
        elif isinstance(cmd, list):
            self.cmd_list = cmd

        if len(self.cmd_list) == 0:
            raise ErrCmd(self, "Error: empty command")

        if error_log:
            self.stderr = error_log

        if self.use_tmp_file:
            self._tmp_stdout = tempfile.TemporaryFile(prefix="stdout_")
            self._tmp_stderr = tempfile.TemporaryFile(prefix="stderr_")

            self.stdout = self._tmp_stdout
            self.stderr = self._tmp_stderr

        else:
            self.stdout = subprocess.PIPE
            self.stderr = subprocess.PIPE

        if target_file:
            self.stdout = target_file

        self.is_process_end_and_output_loaded = False
        self.is_started = False

        _logger.debug("raw command line : %s" % self.raw_cmd)

    def __repr__(self):
        template_str = """input cmd={raw_cmd}
        return value={return_value}
        stdin={stdin}
        stdout={stdout}
        stderr={stderr}"""
        str_self = template_str.format(raw_cmd=str(self.raw_cmd),
                                       return_value=str(self.return_value),
                                       stdin=str(self.stdin),
                                       stdout=str(self.stdout),
                                       stderr=str(self.stderr)
                                       )

        return str_self.strip()

    def start(self):
        if not utilities.which(self.cmd_list[0]):
            raise FileExistsError("unable to find the first term of cli-string: %s" % self.cmd_list[0])

        try:
            self.process = subprocess.Popen(args=self.cmd_list,
                                            stdin=self.stdin,
                                            stdout=self.stdout,
                                            stderr=self.stderr,
                                            env=self.env,
                                            universal_newlines=True
                                            )
        except OSError as e:
            error_description = "Error: error on invoking :{}".format(self.raw_cmd)
            _logger.error(error_description)
            self.os_error = e
            raise ErrCmd(self, error_description)

        except ValueError as ve:
            err_description_value_error = "Error: InValid Arguments in {}".format(self.raw_cmd)
            _logger.error(err_description_value_error)
            self.os_error = ve
            raise ErrCmd(self, err_description_value_error)

        self.is_started = True
        return self

    def to(self, other_cmd):
        self.run()
        other_cmd.receive(self.stdout)
        return other_cmd

    def receive(self, message):
        if not self.is_process_end_and_output_loaded:
            self.pipe_in = message
        else:
            raise ErrCmd(self, "Error: Too late to send message")

    @staticmethod
    def _load_output_file(tmp_file):
        tmp_file.seek(0)
        return tmp_file.read()

    def _wait_process_to_end(self):
        if self.process and self.is_started:
            if self.use_tmp_file:
                self.process.wait()
                self.stdout = self._load_output_file(self._tmp_stdout)
                self.stderr = self._load_output_file(self._tmp_stderr)
            else:
                self.stdout, self.stderr = self.process.communicate(input=self.pipe_in)

            self.return_value = self.process.returncode
            self.is_process_end_and_output_loaded = True
        else:  # check self.process
            _logger.warning("Warning: process not running yet")

    def wait_until_finish(self):
        self._wait_process_to_end()

    def suicide(self):
        if self.process:
            if self.is_process_end_and_output_loaded:
                _logger.warning("Warning: process already finished")
            else:
                self.process.kill()
                _logger.warning("Warning: kill this process :{}".format(self.pid()))

        else:
            _logger.warning("Warning: try to suicide a non-existing process")

    def pid(self):
        if self.process:
            return self.process.pid

    def get_return_code(self):
        if isinstance(self.return_value, int):
            return self.return_value
        else:
            _logger.warning("Warning: process not finish yet")
            return None

    def run(self):
        if not self.is_started:
            self.start()

        if not self.is_process_end_and_output_loaded:
            self.wait_until_finish()

        print(self.stdout, file=sys.stdout)
        print(self.stderr, file=sys.stderr)

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.wait_until_finish()


def _run(cli_str):
    tin_can_cmd = Cmd(cli_str)
    tin_can_cmd.run()
    return tin_can_cmd


ERR_REPORT_WHEN_INVOKE_EXTERNAL_COMMAND = """
Error: this external command do not return 0 :
----start of external command----
{cmd}
----end of external command----

----stdout start----
{stdout}

----stdout ends ----

----stderr start----
{stderr}
----stderr ends----
"""


def run(cli_str):
    p = _run(cli_str)

    if p.get_return_code() != 0:
        raise ChildProcessError(ERR_REPORT_WHEN_INVOKE_EXTERNAL_COMMAND.format(cmd=cli_str,
                                                                               stdout=p.stdout,
                                                                               stderr=p.stderr
                                                                               ))


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print(__doc__)

    else:
        pass


