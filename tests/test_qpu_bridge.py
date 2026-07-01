#!/usr/bin/env python3
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#

# -*- coding: utf-8 -*-

"""Tests for the cross-cluster QPU bridge (Qdislib.utils.qpu_bridge).

qililab (the quantum-control library) is not installable off-cluster, so the
listener's hardware execution is not exercised here; instead the file-based
hand-off protocol is tested against a stand-in listener, plus the timeout and
error-propagation paths and the runcard resolution.
"""

import glob
import os
import threading
import time
import unittest

from tests import BaseTimedTestCase

from Qdislib.utils.qpu_bridge import (
    submit_and_wait,
    subcircuits_dir,
    _resolve_runcard,
    _atomic_write_pickle,
    _CIRCUIT_PREFIX,
    _RESULT_PREFIX,
    _ERROR_KEY,
    QPUTimeoutError,
    QPUExecutionError,
)


class _FakeResult:
    """Stand-in for a qililab result object (module-level, so it pickles)."""

    def __init__(self, counts):
        self._counts = counts

    def counts(self):
        return self._counts


def _fake_listener(directory, payload, deadline=5.0):
    """Minimal listener: turn the first circuit_*.pkl into a result_*.pkl."""
    end = time.monotonic() + deadline
    while time.monotonic() < end:
        pending = glob.glob(os.path.join(directory, f"{_CIRCUIT_PREFIX}*.pkl"))
        if pending:
            result_path = pending[0].replace(_CIRCUIT_PREFIX, _RESULT_PREFIX)
            _atomic_write_pickle(payload, result_path)
            os.remove(pending[0])
            return
        time.sleep(0.02)


class SubcircuitsDirTest(BaseTimedTestCase):
    def test_env_override_and_default(self):
        from Qdislib.utils import qpu_bridge

        saved = os.environ.pop("QDISLIB_SUBCIRCUITS_DIR", None)
        try:
            self.assertEqual(subcircuits_dir(), qpu_bridge.DEFAULT_SUBCIRCUITS_DIR)
            os.environ["QDISLIB_SUBCIRCUITS_DIR"] = "/tmp/whatever"
            self.assertEqual(subcircuits_dir(), "/tmp/whatever")
        finally:
            os.environ.pop("QDISLIB_SUBCIRCUITS_DIR", None)
            if saved is not None:
                os.environ["QDISLIB_SUBCIRCUITS_DIR"] = saved


class ResolveRuncardTest(BaseTimedTestCase):
    def test_resolution(self):
        # Explicit path always wins.
        self.assertEqual(_resolve_runcard("/rc.yml", "red"), "/rc.yml")

        saved_red = os.environ.pop("RUNCARD_RED", None)
        saved_blue = os.environ.pop("RUNCARD_BLUE", None)
        try:
            os.environ["RUNCARD_RED"] = "/red.yml"
            os.environ["RUNCARD_BLUE"] = "/blue.yml"
            self.assertEqual(_resolve_runcard(None, "red"), "/red.yml")
            self.assertEqual(_resolve_runcard(None, "blue"), "/blue.yml")
            with self.assertRaises(ValueError):
                _resolve_runcard(None, "green")  # unknown QPU
            os.environ.pop("RUNCARD_RED")
            with self.assertRaises(ValueError):
                _resolve_runcard(None, "red")  # selected but env unset
        finally:
            for key, val in (("RUNCARD_RED", saved_red), ("RUNCARD_BLUE", saved_blue)):
                os.environ.pop(key, None)
                if val is not None:
                    os.environ[key] = val


class SubmitAndWaitTest(BaseTimedTestCase):
    def test_roundtrip_with_stand_in_listener(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            payload = _FakeResult({"00": 6, "11": 4})
            listener = threading.Thread(target=_fake_listener, args=(tmp, payload))
            listener.start()
            result = submit_and_wait(
                "unmeasured-circuit", directory=tmp, timeout=5, poll_interval=0.02
            )
            listener.join()
            self.assertEqual(result.counts(), {"00": 6, "11": 4})
            # Both the circuit and result files are cleaned up.
            self.assertEqual(os.listdir(tmp), [])

    def test_timeout_when_no_listener(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(QPUTimeoutError):
                submit_and_wait("c", directory=tmp, timeout=0.1, poll_interval=0.02)
            # A timed-out submission leaves no orphan files behind.
            self.assertEqual(os.listdir(tmp), [])

    def test_listener_error_is_propagated(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            payload = {_ERROR_KEY: "RuntimeError: qpu offline"}
            listener = threading.Thread(target=_fake_listener, args=(tmp, payload))
            listener.start()
            with self.assertRaises(QPUExecutionError):
                submit_and_wait("c", directory=tmp, timeout=5, poll_interval=0.02)
            listener.join()


if __name__ == "__main__":
    unittest.main()
