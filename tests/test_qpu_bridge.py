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
import importlib.util
import logging
import os
import sys
import threading
import time
import types
import unittest
from unittest import mock

from tests import BaseTimedTestCase

from Qdislib.utils.qpu_bridge import (
    run_qpu_listener,
    submit_and_wait,
    subcircuits_dir,
    _resolve_runcard,
    _atomic_write_pickle,
    _process_one,
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
    def test_the_directory_is_required(self):
        saved = os.environ.pop("QDISLIB_SUBCIRCUITS_DIR", None)
        try:
            with self.assertRaisesRegex(ValueError, "QDISLIB_SUBCIRCUITS_DIR"):
                subcircuits_dir()
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


class ListenerDirectoryTest(BaseTimedTestCase):
    def test_a_directory_with_glob_characters(self):
        import tempfile
        from unittest import mock

        from Qdislib.utils import qpu_bridge

        qililab = types.ModuleType("qililab")
        qililab.logger = logging.getLogger("qililab-stand-in")
        qililab.build_platform = lambda runcard: mock.MagicMock()
        digital = types.ModuleType("qililab.digital")
        digital.DigitalTranspilationConfig = lambda routing: None
        seen = []

        def process(platform, circuit_path, nshots, transpilation):
            seen.append(os.path.basename(circuit_path))
            os.remove(circuit_path)

        with tempfile.TemporaryDirectory(prefix="run[1]_") as tmp:
            open(os.path.join(tmp, "circuit_1.pkl"), "wb").close()
            with (
                mock.patch.dict(
                    sys.modules, {"qililab": qililab, "qililab.digital": digital}
                ),
                mock.patch.object(qpu_bridge, "_process_one", process),
            ):
                qpu_bridge.run_qpu_listener(
                    tmp, runcard="/r.yml", poll_interval=0.01, idle_timeout=0.05
                )
        self.assertEqual(seen, ["circuit_1.pkl"])


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


class PortableResultTest(BaseTimedTestCase):
    """A listener may write plain data instead of a pickled qililab object.

    The classical cluster has no qililab, so a result pickled as a library
    object cannot even be unpickled there. Both shapes must reach the caller
    through the same ``.counts()`` API.
    """

    def test_plain_dict_is_wrapped_with_the_counts_api(self):
        import tempfile

        from Qdislib.utils.qpu_bridge import QPUResult

        payload = {
            "counts": {"00": 1, "10": 1923, "11": 76},
            "shots": 2000,
            "device": "red",
            "n_qubits": 2,
        }
        with tempfile.TemporaryDirectory() as tmp:
            listener = threading.Thread(target=_fake_listener, args=(tmp, payload))
            listener.start()
            result = submit_and_wait("c", directory=tmp, timeout=5, poll_interval=0.02)
            listener.join()

        self.assertIsInstance(result, QPUResult)
        self.assertEqual(result.counts(), payload["counts"])
        self.assertEqual(result.shots, 2000)
        self.assertEqual(result.device, "red")

    def test_pickled_result_object_still_passes_through_untouched(self):
        import tempfile

        payload = _FakeResult({"01": 500, "10": 500})
        with tempfile.TemporaryDirectory() as tmp:
            listener = threading.Thread(target=_fake_listener, args=(tmp, payload))
            listener.start()
            result = submit_and_wait("c", directory=tmp, timeout=5, poll_interval=0.02)
            listener.join()

        self.assertIsInstance(result, _FakeResult)
        self.assertEqual(result.counts(), {"01": 500, "10": 500})

    def test_shots_default_to_the_sum_of_counts(self):
        from Qdislib.utils.qpu_bridge import QPUResult

        result = QPUResult({"counts": {"0": 3, "1": 7}})
        self.assertEqual(result.shots, 10)

    def test_error_marker_still_wins_over_the_counts_shape(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            payload = {_ERROR_KEY: "boom", "counts": {"0": 1}}
            listener = threading.Thread(target=_fake_listener, args=(tmp, payload))
            listener.start()
            with self.assertRaises(QPUExecutionError):
                submit_and_wait("c", directory=tmp, timeout=5, poll_interval=0.02)
            listener.join()


class _FakePlatform:
    """Stand-in for a qililab platform: records the circuit and returns counts."""

    def execute(self, program, num_bins, num_avg, transpilation_config):
        self.program = program
        return {"counts": {"0" * program.nqubits: num_bins}, "shots": num_bins}


@unittest.skipUnless(importlib.util.find_spec("qibo"), "qibo required")
class ListenerTest(BaseTimedTestCase):
    def _submitted(self, directory):
        from qibo import gates, models

        circuit = models.Circuit(1)
        circuit.add(gates.X(0))
        path = os.path.join(directory, f"{_CIRCUIT_PREFIX}1.pkl")
        _atomic_write_pickle(circuit, path)
        return path

    def test_the_result_goes_next_to_the_circuit(self):
        # the shared directory's own name may contain the file prefix: the writer
        # would wait for a result written somewhere else, until its timeout
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            directory = os.path.join(tmp, "circuit_cutting_run")
            os.makedirs(directory)
            _process_one(_FakePlatform(), self._submitted(directory), 100, None)
            self.assertEqual(os.listdir(directory), [f"{_RESULT_PREFIX}1.pkl"])

    def test_a_circuit_that_cannot_be_pickled_leaves_no_file(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(Exception):
                submit_and_wait(lambda: None, directory=tmp, timeout=1)
            self.assertEqual(os.listdir(tmp), [])

    def test_a_result_that_cannot_be_pickled_is_reported(self):
        import pickle
        import tempfile

        class Unpicklable(_FakePlatform):
            def execute(self, **_):
                return lambda: None

        with tempfile.TemporaryDirectory() as tmp:
            _process_one(Unpicklable(), self._submitted(tmp), 100, None)
            self.assertEqual(os.listdir(tmp), [f"{_RESULT_PREFIX}1.pkl"])
            with open(os.path.join(tmp, f"{_RESULT_PREFIX}1.pkl"), "rb") as handle:
                self.assertIn(_ERROR_KEY, pickle.load(handle))


class _FakeInstruments(_FakePlatform):
    """A qililab platform: records how it was brought up and put down."""

    def __init__(self):
        self.calls = []

    def connect(self):
        self.calls.append("connect")

    def initial_setup(self):
        self.calls.append("setup")

    def turn_on_instruments(self):
        self.calls.append("on")

    def turn_off_instruments(self):
        self.calls.append("off")
        raise RuntimeError("teardown is best effort")

    def disconnect(self):
        self.calls.append("disconnect")


@unittest.skipUnless(importlib.util.find_spec("qibo"), "qibo required")
class ListenerLoopTest(BaseTimedTestCase):
    """The whole hand-off, with a stand-in for qililab."""

    def _qililab(self, platform, runcards):
        def build_platform(runcard):
            runcards.append(runcard)
            return platform

        qililab = types.ModuleType("qililab")
        qililab.logger = logging.getLogger("fake_qililab")
        qililab.build_platform = build_platform
        digital = types.ModuleType("qililab.digital")
        digital.DigitalTranspilationConfig = lambda routing: {"routing": routing}
        return {"qililab": qililab, "qililab.digital": digital}

    def test_a_circuit_goes_there_and_back_and_a_stop_file_ends_it(self):
        import tempfile

        from qibo import gates, models

        platform, runcards = _FakeInstruments(), []
        circuit = models.Circuit(2)
        circuit.add(gates.X(0))
        with tempfile.TemporaryDirectory() as tmp:
            with mock.patch.dict(sys.modules, self._qililab(platform, runcards)):
                listener = threading.Thread(
                    target=run_qpu_listener,
                    kwargs={
                        "directory": tmp,
                        "runcard": "/rc.yml",
                        "nshots": 40,
                        "poll_interval": 0.01,
                    },
                )
                listener.start()
                result = submit_and_wait(
                    circuit, directory=tmp, timeout=30, poll_interval=0.01
                )
                open(os.path.join(tmp, "STOP"), "w").close()
                listener.join(30)
            self.assertFalse(listener.is_alive())
            self.assertEqual(result.counts(), {"00": 40})
            self.assertEqual(os.listdir(tmp), ["STOP"])
        self.assertEqual(runcards, ["/rc.yml"])
        # measured there, on every qubit
        self.assertEqual(type(platform.program.queue[-1]).__name__, "M")
        # the teardown ran, its failure did not get out
        self.assertEqual(platform.calls, ["connect", "setup", "on", "off"])

    def test_idle_timeout_ends_it(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            with mock.patch.dict(sys.modules, self._qililab(_FakeInstruments(), [])):
                run_qpu_listener(
                    directory=tmp, runcard="r", idle_timeout=0.01, poll_interval=0.01
                )


@unittest.skipUnless(importlib.util.find_spec("qibo"), "qibo required")
class BSCBackendShotsTest(BaseTimedTestCase):
    def test_a_listener_with_other_shots_is_pointed_out(self):
        import warnings
        from unittest import mock

        import networkx as nx

        import Qdislib as qd
        from Qdislib.utils import logging_config
        from Qdislib.utils.qpu_bridge import QPUResult

        dag = nx.DiGraph()
        dag.add_node("X_1", gate="x", qubits=(0,), parameters=())
        backend = qd.build_backend("BSC_QPU_BACKEND", qd.BackendConfig(num_qubits=5))
        logging_config._WARNED.clear()
        result = QPUResult({"counts": {"1": 1000}})
        with mock.patch("Qdislib.utils.qpu_bridge.submit_and_wait", lambda _: result):
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                self.assertEqual(backend.evaluate_samples(dag, 1, 1000), {"1": 1000})
                self.assertFalse(caught)
                self.assertEqual(backend.evaluate_samples(dag, 1, 4096), {"1": 1000})
        self.assertTrue(any("nshots" in str(w.message) for w in caught))


if __name__ == "__main__":
    unittest.main()
