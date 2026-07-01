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

"""Cross-cluster QPU bridge for Qdislib.

At BSC the classical (PyCOMPSs) job and the quantum-control job run on separate
clusters that cannot be launched from a single scheduler allocation, so they
hand subcircuits off through a **shared filesystem directory**:

- the classical cut workflow writes ``circuit_<id>.pkl`` and waits (bounded) for
  ``result_<id>.pkl`` -- :func:`submit_and_wait`, used by the QPU evaluators;
- on the quantum cluster a listener executes each circuit on a QPU (BSC **Blue**
  or **Red**) and writes the result back -- :func:`run_qpu_listener`.

Both sides resolve the directory through :func:`subcircuits_dir`
(``$QDISLIB_SUBCIRCUITS_DIR``), so writer and listener always meet in the same
place. Files are written atomically (temp + ``os.replace``) so neither side ever
reads a half-written file, the wait is bounded (:class:`QPUTimeoutError`), and a
listener-side failure is reported back (:class:`QPUExecutionError`) so the writer
fails fast instead of hanging.

The listener needs ``qililab`` (the BSC quantum-control library) and a platform
runcard; the writer side needs neither, so importing Qdislib on the classical
cluster does not require ``qililab``.
"""

import glob
import logging
import os
import pickle
import time

from Qdislib.utils.exceptions import QdislibException

logger = logging.getLogger(__name__)

#: Default shared directory (the original BSC path); override with the
#: ``QDISLIB_SUBCIRCUITS_DIR`` environment variable on both clusters.
DEFAULT_SUBCIRCUITS_DIR = "/home/bsc/bsc019635/ona_proves/subcircuits"

_CIRCUIT_PREFIX = "circuit_"
_RESULT_PREFIX = "result_"
_STOP_FILE = "STOP"
_ERROR_KEY = "__qdislib_qpu_error__"

#: BSC has two QPUs; each has its own qililab runcard behind an env var.
_QPU_RUNCARD_ENV = {"red": "RUNCARD_RED", "blue": "RUNCARD_BLUE"}

__all__ = [
    "subcircuits_dir",
    "run_qpu_listener",
    "QPUTimeoutError",
    "QPUExecutionError",
]


class QPUTimeoutError(QdislibException, TimeoutError):
    """Raised when a QPU result does not arrive within the timeout."""


class QPUExecutionError(QdislibException, RuntimeError):
    """Raised when the listener reports it failed to execute a submitted circuit."""


def subcircuits_dir():
    """Return the shared hand-off directory (``$QDISLIB_SUBCIRCUITS_DIR``).

    Single source of truth shared by the writer (:func:`submit_and_wait`) and the
    listener (:func:`run_qpu_listener`), so they always agree on where to meet.
    """
    return os.environ.get("QDISLIB_SUBCIRCUITS_DIR", DEFAULT_SUBCIRCUITS_DIR)


def _atomic_write_pickle(obj, path):
    """Pickle ``obj`` to ``path`` atomically (temp file + ``os.replace``).

    So a reader on the other cluster never observes a partially written file.
    """
    tmp = f"{path}.{os.getpid()}.tmp"
    with open(tmp, "wb") as handle:
        pickle.dump(obj, handle)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)  # atomic within a single filesystem


def submit_and_wait(circuit, directory=None, timeout=None, poll_interval=1.0):
    """Ship ``circuit`` to the shared dir and wait (bounded) for the result.

    Writes ``circuit_<id>.pkl`` and polls for the listener's ``result_<id>.pkl``.
    The wait is bounded so a down listener surfaces a clear error instead of
    hanging the task forever.

    :param circuit: The (unmeasured) subcircuit to execute; the listener adds the
        measurements and runs it on the QPU.
    :param directory: Shared hand-off dir (default :func:`subcircuits_dir`).
    :param timeout: Max seconds to wait; defaults to ``$QDISLIB_QPU_TIMEOUT`` or
        3600. ``None`` waits forever (not recommended).
    :param poll_interval: Seconds between checks for the result file.
    :raises QPUTimeoutError: If no result appears within ``timeout``.
    :raises QPUExecutionError: If the listener reports an execution failure.
    :return: The result object written by the listener (exposes ``.counts()``).
    """
    directory = directory or subcircuits_dir()
    if timeout is None:
        timeout = float(os.environ.get("QDISLIB_QPU_TIMEOUT", 3600.0))
    os.makedirs(directory, exist_ok=True)

    unique = f"{time.time_ns()}_{os.getpid()}"
    circuit_path = os.path.join(directory, f"{_CIRCUIT_PREFIX}{unique}.pkl")
    result_path = os.path.join(directory, f"{_RESULT_PREFIX}{unique}.pkl")

    _atomic_write_pickle(circuit, circuit_path)
    logger.info(
        "QPU bridge: submitted %s (timeout %ss)",
        os.path.basename(circuit_path),
        timeout,
    )

    deadline = None if timeout is None else time.monotonic() + timeout
    while not os.path.exists(result_path):
        if deadline is not None and time.monotonic() > deadline:
            # Clean up so a late listener does not leave orphaned files behind.
            for path in (circuit_path, result_path):
                try:
                    os.remove(path)
                except OSError:
                    pass
            raise QPUTimeoutError(
                f"No QPU result after {timeout}s in {directory!r}. "
                "Is run_qpu_listener() running against this directory?"
            )
        time.sleep(poll_interval)

    with open(result_path, "rb") as handle:
        result = pickle.load(handle)
    os.remove(result_path)

    if isinstance(result, dict) and _ERROR_KEY in result:
        raise QPUExecutionError(result[_ERROR_KEY])
    return result


def _resolve_runcard(runcard, qpu):
    """Resolve the qililab runcard path from an explicit path or a QPU name."""
    if runcard is not None:
        return runcard
    if qpu is not None:
        env = _QPU_RUNCARD_ENV.get(qpu.lower())
        if env is None:
            raise ValueError(f"Unknown QPU {qpu!r}; expected 'red' or 'blue'.")
        value = os.environ.get(env)
        if not value:
            raise ValueError(
                f"QPU {qpu!r} selected but ${env} is not set to a runcard path."
            )
        return value
    for env in ("RUNCARD_RED", "RUNCARD_BLUE", "RUNCARD"):
        value = os.environ.get(env)
        if value:
            return value
    raise ValueError(
        "No qililab runcard: pass runcard=, or qpu='red'/'blue' with "
        "$RUNCARD_RED / $RUNCARD_BLUE set."
    )


def _process_one(platform, circuit_path, nshots, transpilation):
    """Execute one submitted circuit; write its result (or an error marker).

    A failure on one circuit is logged and reported back to the writer via an
    error marker -- it never crashes the listener loop.
    """
    from qibo import gates

    result_path = circuit_path.replace(_CIRCUIT_PREFIX, _RESULT_PREFIX)
    try:
        with open(circuit_path, "rb") as handle:
            circuit = pickle.load(handle)

        # The writer ships the unmeasured subcircuit; add the measurements here,
        # on the qubits the QPU will read out, right before the hardware run.
        circuit.wire_names = list(range(circuit.nqubits))
        circuit.add(gates.M(*range(circuit.nqubits)))

        result = platform.execute(
            program=circuit,
            num_bins=nshots,
            num_avg=1,
            transpilation_config=transpilation,
        )
        _atomic_write_pickle(result, result_path)
        logger.info("QPU listener: executed %s", os.path.basename(circuit_path))
    except Exception as exc:  # never let one bad circuit kill the loop
        logger.exception("QPU listener: failed on %s", os.path.basename(circuit_path))
        try:
            _atomic_write_pickle(
                {_ERROR_KEY: f"{type(exc).__name__}: {exc}"}, result_path
            )
        except Exception:
            logger.debug("QPU listener: could not write error marker", exc_info=True)
    finally:
        try:
            os.remove(circuit_path)
        except OSError:
            pass


def run_qpu_listener(
    directory=None,
    qpu="red",
    runcard=None,
    nshots=1000,
    poll_interval=1.0,
    idle_timeout=None,
    routing=True,
):
    """Run the QPU-side listener loop on the quantum cluster.

    Watches the shared directory for ``circuit_*.pkl`` submitted by a Qdislib cut
    workflow, executes each on a BSC QPU via ``qililab``, and writes back
    ``result_*.pkl`` (atomically). Runs until a ``STOP`` file appears in the
    directory, a ``KeyboardInterrupt``, or ``idle_timeout`` idle seconds elapse.

    Launch it as its own job on the quantum cluster (see
    ``examples/qpu_listener.py``); the classical cut job submits work to the same
    ``$QDISLIB_SUBCIRCUITS_DIR`` from the other cluster.

    :param directory: Shared hand-off dir (default :func:`subcircuits_dir`).
    :param qpu: ``"red"`` or ``"blue"`` -- selects ``$RUNCARD_RED`` /
        ``$RUNCARD_BLUE``. Ignored if ``runcard`` is given.
    :param runcard: Explicit qililab runcard path (overrides ``qpu``).
    :param nshots: Shots (qililab ``num_bins``) per circuit execution.
    :param poll_interval: Seconds between directory scans when idle.
    :param idle_timeout: Stop after this many idle seconds (``None`` = run until
        a STOP file or KeyboardInterrupt).
    :param routing: Enable qililab digital-transpilation routing.
    :raises ImportError: If ``qililab`` is not installed.
    :raises ValueError: If no runcard can be resolved.
    """
    import qililab as ql
    from qililab.digital import DigitalTranspilationConfig

    directory = directory or subcircuits_dir()
    os.makedirs(directory, exist_ok=True)
    runcard = _resolve_runcard(runcard, qpu)

    ql.logger.setLevel(logging.ERROR)  # keep qililab quiet unless it errors
    logger.info("QPU listener: connecting to the %s QPU (%s)", qpu, runcard)
    platform = ql.build_platform(runcard=runcard)
    platform.connect()
    platform.initial_setup()
    platform.turn_on_instruments()
    transpilation = DigitalTranspilationConfig(routing=routing)
    logger.info("QPU listener: watching %s", directory)

    stop_path = os.path.join(directory, _STOP_FILE)
    last_activity = time.monotonic()
    try:
        while True:
            if os.path.exists(stop_path):
                logger.info("QPU listener: STOP file found; shutting down.")
                break
            pending = sorted(
                glob.glob(os.path.join(directory, f"{_CIRCUIT_PREFIX}*.pkl"))
            )
            if not pending:
                idle = time.monotonic() - last_activity
                if idle_timeout is not None and idle > idle_timeout:
                    logger.info("QPU listener: idle %.0fs; shutting down.", idle)
                    break
                time.sleep(poll_interval)
                continue
            for circuit_path in pending:
                _process_one(platform, circuit_path, nshots, transpilation)
            last_activity = time.monotonic()
    except KeyboardInterrupt:
        logger.info("QPU listener: interrupted; shutting down.")
    finally:
        try:
            platform.turn_off_instruments()
            platform.disconnect()
        except Exception:  # best-effort teardown
            logger.debug("QPU listener: platform teardown raised", exc_info=True)
    logger.info("QPU listener: stopped.")
