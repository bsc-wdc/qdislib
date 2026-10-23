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

"""The calling convention every cutting algorithm shares (since 2.0.0).

All-in-one::

    X(circuit, cuts, observable=None, *, backend=None, shots=1024, seed=None,
      sync=True, <options>)
    X_sampling(circuit, cuts, *, backend=None, shots=1024, seed=None, sync=True,
               <options>)

Three-step::

    X_subcircuits(circuit, cuts, observable=None, *, software=None, <options>)
    X_sampling_subcircuits(circuit, cuts, *, software=None, <options>)

``backend`` is where the subcircuits run, ``software`` the SDK the three-step
subcircuits are written in. The pre-2.0 names keep working for one release with a
:class:`DeprecationWarning` (:func:`standard_arguments`).
"""

import functools
import inspect
import warnings

#: The SDKs a circuit can be written in.
SDKS = ("qibo", "qiskit", "cudaq", "pennylane")

#: Qibo platforms, which ``wire_cutting(backend=...)`` used to take.
_QIBO_PLATFORMS = {
    "numpy",
    "qibojit",
    "numba",
    "cupy",
    "cuquantum",
    "tensorflow",
    "pytorch",
    "jax",
}

#: Renamed arguments: old name -> standard name.
_RENAMED = {
    "observables": "observable",
    "gates_cut": "cuts",
    "cut": "cuts",
    "dag": "circuit",
    "rand_qc": "circuit",
}

#: Execution options that ``backend=`` replaces; still honoured for one release.
_LEGACY_EXECUTION = ("qpu", "qpu_dict", "qpu_config", "gpu", "gpu_min_qubits", "method")


def _deprecated(func, message):
    warnings.warn(
        f"{func.__name__}(): {message}; the old name is removed in the next release.",
        DeprecationWarning,
        stacklevel=3,
    )


def standard_arguments(func):
    """Accept the pre-2.0 argument names of a cutting function, with a warning.

    - ``observables=``, ``gates_cut=``/``cut=``, ``dag=``/``rand_qc=`` are renamed;
    - ``software=`` becomes ``backend=`` on functions that run circuits (the
      three-step generators keep ``software``: the SDK of their subcircuits);
    - a Qibo platform passed as ``backend="numpy"`` becomes ``backend="qibo"``;
    - ``qpu``, ``qpu_dict``, ``qpu_config``, ``gpu``, ``gpu_min_qubits``, ``method``
      are handed to the function as ``_legacy`` (when it still honours them);
    - the ``cuts`` may name a gate as any SDK does (``CX_7``, ``CNOT_7``, ``cx_7``):
      :func:`~Qdislib.utils.gates.resolve_cuts` maps them to the circuit's names.

    Positional calls beyond ``(circuit, cuts, observable)`` raise a ``TypeError``
    that says every other argument is keyword-only.
    """
    signature = inspect.signature(func)
    params = signature.parameters
    positional = [p.name for p in params.values() if p.kind is p.POSITIONAL_OR_KEYWORD]

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        if len(args) > len(positional):
            raise TypeError(
                f"{func.__name__}() takes ({', '.join(positional)}) by position; "
                "since 2.0.0 every other argument is keyword-only, e.g. "
                "shots=..., backend=..., sync=..."
            )
        if (
            "observable" in positional
            and len(args) > positional.index("observable")
            and isinstance(args[positional.index("observable")], bool)
        ):
            raise TypeError(
                f"{func.__name__}(): the third positional argument is the observable "
                "since 2.0.0; pass sync= by keyword."
            )
        for old, new in _RENAMED.items():
            if old in kwargs:
                if new in kwargs:
                    raise TypeError(f"{func.__name__}() got both {old}= and {new}=")
                _deprecated(func, f"{old}= is now {new}=")
                kwargs[new] = kwargs.pop(old)
        legacy = {}
        if "software" in kwargs and "software" not in params:
            if "backend" in kwargs:
                raise TypeError(f"{func.__name__}() got both software= and backend=")
            _deprecated(func, "software= is now backend=")
            kwargs["backend"] = kwargs.pop("software")
        backend = kwargs.get("backend")
        if isinstance(backend, str) and backend in _QIBO_PLATFORMS:
            _deprecated(
                func, f'backend="{backend}" (a Qibo platform) is now backend="qibo"'
            )
            legacy["platform"] = backend
            kwargs["backend"] = "qibo"
        used = [name for name in _LEGACY_EXECUTION if name in kwargs]
        if used:
            _deprecated(func, f"{'=, '.join(used)}= are replaced by backend=")
            for name in used:
                legacy[name] = kwargs.pop(name)
        if legacy:
            if "_legacy" not in params:
                raise TypeError(
                    f"{func.__name__}() does not take {sorted(legacy)}; use backend="
                )
            kwargs["_legacy"] = legacy
        if "cuts" in positional and "circuit" in positional:
            # A cut may name its gate as any SDK does (CX_7 / CNOT_7 / cx_7).
            from Qdislib.utils.gates import resolve_cuts

            bound = signature.bind_partial(*args, **kwargs)
            if "cuts" in bound.arguments and "circuit" in bound.arguments:
                bound.arguments["cuts"] = resolve_cuts(
                    bound.arguments["circuit"], bound.arguments["cuts"]
                )
                args, kwargs = bound.args, bound.kwargs
        return func(*args, **kwargs)

    wrapper.__signature__ = signature.replace(
        parameters=[p for p in params.values() if p.name != "_legacy"]
    )
    return wrapper


def resolve_backend(backend, circuit):
    """What the executors take as ``software``: an SDK name, a registered backend
    name or a backend object. ``None`` is the circuit's own SDK (Qiskit Aer for a
    DAG, the registry's default)."""
    if backend is None:
        from Qdislib.utils.circuit_cache import _detect_software

        return _detect_software(circuit) or "qiskit"
    if isinstance(backend, str) and backend not in SDKS:
        from Qdislib.core.backend.registry import available_backends

        if backend not in available_backends():
            from Qdislib.utils.exceptions import UnsupportedBackendError

            raise UnsupportedBackendError(
                f"Unknown backend {backend!r}: use one of {list(SDKS)}, a registered "
                f"backend ({available_backends()}) or a build_backend(...) object."
            )
    return backend


def resolve_software(software, circuit):
    """The SDK the three-step subcircuits are written in (``None``: the circuit's)."""
    if software is None:
        from Qdislib.utils.circuit_cache import _detect_software

        return _detect_software(circuit) or "qiskit"
    if software not in SDKS:
        from Qdislib.utils.exceptions import UnsupportedBackendError

        raise UnsupportedBackendError(
            f"software={software!r}: the subcircuits can be written in {list(SDKS)}. "
            "To run them somewhere else, pass backend= to expectation_value / "
            "sample_subcircuit."
        )
    return software


def default_observable(observable, circuit):
    """``None`` means ``Z`` on every qubit."""
    if observable is not None:
        return observable
    from Qdislib.core.cutting_algorithms.hadamard_cut import _num_circuit_qubits

    return "Z" * _num_circuit_qubits(circuit)


def is_simulator(backend):
    """True for an SDK name (its local simulator)."""
    return isinstance(backend, str) and backend in SDKS
