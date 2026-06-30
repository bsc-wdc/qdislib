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

"""
Qdislib backend-agnostic helpers.

Utilities to detect which quantum SDK (Qiskit, Qibo or CUDA-Q) a circuit comes
from and to convert it into Qdislib's internal DAG (``networkx``) representation.
Keeping this dispatch in one place avoids duplicating the type-to-converter map
across the cutting algorithms.
"""

import logging

from Qdislib.utils.graph_qibo import circuit_qibo_to_dag
from Qdislib.utils.graph_qiskit import circuit_qiskit_to_dag
from Qdislib.utils.graph_cudaq import circuit_cudaq_to_dag

# Each SDK is optional. When one is absent its circuit class is replaced by a
# placeholder sentinel so the type dispatch below stays valid (no real object is
# ever an instance of the sentinel, so it simply never matches).
try:
    from qiskit.circuit import QuantumCircuit
except ImportError:
    QuantumCircuit = type("QuantumCircuit", (), {})
try:
    from qibo.models import Circuit as QiboCircuit
except ImportError:
    QiboCircuit = type("QiboCircuit", (), {})
try:
    from cudaq.kernel.kernel_decorator import PyKernelDecorator
except ImportError:
    PyKernelDecorator = type("PyKernelDecorator", (), {})

logger = logging.getLogger(__name__)

__all__ = [
    "convert_to_dag",
    "iter_gates",
    "num_qubits",
    "expectation_value",
    "normalize_observables",
    "CutSubcircuits",
]

#: Mapping from a circuit's SDK class to the converter that turns it into a DAG.
_CIRCUIT_TO_DAG = {
    QuantumCircuit: circuit_qiskit_to_dag,
    QiboCircuit: circuit_qibo_to_dag,
    PyKernelDecorator: circuit_cudaq_to_dag,
}


def convert_to_dag(circuit, observables=None):
    """Convert a Qiskit/Qibo/CUDA-Q circuit into a Qdislib DAG.

    Dispatches on the circuit type. If ``circuit`` is already a DAG (or any
    unrecognised type) it is returned unchanged.

    :param circuit: A Qiskit ``QuantumCircuit``, a Qibo ``Circuit``, a CUDA-Q
        kernel, or an already-built ``networkx`` DAG.
    :param observables: Optional observable string; when it contains ``"I"`` the
        identity observable nodes are appended to the resulting DAG.
    :return: The corresponding ``networkx`` DAG, or ``circuit`` unchanged if the
        type is not recognised.
    """
    converter = next(
        (fn for t, fn in _CIRCUIT_TO_DAG.items() if isinstance(circuit, t)), None
    )
    if converter is None:
        logger.warning("Unsupported circuit type: %s", type(circuit).__name__)
        return circuit  # Unknown type (e.g. already a DAG), return unchanged.

    kwargs = {}
    if observables and "I" in observables:
        kwargs["obs_I"] = observables

    return converter(circuit, **kwargs)


def normalize_observables(observables):
    """Normalize an observable specification into ``[(coefficient, pauli_string), ...]``.

    Lets the cutting routines accept a **weighted sum of Pauli operators** (a
    Hamiltonian), not just a single per-qubit Pauli string, since
    ``<sum_i c_i P_i> = sum_i c_i <P_i>``. Accepted forms:

    - a single Pauli string, e.g. ``"ZZI"`` -> ``[(1.0, "ZZI")]``;
    - a mapping ``{pauli_string: coefficient}``, e.g. ``{"ZZ": 0.5, "XI": -1.0}``;
    - an iterable of ``(coefficient, pauli_string)`` pairs, e.g.
      ``[(0.5, "ZZ"), (-1.0, "XI")]``.

    :param observables: A Pauli string, mapping, or iterable of weighted terms.
    :return: A list of ``(float, str)`` ``(coefficient, pauli_string)`` terms.
    :raises ValueError: If a Pauli string contains a character other than I/X/Y/Z,
        or if the terms have inconsistent lengths.
    """
    if isinstance(observables, str):
        terms = [(1.0, observables)]
    elif isinstance(observables, dict):
        terms = [(float(coeff), pauli) for pauli, coeff in observables.items()]
    else:
        terms = [(float(coeff), str(pauli)) for coeff, pauli in observables]

    if not terms:
        raise ValueError("Empty observable specification.")
    lengths = {len(pauli) for _, pauli in terms}
    for _, pauli in terms:
        if any(ch not in "IXYZ" for ch in pauli):
            raise ValueError(
                f"Invalid Pauli string {pauli!r}: only 'I', 'X', 'Y', 'Z' are allowed."
            )
    if len(lengths) > 1:
        raise ValueError(
            f"All Pauli terms must act on the same number of qubits; got lengths {sorted(lengths)}."
        )
    return terms


def change_basis_dag(dag, observables):
    """Append measurement-basis-rotation gates to a circuit DAG (backend-neutral).

    To measure a non-Z Pauli, the qubit is rotated into the Z basis before
    measurement: ``X`` -> append ``H``; ``Y`` -> append ``SDG`` then ``H``. ``Z``
    and ``I`` qubits are left unchanged. The rotation is appended at the qubit's
    current sink (its last gate), so a subsequent Z-basis measurement yields the
    X/Y expectation value.

    This operates on the DAG rather than the source circuit so it works for any
    input backend; the SDK circuits (Qiskit ``QuantumCircuit``, CUDA-Q kernels)
    cannot be appended to generically the way a Qibo ``Circuit`` can.

    :param dag: A circuit DAG produced by :func:`convert_to_dag`.
    :param observables: Per-qubit Pauli string (e.g. ``"XZY"``).
    :return: The same DAG with rotation nodes appended (mutated in place).
    """
    if not observables:
        return dag

    rotations = {"X": ("H",), "Y": ("SDG", "H")}
    for qubit, pauli in enumerate(observables):
        chain = rotations.get(pauli)
        if not chain:
            continue
        for offset, gate_name in enumerate(chain):
            node_id = f"BASIS{pauli}{qubit}_{offset}"
            dag.add_node(node_id, gate=gate_name, qubits=(qubit,), parameters=())
            # Reproduce the dependency edges that ``circuit_*_to_dag`` would build
            # for this gate: a single blue edge to the nearest predecessor on the
            # qubit, and red edges to every earlier gate on that qubit. The red
            # edges matter for wire cutting's component generation.
            for pred in reversed(list(dag.nodes)):
                if pred != node_id and qubit in (dag.nodes[pred].get("qubits") or ()):
                    dag.add_edge(pred, node_id, color="blue")
                    break
            for pred in reversed(list(dag.nodes)):
                if (
                    pred != node_id
                    and qubit in (dag.nodes[pred].get("qubits") or ())
                    and not dag.has_edge(pred, node_id)
                ):
                    dag.add_edge(pred, node_id, color="red")

    return dag


def iter_gates(circuit):
    """Yield each gate of a circuit as a backend-neutral ``(name, qubits)`` pair.

    This normalises the very different gate APIs of the supported SDKs into a
    single ordered list, so callers don't have to special-case the backend:

    - ``name`` is the upper-cased gate name following each SDK's own convention,
      kept consistent with ``circuit_{qibo,qiskit,cudaq}_to_dag`` so the
      ``f"{name}_{index}"`` identifiers line up (e.g. Qibo ``CNOT`` vs Qiskit
      ``CX``).
    - ``qubits`` is a tuple of integer qubit indices the gate acts on.

    One entry is produced per operation, in circuit order, matching the 1-based
    indexing used by the DAG converters.

    :param circuit: A Qibo ``Circuit``, Qiskit ``QuantumCircuit`` or CUDA-Q kernel.
    :raises TypeError: If the circuit type is not recognised.
    :return: List of ``(name, qubits)`` tuples.
    """
    if isinstance(circuit, QiboCircuit):
        return [
            (gate.__class__.__name__.upper(), tuple(gate.qubits))
            for gate in circuit.queue
        ]
    if isinstance(circuit, QuantumCircuit):
        return [
            (
                instruction.operation.name.upper(),
                tuple(q._index for q in instruction.qubits),
            )
            for instruction in circuit.data
        ]
    if isinstance(circuit, PyKernelDecorator):
        # CUDA-Q has no gate-iteration API; reuse the source parser that
        # circuit_cudaq_to_dag already relies on. Node insertion order is gate
        # order, so this stays consistent with that converter's indexing.
        dag = circuit_cudaq_to_dag(circuit)
        return [
            (data["gate"].upper(), tuple(data["qubits"]))
            for _, data in dag.nodes(data=True)
        ]
    raise TypeError(
        f"Unsupported circuit type for gate iteration: {type(circuit).__name__}"
    )


def num_qubits(circuit):
    """Return the number of qubits of a Qiskit/Qibo/CUDA-Q circuit.

    Uses the SDK's native attribute when available (so idle qubits are counted);
    for CUDA-Q kernels, which expose no such attribute, the count is inferred
    from the highest qubit index touched by any gate.

    :param circuit: A Qibo ``Circuit``, Qiskit ``QuantumCircuit`` or CUDA-Q kernel.
    :return: The number of qubits as an ``int``.
    """
    if isinstance(circuit, QiboCircuit):
        return circuit.nqubits
    if isinstance(circuit, QuantumCircuit):
        return circuit.num_qubits
    indices = [qubit for _, qubits in iter_gates(circuit) for qubit in qubits]
    return max(indices) + 1 if indices else 0


def expectation_value(subcircuit, shots=1024):
    """Evaluate one cut subcircuit's expectation value, dispatching on its backend.

    Convenience wrapper for the recommended workflow, so the caller does not have
    to pick the backend-specific evaluator by hand::

        subcircuits = qd.gate_cutting_subcircuits(circuit, cut, software="qibo")
        values = [qd.expectation_value(s) for s in subcircuits]
        result = subcircuits.reconstruct(values)

    The backend is detected from the subcircuit produced by
    :func:`gate_cutting_subcircuits` / :func:`wire_cutting_subcircuits`: a Qibo
    ``Circuit`` is run with the Qibo simulator, a Qiskit ``QuantumCircuit`` with
    Qiskit Aer, and a ``(dag, num_qubits)`` payload with CUDA-Q.

    :param subcircuit: A subcircuit ``[payload, observables]`` (or a bare circuit).
    :param shots: Number of measurement shots.
    :return: The subcircuit's expectation value (a PyCOMPSs future under COMPSs).
    """
    # Evaluators live in the cutting modules; import lazily to avoid an import
    # cycle (those modules import this one).
    from Qdislib.core.cutting_algorithms.wire_cutting import (
        _expec_value_qibo,
        _expec_value_cudaq,
    )
    from Qdislib.core.cutting_algorithms.gate_cutting import _expec_value_qiskit

    if subcircuit is None:  # empty component produced by a cut; contributes nothing
        return None

    payload = subcircuit
    if (
        isinstance(subcircuit, (list, tuple))
        and len(subcircuit) == 2
        and isinstance(subcircuit[1], str)
    ):
        payload = subcircuit[0]

    if isinstance(payload, QiboCircuit):
        return _expec_value_qibo(subcircuit, shots)
    if isinstance(payload, QuantumCircuit):
        return _expec_value_qiskit(subcircuit, shots)
    # CUDA-Q subcircuits carry a ``(dag, num_qubits)`` payload.
    return _expec_value_cudaq(subcircuit, shots)


class CutSubcircuits(list):
    """The subcircuits produced by a cut, that also know how to be recombined.

    This behaves exactly like the plain ``list`` of subcircuits returned before
    (iterate, index, ``len`` as usual), but additionally records the metadata
    needed to recombine the per-subcircuit expectation values and exposes a
    convenient :meth:`reconstruct`::

        subcircuits = qd.wire_cutting_subcircuits(circuit, cut, software="qibo")
        values = [qd.expectation_value(s) for s in subcircuits]
        result = subcircuits.reconstruct(values)        # no manual bookkeeping

    The standalone ``*_subcircuit_reconstruction`` functions remain available for
    callers that prefer to pass ``number_cuts``/``number_components`` explicitly.
    """

    def __init__(self, subcircuits, *, number_cuts, number_components, reconstructor):
        super().__init__(subcircuits)
        #: Number of (single-wire / single-gate) cuts applied.
        self.number_cuts = number_cuts
        #: Number of independent components each cut configuration splits into.
        self.number_components = number_components
        self._reconstructor = reconstructor

    def reconstruct(self, expected_values):
        """Recombine per-subcircuit expectation values into the original one.

        :param expected_values: Expectation values, one per subcircuit, in the
            same order as the subcircuits in this object.
        :return: The reconstructed expectation value of the uncut circuit.
        """
        return self._reconstructor(
            expected_values, self.number_cuts, self.number_components
        )
