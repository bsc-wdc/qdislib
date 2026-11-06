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
Qdislib circuit utils.

This file contains all auxiliary circuit classes and functions.
"""

from __future__ import annotations  # defer annotation eval so qibo can be optional

import logging

logger = logging.getLogger(__name__)

__all__ = [
    "analytical_solution",
    "exact_statevector",
    "draw_to_circuit",
    "parse_qsim",
    "from_qasm",
]

import inspect
import numpy as np
import random
import typing
from collections import defaultdict

try:
    from qibo import models, gates

    qibo_available = True
except ImportError:  # qibo is optional; functions that need it fail when called
    models = gates = None
    qibo_available = False

from Qdislib.utils.exceptions import QdislibException
from Qdislib.utils.backends import normalize_observables


def from_qasm(qasm, software="qibo"):
    """Build a circuit from an OpenQASM 2.0 string.

    A lightweight interoperability entry point: any framework that exports
    OpenQASM 2.0 (Cirq, PennyLane, Amazon Braket, Qiskit, ...) can feed circuits
    into Qdislib's cutting workflows through it, without a dedicated backend.

    :param qasm: An OpenQASM 2.0 program (string).
    :param software: The SDK of the circuit to return (``"qibo"`` by default; any
        name of :func:`~Qdislib.utils.sdks.sdk_names`, e.g. ``"qiskit"``,
        ``"pennylane"``, ``"qilisdk"``); pick whichever you will run on.
    :raises ImportError: If the requested backend SDK is not installed.
    :raises ValueError: For an unknown ``software``.
    :return: A circuit of ``software``: a Qibo ``Circuit``, a Qiskit
        ``QuantumCircuit``, a PennyLane ``QuantumTape``, ...

    Example:
        .. code-block:: python

            import cirq
            import Qdislib as qd

            qc = qd.from_qasm(cirq.qasm(cirq_circuit), software="qiskit")
            cut = qd.find_cut(qc, gate_cut=True, wire_cut=False)
    """
    if software == "qibo":
        if not qibo_available:
            raise ImportError(
                "from_qasm(software='qibo') requires qibo (install Qdislib[qibo])."
            )
        return models.Circuit.from_qasm(qasm)
    if software == "qiskit":
        try:
            from qiskit import QuantumCircuit
        except ImportError as exc:
            raise ImportError(
                "from_qasm(software='qiskit') requires qiskit "
                "(install Qdislib[qiskit])."
            ) from exc
        return QuantumCircuit.from_qasm_str(qasm)
    if software == "pennylane":
        try:
            import pennylane as qml
        except ImportError as exc:
            raise ImportError(
                "from_qasm(software='pennylane') requires pennylane "
                "(install Qdislib[pennylane])."
            ) from exc
        # PennyLane's own QASM loader needs the extra pennylane-qiskit plugin, so
        # instead parse the QASM with whichever standard parser is present (Qiskit,
        # else Qibo) and rebuild it as a PennyLane QuantumTape through Qdislib's DAG.
        from Qdislib.utils.backends import convert_to_dag, num_qubits as _num_qubits
        from Qdislib.utils.graph_pennylane import (
            _apply_pennylane_gate,
            dag_to_circuit_pennylane,
        )

        try:
            base = from_qasm(qasm, software="qiskit")
        except ImportError:
            base = from_qasm(qasm, software="qibo")
        ops = dag_to_circuit_pennylane(convert_to_dag(base), _num_qubits(base))[0].ops
        with qml.tape.QuantumTape() as tape:
            for gate, qubits, params in ops:
                _apply_pennylane_gate(gate, qubits, params)
        return tape
    from Qdislib.utils.sdks import sdk, sdk_names

    entry = sdk(software)
    if entry is None:
        raise ValueError(
            f"Unknown software {software!r}; expected one of {list(sdk_names())}."
        )
    if not entry.installed():
        raise ImportError(
            f"from_qasm(software={software!r}) requires {entry.package} "
            f"(install Qdislib[{software}])."
        )
    # Parse with a standard parser, then build the circuit through the DAG.
    from Qdislib.utils.backends import convert_to_dag, num_qubits as _num_qubits

    try:
        base = from_qasm(qasm, software="qiskit")
    except ImportError:
        base = from_qasm(qasm, software="qibo")
    built = entry.call("from_dag", convert_to_dag(base), _num_qubits(base))
    return built[0] if isinstance(built, (list, tuple)) else built


# Single-qubit Pauli matrices for the exact expectation-value oracle.
_PAULI_MATRIX = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
    "Y": np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex),
    "Z": np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex),
}


def exact_statevector(circuit: typing.Any) -> np.ndarray:
    """Exact statevector of ``circuit`` as a flat numpy array, qubit-0-first.

    Backend-agnostic oracle building block: accepts a circuit of any SDK (a Qibo
    ``Circuit``, a Qiskit ``QuantumCircuit``, a CUDA-Q kernel, a PennyLane
    ``QuantumTape`` / ``QNode``, a QiliSDK ``Circuit``), and returns the state in a single convention -- amplitude-index bit ``n-1-q`` is
    qubit ``q`` (qubit 0 is the most-significant bit) -- so it lines up with an
    ``observables`` string indexed qubit-0-first regardless of each SDK's native
    endianness (Qibo and PennyLane are already qubit-0-first; Qiskit and CUDA-Q are
    qubit-0-last and get reversed here).

    :param circuit: A circuit, without mid-circuit measurements.
    :raises TypeError: If ``circuit`` is not a supported circuit type.
    :return: The statevector as a length-``2**n`` complex numpy array.

    Example::

        qc = QuantumCircuit(2)
        qc.x(0)
        qd.exact_statevector(qc)  # amplitude 1 at index 0b10: qubit 0 is the first bit
    """
    from Qdislib.utils.sdks import sdk_names, sdk_of

    entry = sdk_of(circuit)
    if entry is None or not entry.statevector:
        raise TypeError(
            f"Unsupported circuit type {type(circuit)!r}; "
            f"expected a circuit of {list(sdk_names())}."
        )
    return np.asarray(entry.call("statevector", circuit), dtype=complex).reshape(-1)


def _pauli_expectation(state: np.ndarray, pauli: str) -> float:
    """Exact ``<psi|P|psi>`` for a Pauli string ``P`` on a qubit-0-first ``state``.

    Applies each single-qubit Pauli to the statevector tensor (``O(2**n)``) instead
    of forming the dense ``2**n x 2**n`` operator (``O(4**n)``), so the oracle scales
    to many more qubits than a full-matrix approach.
    """
    n = len(pauli)
    psi = np.asarray(state, dtype=complex).reshape([2] * n)
    out = psi
    for qubit, ch in enumerate(pauli):
        if ch == "I":
            continue
        out = np.moveaxis(
            np.tensordot(_PAULI_MATRIX[ch], out, axes=([1], [qubit])), 0, qubit
        )
    return float(np.real(np.vdot(psi, out)))


def analytical_solution(circuit: typing.Any, observables: str) -> float:
    """
    Compute the analytical expectation value of a given quantum circuit.

    This function evaluates the exact (statevector) expectation value of a quantum
    circuit with respect to a set of observables. It is backend-agnostic: the circuit
    may be of any SDK -- each one's native endianness is handled so ``observables`` is
    always indexed qubit-0-first (see :func:`exact_statevector`).

    :param circuit:
        A circuit of any SDK, without mid-circuit measurements.
    :param observables:
        A Pauli string (e.g. ``"XZI"``), one letter per qubit, qubit 0 first; or a
        weighted sum of them (anything
        :func:`~Qdislib.utils.backends.normalize_observables` takes).
    :raises TypeError: If ``circuit`` is not a supported circuit type.
    :raises ValueError: For a letter other than I, X, Y, Z, or a string without one
        letter per qubit.
    :return:
        The computed expectation value as a floating-point number.

    Per-term expectation values are logged at DEBUG level; enable them with
    ``qd.set_log_level("DEBUG")``.

    Example:

    .. code-block:: python

        from qibo import models, gates
        import Qdislib as qd

        circuit = models.Circuit(3)
        circuit.add(gates.H(0))
        circuit.add(gates.CNOT(0, 1))

        qd.analytical_solution(circuit, "ZZI")  # 1.0: the two qubits always agree
    """

    # Backend-agnostic exact reference: get the statevector in a single qubit-0-first
    # convention (see exact_statevector) and evaluate <psi|P|psi> per Pauli term with
    # a scalable per-qubit contraction. ``observables`` may be a single Pauli string
    # or a weighted sum of them; normalize_observables handles both.
    psi = exact_statevector(circuit)
    total = 0.0
    for coeff, pauli in normalize_observables(observables):
        term_value = _pauli_expectation(psi, pauli)
        logger.debug(
            "Term %s (coefficient %s) expectation value: %s",
            pauli,
            coeff,
            term_value,
        )
        total += coeff * term_value

    return total


def _random_circuit(
    qubits: int, gate_max, num_cz: typing.Optional[int], p: typing.Any
) -> models.Circuit:
    """Generate a random circuit.

    :param qubits: Number of qubits.
    :param gate_max: Maximum number of single-qubit gates between CZ.
    :param num_cz: Number of CZ qbit gates. If none, p will be taken.
    :param p: Probability of adding an edge.
    :raise QdislibException: Unless exactly one of num_cz and p is given.
    :return: New random circuit.
    """
    import networkx

    # Seeded from `random`, so random.seed() keeps the whole circuit reproducible.
    seed = random.randrange(2**32)
    if (p is None) == (num_cz is None):
        raise QdislibException(
            "Error: exactly one of the number of edges and the probability of adding "
            "an edge must be specified"
        )
    if p is None:
        graph = networkx.gnm_random_graph(qubits, num_cz, seed=seed)
    else:
        graph = networkx.gnp_random_graph(qubits, p, seed=seed)

    edge_list = list(graph.edges())

    gates_pull = [gates.X, gates.H, gates.S, gates.T]  # pull of single-qubit gates
    circuit = models.Circuit(qubits)
    for edge in edge_list:
        # Number of single-qubit gates between CZ
        rand_tmp = random.randint(0, gate_max)
        for _ in range(rand_tmp):
            # Gate selected from the pull
            sel_gate = random.choice(gates_pull)
            # Qubit selected to apply the gate
            sel_qubit = random.randint(0, qubits - 1)
            circuit.add(sel_gate(sel_qubit))
        # 2-qubit gate from graph
        circuit.add(gates.CZ(edge[0], edge[1]))

    return circuit


def draw_to_circuit(
    text_draw: str, parameters: typing.Optional[typing.List[typing.Any]] = None
) -> models.Circuit:
    """Convert text circuit to circuit object.

    This function takes a string representation of a quantum circuit and converts it
    into a Qibo `models.Circuit` object. The input string should represent the circuit in a
    "text drawing" format, where qubits are represented by lines and gates are represented
    by various symbols. The function processes these text representations and constructs
    the equivalent quantum circuit object.

    :param text_draw:
        A string representing the quantum circuit in a text drawing format, where
        qubits are represented by lines and gates are represented by symbols.
    :param parameters:
        An optional list of parameters to be passed to the gates, defaults to None.
    :return:
        A Qibo `models.Circuit` object corresponding to the input text drawing.

    :example:

    .. code-block:: python

        from qibo.models import Circuit
        import Qdislib as qd

        # Define the text representation of the quantum circuit
        text_rep = '''
        q0: ─H─S─X───SDG───H───
        q1: ─X─S───o─X─T─o─S───
        q2: ───────Z─T───|───o─
        q3: ─────────────|─X─|─
        q4: ──RX─────RX──Z───Z─
        '''

        # Convert the text drawing into a quantum circuit
        circuit = qd.draw_to_circuit(text_rep)

        # Print the resulting circuit
        print(circuit.draw())
    """
    split = text_draw.splitlines()
    qubits_lst = []
    split = [element for element in split if element.strip()]
    split = [element for element in split if element != ""]
    for line in split:
        index = line.index("─")
        qubits_lst.append(line[index:])

    list_multiple_gates = defaultdict(list)
    # Now we will process each line to identify multi-qubit gates
    for idx, qubit_line in enumerate(qubits_lst):
        qubit_state = list(qubit_line)

        # Boolean to track if we are inside a multi-qubit gate
        for i, symbol in enumerate(qubit_state):
            if symbol == "o":
                index = i
                for idx2, qubit in enumerate(qubits_lst[idx + 1 :]):
                    if list(qubit)[index] != "|":
                        name = list(qubit)[index]
                        if name == "Z":
                            name = "CZ"
                        elif name == "X":
                            name = "CNOT"
                        list_multiple_gates[idx].append((name, (idx, idx2 + idx + 1)))
                        qubits_lst[idx2 + idx + 1] = (
                            qubits_lst[idx2 + idx + 1][:index]
                            + "─"
                            + qubits_lst[idx2 + idx + 1][index + 1 :]
                        )
                        break

    circuit = models.Circuit(len(qubits_lst))
    num_steps = len(list(qubits_lst[0]))  # Total number of time steps (columns)

    for step in range(num_steps):
        saved_qubit = []
        for idx, qubit_line in enumerate(qubits_lst):
            qubit_state = list(qubit_line)
            parameter_tracker = 0

            char = qubit_state[step]
            # Treat out-of-range neighbours (first/last column) as empty wire so
            # a gate at the circuit boundary does not raise IndexError.
            prev_char = qubit_state[step - 1] if step > 0 else "─"
            next_char = qubit_state[step + 1] if step + 1 < len(qubit_state) else "─"
            if char != "─" and char != "|":
                if char != "o":
                    if next_char == "─" and prev_char == "─":
                        tmp = char
                        gate_name = tmp
                        qubits = idx

                        # Get the gate class from the qibo.gates module
                        gate_class = getattr(gates, gate_name)

                        # Get the signature of the gate's __init__ method
                        signature = inspect.signature(gate_class.__init__)

                        # Count the number of required positional arguments (excluding 'self')
                        param_count = len(signature.parameters) - 1  # exclude 'self'

                        # Check if parameters are provided and the gate requires them
                        if parameters is not None and param_count > 1:
                            param = parameters[idx][parameter_tracker][1]
                            # Pass qubits and parameters if the gate requires both
                            circuit.add(gate_class(qubits, param))
                            parameter_tracker += 1
                        else:
                            # Otherwise, pass only the qubits
                            circuit.add(gate_class(qubits))

                    elif prev_char == "─" and next_char != "─":
                        tmp = ""
                        for i in range(step, num_steps):
                            if (
                                qubit_state[i + 1] if i + 1 < len(qubit_state) else "─"
                            ) == "─":
                                tmp = tmp + qubit_state[i]

                                gate_name = tmp
                                qubits = idx
                                # Get the gate class from the qibo.gates module
                                gate_class = getattr(gates, gate_name)

                                # Get the signature of the gate's __init__ method
                                signature = inspect.signature(gate_class.__init__)

                                # Count the number of required positional arguments (excluding 'self')
                                param_count = (
                                    len(signature.parameters) - 1
                                )  # exclude 'self'

                                # Check if parameters are provided and the gate requires them
                                if parameters is not None and param_count > 1:
                                    param = parameters[idx][parameter_tracker][1]
                                    # Pass qubits and parameters if the gate requires both
                                    circuit.add(gate_class(qubits, param))
                                    parameter_tracker += 1
                                    break
                                else:
                                    circuit.add(gate_class(qubits))
                                    break

                            else:
                                tmp = tmp + qubit_state[i]

                elif char == "o":
                    saved_qubit.append(idx)

        for idx in saved_qubit:
            circuit.add(
                getattr(gates, list_multiple_gates[idx][0][0])(
                    *list_multiple_gates[idx][0][1]
                )
            )
            list_multiple_gates[idx].remove(list_multiple_gates[idx][0])
    return circuit


#: qsim's gates (its ``lib/gates_qsim.h``) as qibo gates, up to a global phase:
#: ``name -> (qubits, angles, the gates for those qubits and angles)``.
_QSIM_GATES = {
    "id1": (1, 0, lambda q, _: [gates.I(q[0])]),
    "h": (1, 0, lambda q, _: [gates.H(q[0])]),
    "t": (1, 0, lambda q, _: [gates.T(q[0])]),
    "s": (1, 0, lambda q, _: [gates.S(q[0])]),
    "x": (1, 0, lambda q, _: [gates.X(q[0])]),
    "y": (1, 0, lambda q, _: [gates.Y(q[0])]),
    "z": (1, 0, lambda q, _: [gates.Z(q[0])]),
    "x_1_2": (1, 0, lambda q, _: [gates.RX(q[0], np.pi / 2)]),
    "y_1_2": (1, 0, lambda q, _: [gates.RY(q[0], np.pi / 2)]),
    # the square root of (X + Y) / sqrt(2): a pi/2 turn about that axis
    "hz_1_2": (1, 0, lambda q, _: [gates.U3(q[0], np.pi / 2, -np.pi / 4, np.pi / 4)]),
    "rx": (1, 1, lambda q, a: [gates.RX(q[0], a[0])]),
    "ry": (1, 1, lambda q, a: [gates.RY(q[0], a[0])]),
    "rz": (1, 1, lambda q, a: [gates.RZ(q[0], a[0])]),
    # a turn of a[1] about the axis at angle a[0] in the XY plane
    "rxy": (
        1,
        2,
        lambda q, a: [gates.U3(q[0], a[1], a[0] - np.pi / 2, np.pi / 2 - a[0])],
    ),
    "id2": (2, 0, lambda q, _: [gates.I(q[0]), gates.I(q[1])]),
    "cz": (2, 0, lambda q, _: [gates.CZ(q[0], q[1])]),
    "cnot": (2, 0, lambda q, _: [gates.CNOT(q[0], q[1])]),
    "sw": (2, 0, lambda q, _: [gates.SWAP(q[0], q[1])]),
    "is": (2, 0, lambda q, _: [gates.iSWAP(q[0], q[1])]),
    "fs": (2, 2, lambda q, a: [gates.fSim(q[0], q[1], a[0], a[1])]),
    "cp": (2, 1, lambda q, a: [gates.CU1(q[0], q[1], -a[0])]),  # qsim's is exp(-i phi)
}


def parse_qsim(fname: str, depth: int = 22):
    """A Qibo :class:`Circuit` from a qsim input file, the same unitary up to a
    global phase. The angles are qsim's, in radians.

    The file starts with the number of qubits; each other line is ``<time> <gate>
    <qubits> [angles]``. Gates at time ``depth`` or later are left out. Its ``fs``
    gates stay Qibo's ``fSim``, which a cut runs on Qibo and Qiskit.

    :param fname: Path to the .qsim file.
    :param depth: Gates from this time on are left out. Default is 22.
    :raises ValueError: If the first line is not the number of qubits, a line does not
        have its gate's qubits and angles, or a gate is not one of ``id1 h t s x y z
        x_1_2 y_1_2 hz_1_2 rx ry rz rxy id2 cz cnot sw is fs cp`` (measurements,
        ``mat1``, ``mat2`` and controlled gates are refused, not dropped).
    :return: A :class:`models.Circuit` from Qibo.

    :example:

    Suppose you have a `.qsim` file with the following content:

    .. code-block:: text

        3
        0 rz 0 0.5
        1 x_1_2 1
        2 rz 2 1.0
        3 hz_1_2 0
        4 fs 1 2 0.5 0.1

    You can parse this file and create a circuit like this:

    .. code-block:: python

        import Qdislib as qd

        circuit = qd.parse_qsim("path_to_qsim_file.qsim", depth=10)
        print(circuit.draw())
    """

    with open(fname, "r", encoding="utf-8") as f:
        data = f.read()
    logger.info(data)
    lines = data.strip().splitlines()

    try:
        qcount = int(lines.pop(0))
    except (IndexError, ValueError) as exc:
        raise ValueError("First line should be qubit count") from exc

    c = models.Circuit(qcount)

    for line in lines:
        fields = line.split()
        if not fields:
            continue
        # qsim only orders the gates of each qubit: a later line can have an earlier time
        if len(fields) < 2:
            raise ValueError(f"{fname}: cannot read the line {line!r}")
        if int(fields[0]) >= depth:
            continue
        if fields[1] not in _QSIM_GATES:
            raise ValueError(f"{fname}: cannot read the qsim gate in {line!r}")
        qubits, angles, build = _QSIM_GATES[fields[1]]
        if len(fields) != 2 + qubits + angles:
            raise ValueError(
                f"{fname}: {fields[1]} takes {qubits} qubit(s) and {angles} angle(s), "
                f"in {line!r}"
            )
        on = [int(q) for q in fields[2 : 2 + qubits]]
        c.add(build(on, [float(a) for a in fields[2 + qubits :]]))
    return c
