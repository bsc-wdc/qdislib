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
Qdislib CUDA-Q graph utils.

This file contains all auxiliary CUDA-Q graph classes and functions.
"""

import logging

logger = logging.getLogger(__name__)

__all__ = [
    "circuit_cudaq_to_dag",
    "dag_to_circuit_cudaq",
    "dag_to_circuit_cudaq_kernel",
]

import ast
import json
import math
import networkx as nx
import numpy as np

try:
    import cudaq

    cudaq_available = True
except ImportError:
    cudaq = None
    cudaq_available = False

from Qdislib.utils.runtime import task
from Qdislib.utils.circuit_dag import ops_to_dag
from Qdislib.utils.gates import lowered

#: The adjoint of a gate, as the gate it is (``t.adj(q)`` -> ``"TDG"``).
_CUDAQ_ADJOINT = {"S": "SDG", "T": "TDG", "H": "H", "X": "X", "Y": "Y", "Z": "Z"}
#: A measurement Qdislib's final read in Z cannot stand for.
_CUDAQ_OTHER_BASES = ("MX", "MY")
#: Statements that change which gates run: the source alone does not say.
_CUDAQ_CONTROL_FLOW = (ast.For, ast.While, ast.If, ast.With, ast.Try)


def _unreadable(node, why):
    """The error for a kernel statement the DAG cannot hold. Leaving it out would
    cut, and answer for, a different circuit."""
    return ValueError(
        f"Cannot read this CUDA-Q kernel: {why} ({ast.unparse(node)!r}). Cutting reads "
        "a kernel's source, so write each gate on its own line, on a literal qubit "
        "(q[0]) with a literal angle, or build the circuit in another SDK."
    )


def _cudaq_gate_name(call):
    """Upper-cased gate name for a CUDA-Q gate call AST node, or ``None``.

    Handles plain calls (``h(...)`` -> ``"H"``), the controlled-gate attribute
    syntax (``x.ctrl(c, t)`` -> ``"CX"``, ``z.ctrl(c, t)`` -> ``"CZ"``) and adjoints
    (``t.adj(q)`` -> ``"TDG"``). Returns ``None`` for what only allocates qubits
    (``cudaq.qvector(...)``).

    :raises ValueError: For any other attribute call (``cudaq.control(...)``, ...).
    """
    func = call.func
    if isinstance(func, ast.Name):
        return func.id.upper()
    if isinstance(func, ast.Attribute) and isinstance(func.value, ast.Name):
        base = func.value.id.upper()
        if func.attr == "ctrl":
            return {"X": "CX", "Z": "CZ"}.get(base, "C" + base)
        if func.attr == "adj" and base in _CUDAQ_ADJOINT:
            return _CUDAQ_ADJOINT[base]
        if func.value.id == "cudaq" and func.attr in ("qvector", "qubit"):
            return None
    raise _unreadable(call, "a call that is not a plain gate")


def _eval_cudaq_param(node):
    """Best-effort evaluation of a numeric gate-parameter expression node.

    Supports literals and expressions such as ``np.pi / 2`` / ``math.pi``;
    returns ``None`` if the value cannot be resolved to a float.
    """
    try:
        value = eval(  # noqa: S307 - trusted kernel source, restricted namespace
            compile(ast.Expression(node), "<cudaq-param>", "eval"),
            {"np": np, "numpy": np, "math": math, "pi": math.pi},
        )
        return float(value)
    except Exception:
        return None


def _cudaq_allocation(statement):
    """``(register name, its size)`` for ``a = cudaq.qvector(n)`` / ``cudaq.qubit()``,
    else ``None``. The name is ``None`` for a bare allocation; the size is ``None``
    where the source does not state it as a number."""
    target = None
    if isinstance(statement, ast.Assign) and len(statement.targets) == 1:
        target = statement.targets[0]
    elif isinstance(statement, ast.AnnAssign):
        target = statement.target
    elif not isinstance(statement, ast.Expr):
        return None
    call = statement.value
    if not isinstance(call, ast.Call):
        return None
    func = call.func
    if isinstance(func, ast.Attribute) and isinstance(func.value, ast.Name):
        kind = func.attr if func.value.id == "cudaq" else None
    else:
        kind = func.id if isinstance(func, ast.Name) else None
    if kind not in ("qvector", "qubit"):
        return None
    name = target.id if isinstance(target, ast.Name) else None
    if kind == "qubit":
        return name, 1
    size = _eval_cudaq_param(call.args[0]) if call.args else None
    return name, int(size) if size is not None and size == int(size) else None


def _cudaq_qubit(call, arg, registers):
    """The global qubit of ``arg``, a register's ``q[i]``: registers are numbered in
    the order they are allocated, as the kernel's state is."""
    index = arg.slice.value if isinstance(arg.slice, ast.Constant) else None
    if not isinstance(index, int):
        raise _unreadable(call, "a qubit that is not a literal index")
    if not isinstance(arg.value, ast.Name) or arg.value.id not in registers:
        raise _unreadable(call, "a qubit of a register the kernel does not allocate")
    first, size = registers[arg.value.id]
    if size is not None and not 0 <= index < size:
        raise _unreadable(call, "a qubit outside its register")
    return first + index


def _cudaq_source(kernel):
    """The Python source of a ``@cudaq.kernel``, ``None`` for a builder's kernel."""
    source = getattr(kernel, "funcSrc", None)  # CUDA-Q >= 0.16 has no to_json()
    if source is None and hasattr(kernel, "to_json"):
        payload = kernel.to_json()
        if isinstance(payload, str):
            payload = json.loads(payload)
        source = payload["funcSrc"]
    return source


def _cudaq_listing(kernel):
    """A kernel's source, statement by statement.

    :param kernel: A CUDA-Q kernel with source (``PyKernelDecorator``).
    :raises ValueError: See :func:`circuit_cudaq_to_dag`.
    :return: ``(operations, width)``: ``(label, name, qubits, params, is_measurement)``
        per statement, in order, and the qubits allocated (``None`` if not stated).
    """
    tree = ast.parse(_cudaq_source(kernel))
    func = next(
        (node for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)), None
    )
    if func is None:
        return [], None

    operations = []
    registers = {}  # name -> (its first qubit, its size or None)
    next_qubit = 0  # None after a register of no stated size
    measured = []  # (first, end) of the qubits measured so far
    for statement in func.body:
        if isinstance(statement, _CUDAQ_CONTROL_FLOW):
            raise _unreadable(statement, "a loop or a branch")
        allocation = _cudaq_allocation(statement)
        if allocation is not None:
            if next_qubit is None:
                raise _unreadable(statement, "a register after one of no stated size")
            name, size = allocation
            if name is not None:
                registers[name] = (next_qubit, size)
            next_qubit = None if size is None else next_qubit + size
            continue
        # Gate operations are bare calls (``h(q[0])``) or measurement
        # assignments (``m = mz(q[0])``); everything else is skipped.
        if isinstance(statement, ast.Expr) and isinstance(statement.value, ast.Call):
            call = statement.value
        elif isinstance(statement, ast.Assign) and isinstance(
            statement.value, ast.Call
        ):
            call = statement.value
        else:
            continue

        name = _cudaq_gate_name(call)
        if name is None:
            continue
        if name in _CUDAQ_OTHER_BASES:
            raise _unreadable(call, "a measurement not in the Z basis")
        is_measurement = name == "MZ"

        # Subscript arguments (``q[0]``) are qubits; the rest are parameters.
        qubits, parameters = [], []
        for arg in call.args:
            if isinstance(arg, ast.Subscript):
                qubits.append(_cudaq_qubit(call, arg, registers))
            elif is_measurement and isinstance(arg, ast.Name):  # mz(q): the register
                if arg.id not in registers:
                    raise _unreadable(call, "a register the kernel does not allocate")
                first, size = registers[arg.id]
                end = math.inf if size is None else first + size
                measured.append((first, end))
                qubits.extend(range(first, first + (size or 0)))
            else:
                value = _eval_cudaq_param(arg)
                if value is None:
                    raise _unreadable(call, "an argument that is not a number or q[i]")
                parameters.append(value)

        if is_measurement:
            measured.extend((qubit, qubit + 1) for qubit in qubits)
        elif not qubits:  # no qubit named: not a gate of this circuit
            continue
        elif any(first <= q < end for q in qubits for first, end in measured):
            raise _unreadable(call, "a gate after a measurement of its qubit")
        label = f"{name}_{len(operations) + 1}"
        operations.append((label, name, qubits, parameters, is_measurement))
    return operations, next_qubit


def circuit_cudaq_to_dag(kernel, obs_I=None) -> nx.DiGraph:
    """Convert a CUDA-Q kernel into a Qdislib DAG.

    The kernel's Python source (``kernel.funcSrc``, or ``kernel.to_json()["funcSrc"]``
    before CUDA-Q 0.16) is parsed with
    the :mod:`ast` module rather than scanned with regular expressions, so it is
    robust to formatting/comments, the qubit-register variable name, expression
    parameters (e.g. ``np.pi / 2``) and the ``gate.ctrl(...)`` controlled syntax.

    Measurements that end their qubits are left out, keeping the numbering, as for the
    other SDKs: Qdislib measures every qubit at the end anyway.

    :param kernel: A CUDA-Q kernel (``PyKernelDecorator``).
    :param obs_I: Optional observable string; each ``"I"`` appends an identity
        observable sink node on that qubit (internal use when cutting).
    :raises ValueError: For what the source does not say gate by gate: a loop or a
        branch, a gate on a whole register or on a computed qubit, a register after
        one of no stated size, an angle that is not a number there (a kernel argument,
        a variable); and for what a final read in Z cannot stand for: a gate after a
        measurement of its qubit, a measurement in another basis.
    :return: A :class:`networkx.DiGraph` whose nodes carry ``gate`` (upper-cased
        name), ``qubits`` (tuple of ints) and ``parameters`` (list of floats),
        with blue (nearest-predecessor) and red (other-predecessor) edges.
    """
    operations, _ = _cudaq_listing(kernel)
    ops = [op[:4] for op in operations if not op[4]]
    return ops_to_dag(ops, obs_I, sdk="cudaq")


def _measures(kernel):
    """Whether a kernel's own source has a measurement (a builder's kernel has no
    source: only Qdislib builds those, with the measurements a cut needs)."""
    source = _cudaq_source(kernel)
    if source is None:
        return False
    return any(
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id.upper() in ("MZ", *_CUDAQ_OTHER_BASES)
        for node in ast.walk(ast.parse(source))
    )


def cudaq_without_measurements(kernel):
    """``kernel`` as CUDA-Q should run it for Qdislib: without its final measurements.

    ``cudaq.get_state`` collapses at a measurement, and ``cudaq.sample`` reports only
    the qubits measured (when some are) and does not collapse a mid-circuit one. A
    kernel with none is returned as it is.

    :param kernel: A CUDA-Q kernel.
    :raises ValueError: See :func:`circuit_cudaq_to_dag`.
    :return: A kernel on the same qubits.
    """
    if not _measures(kernel):
        return kernel
    operations, width = _cudaq_listing(kernel)
    if width is None:
        width = max((q for op in operations for q in op[2]), default=-1) + 1
    return dag_to_circuit_cudaq_kernel(circuit_cudaq_to_dag(kernel), width)[0]


@task(returns=1)
def dag_to_circuit_cudaq(dag: nx.DiGraph, num_qubits: int):
    if dag is None:
        return None

    topo_order = list(nx.topological_sort(dag))

    # Extract "Observable I" nodes if present
    obs_I = []
    for node in topo_order:
        if dag.nodes[node]["gate"].upper() == "OBSERVABLE I":
            obs_I.append(dag.nodes[node]["qubits"][0])
            dag.remove_node(node)

    observables = ["Z"] * num_qubits
    for element in obs_I:
        observables[element] = "I"
    return [(dag, num_qubits), "".join(observables)]


# @task(returns=1)
def dag_to_circuit_cudaq_kernel(dag: nx.DiGraph, num_qubits: int, explicit_final=False):
    """
    Convert a DAG into a CUDA-Q kernel using the modern builder API.

    :param dag: networkx.DiGraph whose nodes carry ``gate`` (str), ``qubits``
        (tuple of indices) and ``parameters`` (list) attributes.
    :param num_qubits: Total number of qubits.
    :param explicit_final: When True, append an explicit ``mz`` on every qubit at
        the end and return them so a ``cudaq.sample(..., explicit_measurements=True)``
        run yields, per shot, the mid-circuit cut outcomes followed by the final
        measurements. Needed for the gate-cut identity-observable sign recovery
        (see :func:`~Qdislib.core.cutting_algorithms.wire_cutting._expec_value_cudaq`).
    :return: ``[kernel, observables]``, or ``[kernel, observables, measured_order]``
        when ``explicit_final`` is True (``measured_order`` = cut-measured qubits
        in ``mz``-call order).
    """
    if dag is None:
        return None

    # Extract "Observable I" nodes (identity-observed qubits) before building the
    # kernel: they are not gates, and their positions become "I" in the per-qubit
    # observable string returned alongside the kernel.
    obs_I = []
    for node in list(nx.topological_sort(dag)):
        if dag.nodes[node]["gate"].upper() == "OBSERVABLE I":
            obs_I.append(dag.nodes[node]["qubits"][0])
            dag.remove_node(node)
    observables = ["Z"] * num_qubits
    for qubit in obs_I:
        observables[qubit] = "I"
    observables = "".join(observables)

    if explicit_final:
        # cudaq.sample does NOT collapse mid-circuit measurements, so the gate-cut
        # measure-and-prepare path (where a cut qubit must collapse to its outcome
        # and keep evolving) is emitted as a @cudaq.kernel and executed with
        # cudaq.run, which does collapse. See _build_cudaq_run_kernel.
        return _build_cudaq_run_kernel(dag, num_qubits, observables)

    # Builder path (used only for subcircuits with no mid-circuit measurement, run
    # under cudaq.sample where deferred measurement is fine).
    # Initialize CUDA-Q kernel
    kernel = cudaq.make_kernel()

    # Allocate qubits explicitly
    q = [kernel.qalloc() for _ in range(num_qubits)]
    # Traverse DAG
    for node in nx.topological_sort(dag):
        node_data = dag.nodes[node]
        for gate, positions, params in lowered(
            node_data["gate"], node_data["parameters"] or ()
        ):
            gate = gate.upper()
            qubits = node_data["qubits"]
            if positions is not None:
                qubits = tuple(qubits[k] for k in positions)

            # --- Common 1-qubit gates ---
            if gate in ["H", "X", "Y", "Z", "T", "S", "SDG"]:
                for qb in qubits:
                    getattr(kernel, gate.lower())(q[qb])

            # --- Rotation gates ---
            elif gate in ["RX", "RY", "RZ"]:
                angle = params[0] if params else 0.0
                for qb in qubits:
                    getattr(kernel, gate.lower())(angle, q[qb])

            # --- General single-qubit unitary / phase ---
            elif gate == "U3":
                for qb in qubits:
                    kernel.u3(params[0], params[1], params[2], q[qb])

            elif gate in ["U1", "CPHASE_SINGLE"]:  # r1 = phase gate diag(1, e^{i lam})
                for qb in qubits:
                    kernel.r1(params[0], q[qb])

            # --- Controlled gates ---
            elif gate in ["CNOT", "CX"]:
                ctrl, tgt = qubits
                kernel.cx(q[ctrl], q[tgt])

            elif gate == "CZ":
                ctrl, tgt = qubits
                kernel.cz(q[ctrl], q[tgt])

            elif gate == "CU3":
                # Controlled SU(2) U3. cudaq's cu3 applies the OpenQASM u3 (an extra
                # e^{i(phi+lam)/2} global phase) on the control-1 block; cancel it with a
                # phase (r1) on the control so the block is exactly qibo's SU(2) U3.
                ctrl, tgt = qubits
                theta, phi, lam = params[0], params[1], params[2]
                kernel.cu3(theta, phi, lam, q[ctrl], q[tgt])
                kernel.r1(-(phi + lam) / 2, q[ctrl])

            elif gate in ["CU1", "CPHASE"]:  # controlled phase on |11>
                ctrl, tgt = qubits
                kernel.cr1(params[0], q[ctrl], q[tgt])

            elif gate == "SWAP":
                kernel.swap(q[qubits[0]], q[qubits[1]])

            elif (
                gate == "RZZ"
            ):  # exp(-i theta/2 ZZ) = CX . RZ(theta) . CX (no native rzz)
                ctrl, tgt = qubits
                kernel.cx(q[ctrl], q[tgt])
                kernel.rz(params[0], q[tgt])
                kernel.cx(q[ctrl], q[tgt])

            elif gate == "ISWAP":  # no native iswap; standard 6-gate decomposition
                a, b = qubits
                kernel.s(q[a])
                kernel.s(q[b])
                kernel.h(q[a])
                kernel.cx(q[a], q[b])
                kernel.cx(q[b], q[a])
                kernel.h(q[b])

            # --- Measurement gates ---
            elif gate in ["MZ", "MX", "MY"]:
                for qb in qubits:
                    getattr(kernel, gate.lower())(q[qb])

            # --- Gate-cut MEASURE: a projective Z measurement that collapses the
            # --- qubit to its outcome |s> and lets it CONTINUE (the CZ cut's
            # --- (I +/- Z)/2 projector, Mitarai-Fujii arXiv:1909.07534). The outcome
            # --- is recovered from the explicit-measurements record for the
            # --- reconstruction sign; the qubit is NOT reset (resetting it discarded
            # --- the collapse and broke reconstruction for gates after the cut).
            elif gate == "MEASURE":
                for qb in qubits:
                    kernel.mz(q[qb])

            else:
                raise ValueError(f"Unsupported gate type: {gate}")

    return [kernel, observables]


# Compiled run-kernels cached by source hash: cudaq JIT-compiles each kernel and
# many gate-cut configurations share the same (renumbered) subcircuit.
_CUDAQ_RUN_KERNEL_CACHE = {}
#: A long-lived worker builds many: the oldest are dropped past this many.
_CUDAQ_RUN_KERNEL_CACHE_SIZE = 256
_CUDAQ_ONE_QUBIT = {"H": "h", "X": "x", "Y": "y", "Z": "z", "S": "s", "T": "t"}


def _compile_cudaq_source(source, func_name="cut_kernel"):
    """Compile a generated ``@cudaq.kernel`` source string into a callable kernel.

    cudaq's decorator needs real source (``inspect.getsource``), so the generated
    kernel is written to a temp module file and imported. Cached by source hash to
    avoid recompiling identical subcircuits.
    """
    import hashlib
    import importlib.util
    import os
    import tempfile

    key = hashlib.sha256(source.encode()).hexdigest()
    cached = _CUDAQ_RUN_KERNEL_CACHE.get(key)
    if cached is not None:
        return cached
    # The decorator reads the source when it runs: the file is not needed after that,
    # and a worker builds many kernels.
    with tempfile.TemporaryDirectory(prefix="qdislib_cudaq_") as directory:
        path = os.path.join(directory, f"{func_name}.py")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(source)
        spec = importlib.util.spec_from_file_location(f"{func_name}_{key}", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    kernel = getattr(module, func_name)
    if len(_CUDAQ_RUN_KERNEL_CACHE) >= _CUDAQ_RUN_KERNEL_CACHE_SIZE:
        _CUDAQ_RUN_KERNEL_CACHE.pop(next(iter(_CUDAQ_RUN_KERNEL_CACHE)))
    _CUDAQ_RUN_KERNEL_CACHE[key] = kernel
    return kernel


def _build_cudaq_run_kernel(dag, num_qubits, observables):
    """Build a ``@cudaq.kernel`` (for ``cudaq.run``) returning the gate-cut record.

    ``cudaq.sample`` does not collapse mid-circuit measurements; ``cudaq.run``
    does. The measure-and-prepare gate-cut path needs the cut qubit to collapse to
    its outcome and keep evolving, so the subcircuit is emitted as a kernel that
    returns the per-shot record ``[cut-point outcomes..., final outcomes
    q0..q(n-1)]`` packed into an int (bit k = outcome k -- cudaq 0.12+ rejects
    ``list`` returns from ``cudaq.run`` entry-point kernels) and is executed with
    ``cudaq.run``. Returns ``[kernel, observables, measured_order]``
    (``measured_order`` = cut-measured qubits in ``mz``-call order); the host
    unpacks bit k with ``(record >> k) & 1`` so the reconstruction is unchanged.
    """
    lines = []
    measured_vars = []
    measured_order = []
    for node in nx.topological_sort(dag):
        data = dag.nodes[node]
        for gate, positions, params in lowered(data["gate"], data["parameters"] or ()):
            gate = gate.upper()
            qubits = data["qubits"]
            if positions is not None:
                qubits = tuple(qubits[k] for k in positions)
            if gate in _CUDAQ_ONE_QUBIT:
                for qb in qubits:
                    lines.append(f"    {_CUDAQ_ONE_QUBIT[gate]}(q[{qb}])")
            elif gate == "SDG":
                for qb in qubits:
                    lines.append(f"    s.adj(q[{qb}])")
            elif gate == "TDG":
                for qb in qubits:
                    lines.append(f"    t.adj(q[{qb}])")
            elif gate in ("RX", "RY", "RZ"):
                angle = float(params[0]) if params else 0.0
                for qb in qubits:
                    lines.append(f"    {gate.lower()}({angle!r}, q[{qb}])")
            elif gate == "U3":
                theta, phi, lam = (float(p) for p in params)
                for qb in qubits:
                    lines.append(f"    u3({theta!r}, {phi!r}, {lam!r}, q[{qb}])")
            elif gate == "U1":
                for qb in qubits:
                    lines.append(f"    r1({float(params[0])!r}, q[{qb}])")
            elif gate in ("CNOT", "CX"):
                ctrl, tgt = qubits
                lines.append(f"    x.ctrl(q[{ctrl}], q[{tgt}])")
            elif gate == "CZ":
                ctrl, tgt = qubits
                lines.append(f"    z.ctrl(q[{ctrl}], q[{tgt}])")
            elif gate == "MEASURE":
                for qb in qubits:
                    var = f"m{len(measured_vars)}"
                    lines.append(f"    {var} = mz(q[{qb}])")
                    measured_vars.append(var)
                    measured_order.append(qb)
            else:
                raise ValueError(f"Unsupported gate type: {gate}")

    final_vars = [f"f{qb}" for qb in range(num_qubits)]
    for qb in range(num_qubits):
        lines.append(f"    f{qb} = mz(q[{qb}])")

    # cudaq 0.12+ dropped list returns from cudaq.run entry-point kernels, so pack
    # the record into an int: bit k holds outcome k of [cut-points..., finals...].
    # The host unpacks with (record >> k) & 1 (see _joint_counts_cudaq /
    # _expectation_cudaq_decoupled); measured_order maps the cut-point bits back to
    # qubits, exactly as the previous list layout did.
    ordered_vars = measured_vars + final_vars
    lines.append("    r = 0")
    for k, var in enumerate(ordered_vars):
        lines.append(f"    if {var}:")
        lines.append(f"        r = r + {1 << k}")

    body = "\n".join(lines) if lines else "    pass"
    source = (
        "import cudaq\n\n"
        "@cudaq.kernel\n"
        "def cut_kernel() -> int:\n"
        f"    q = cudaq.qvector({num_qubits})\n"
        f"{body}\n"
        "    return r\n"
    )
    kernel = _compile_cudaq_source(source)
    return [kernel, observables, measured_order]


# --- SDK table adapters (see Qdislib.utils.sdks) ------------------------------


def cudaq_gates(kernel):
    """``[(label, qubits, is_measurement), ...]`` in kernel order (CUDA-Q has no gate
    iteration API: the kernel source is parsed, as for the DAG)."""
    return [
        (name, tuple(qubits), is_measurement)
        for _, name, qubits, _, is_measurement in _cudaq_listing(kernel)[0]
    ]


def cudaq_num_qubits(kernel):
    """The highest qubit a gate touches, plus one (a kernel states no width)."""
    qubits = [q for _, qs, _ in cudaq_gates(kernel) for q in qs]
    return max(qubits) + 1 if qubits else 0


def cudaq_statevector(kernel):
    """The exact state, qubit 0 the most significant bit (CUDA-Q's is reversed)."""
    import cudaq

    kernel = cudaq_without_measurements(kernel)
    data = np.asarray(cudaq.get_state(kernel), dtype=complex).reshape(-1)
    # From the state: a kernel can allocate more qubits than its gates touch.
    n = int(round(math.log2(len(data))))
    return data.reshape([2] * n).transpose(list(range(n))[::-1]).reshape(-1)
