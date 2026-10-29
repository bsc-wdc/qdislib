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

"""Free parameters of a cut circuit (Qiskit ``Parameter`` objects): the cut is planned
once with them symbolic, and each set of values is bound where a fragment circuit
is built (``subs.bind(values)``, ``parameter_values=``).

Values are kept by parameter *name*, so they survive the trip to a worker.
"""


def is_symbolic(value):
    """True for a parameter expression that still has free parameters."""
    return bool(getattr(value, "parameters", None))


def _names(value):
    return {p.name for p in value.parameters}


def free_parameters(dag):
    """The names of the free parameters in a DAG's gates."""
    names = set()
    for _, data in dag.nodes(data=True):
        for value in data.get("parameters") or ():
            if is_symbolic(value):
                names |= _names(value)
    return names


def circuit_parameters(circuit, dag):
    """The circuit's free parameters, by name, in the order a list of values
    follows: Qiskit's ``circuit.parameters`` order, else sorted by name."""
    names = free_parameters(dag)
    if not names:
        return ()
    ordered = getattr(circuit, "parameters", None)
    try:
        order = [p.name for p in ordered]
    except TypeError:
        order = []
    if set(order) == names:
        return tuple(order)
    return tuple(sorted(names))


def check_numeric(dag, gate_names):
    """Raise unless every cut gate has numeric parameters: a cut's terms and
    weights are fixed when the cut is planned."""
    symbolic = [
        name
        for name in gate_names
        if name in dag
        and any(is_symbolic(v) for v in dag.nodes[name].get("parameters") or ())
    ]
    if symbolic:
        raise ValueError(
            f"Cut gate(s) {symbolic} have free parameters: a cut gate's angle sets "
            "the cut's terms and weights, so it must be a number. Free parameters "
            "can be anywhere else in the circuit."
        )


def normalise_values(values, parameters):
    """``{name: float}`` for every parameter of ``parameters`` (names, in order),
    from a mapping (``Parameter`` or name keys) or a sequence in that order."""
    if not parameters:
        raise ValueError("this circuit has no free parameters to bind")
    if hasattr(values, "items"):
        bound = {getattr(k, "name", k): float(v) for k, v in values.items()}
    else:
        values = list(values)
        if len(values) != len(parameters):
            raise ValueError(
                f"expected {len(parameters)} values, one per parameter "
                f"{list(parameters)}, got {len(values)}"
            )
        bound = {name: float(v) for name, v in zip(parameters, values)}
    missing = [p for p in parameters if p not in bound]
    unknown = [p for p in bound if p not in parameters]
    if missing or unknown:
        raise ValueError(
            "values must bind every free parameter and nothing else: "
            + "; ".join(
                part
                for part in (
                    f"missing {missing}" if missing else "",
                    f"unknown {unknown}" if unknown else "",
                )
                if part
            )
        )
    return bound


def _bound(value, values):
    if not is_symbolic(value):
        return value
    from qiskit.circuit import Parameter

    if isinstance(value, Parameter):  # a plain parameter: no expression to bind
        return values[value.name]
    result = value.bind({p: values[p.name] for p in value.parameters})
    return float(result)


def bind_dag(dag, values):
    """Bind ``values`` (``{name: float}``) in place, in every gate of ``dag``."""
    for _, data in dag.nodes(data=True):
        params = data.get("parameters")
        if params and any(is_symbolic(v) for v in params):
            data["parameters"] = tuple(_bound(v, values) for v in params)


def require_bound(dag):
    """Raise if a gate of ``dag`` still has free parameters (before running it)."""
    left = free_parameters(dag)
    if left:
        raise ValueError(
            f"this subcircuit has free parameters {sorted(left)}: give their values "
            "with subs.bind(values), or parameter_values= on the cutting function"
        )
