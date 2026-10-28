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

"""A strict stand-in for PyCOMPSs, to test the ``sync=False`` paths without it.

Without PyCOMPSs, Qdislib's tasks are plain functions, so a value computed on the
master from futures (``sum(results)``) goes unnoticed. Here a task returns an opaque
:class:`Future` -- no arithmetic, no iteration, no ``float`` -- that only
``compss_wait_on`` or another task can open, as under PyCOMPSs:

- a task's arguments are resolved when they are futures; futures *inside* a list
  only when that parameter is ``COLLECTION_IN``;
- a task called from inside a task runs inline (a worker has no nested tasks).

Call :func:`install` before importing Qdislib.
"""

import functools
import inspect
import sys
import types

COLLECTION_IN = "COLLECTION_IN"
COLLECTION_OUT = "COLLECTION_OUT"
INOUT = "INOUT"

#: How many tasks ran, for the tests.
calls = {"tasks": 0}
_depth = [0]


class Future:
    """An unresolved task result."""

    __slots__ = ("_value",)

    def __init__(self, value):
        self._value = value

    def __repr__(self):
        return "<Future>"


def _open(value, deep=False):
    if isinstance(value, Future):
        value = value._value
    if deep and isinstance(value, (list, tuple)):
        # rebuilt as a plain list / tuple, as PyCOMPSs does
        return (list if isinstance(value, list) else tuple)(
            _open(v, deep) for v in value
        )
    if deep and isinstance(value, dict):
        # PyCOMPSs takes a COLLECTION_IN apart, nested containers too, and builds
        # a dict back as a plain dict: a dict subclass loses its class (seen on
        # MareNostrum 5, where a lazily computed dict arrived empty and plain)
        return {k: _open(v, deep) for k, v in value.items()}
    return value


def compss_wait_on(obj):
    if isinstance(obj, list):
        return [compss_wait_on(v) for v in obj]
    return _open(obj)


def task(**options):
    returns = options.pop("returns", 1)
    collections = {k for k, v in options.items() if v == COLLECTION_IN}

    def decorator(func):
        signature = inspect.signature(func)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            if _depth[0]:  # inside a task: run inline, as on a worker
                return func(*args, **kwargs)
            bound = signature.bind(*args, **kwargs)
            for name, value in bound.arguments.items():
                bound.arguments[name] = _open(value, deep=name in collections)
            calls["tasks"] += 1
            _depth[0] += 1
            try:
                value = func(*bound.args, **bound.kwargs)
            finally:
                _depth[0] -= 1
            count = returns
            if isinstance(returns, str) and returns.startswith("value_of("):
                # PyCOMPSs: the number of returns read from a parameter at call time
                count = int(bound.arguments[returns[len("value_of("):-1].strip()])
            if isinstance(count, int) and count > 1:
                return tuple(Future(v) for v in value)
            return Future(value)

        return wrapper

    return decorator


def _passthrough(*args, **kwargs):
    return lambda func: func


def install():
    """Register the fake ``pycompss`` modules (before Qdislib is imported)."""
    if "Qdislib" in sys.modules:
        raise RuntimeError("install() must run before Qdislib is imported")
    modules = {}
    for name in (
        "pycompss",
        "pycompss.api",
        "pycompss.api.task",
        "pycompss.api.api",
        "pycompss.api.parameter",
        "pycompss.api.constraint",
        "pycompss.api.implement",
    ):
        modules[name] = sys.modules[name] = types.ModuleType(name)
    modules["pycompss.api.task"].task = task
    modules["pycompss.api.api"].compss_wait_on = compss_wait_on
    parameter = modules["pycompss.api.parameter"]
    parameter.COLLECTION_IN = COLLECTION_IN
    parameter.COLLECTION_OUT = COLLECTION_OUT
    parameter.INOUT = INOUT
    modules["pycompss.api.constraint"].constraint = _passthrough
    modules["pycompss.api.implement"].implement = _passthrough
