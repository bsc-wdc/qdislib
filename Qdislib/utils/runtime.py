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
PyCOMPSs runtime shim.

Exposes the PyCOMPSs decorators and helpers (``task``, ``compss_wait_on``,
``constraint``, ``implement``, the parameter directions) when PyCOMPSs is
installed, and transparent no-op fallbacks for serial execution otherwise.

Centralising this here means the "PyCOMPSs not installed" handling lives in one
place (instead of being duplicated in every module), the warning is emitted once,
and tests can obtain a working ``compss_wait_on`` without importing ``pycompss``
directly.
"""

import logging

logger = logging.getLogger(__name__)

__all__ = [
    "task",
    "compss_wait_on",
    "constraint",
    "implement",
    "COLLECTION_IN",
    "COLLECTION_OUT",
    "INOUT",
    "pycompss_available",
]

try:
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
    from pycompss.api.parameter import COLLECTION_IN, COLLECTION_OUT, INOUT
    from pycompss.api.constraint import constraint
    from pycompss.api.implement import implement

    pycompss_available = True
except ImportError:
    logger.warning("PyCOMPSs is not installed; running serially (no parallelism).")

    pycompss_available = False

    def task(*args, **kwargs):
        """No-op stand-in for ``pycompss.api.task.task`` (serial execution)."""

        def decorator(func):
            return func

        return decorator

    def compss_wait_on(obj):
        """No-op stand-in for ``compss_wait_on`` — returns the object unchanged."""
        return obj

    def constraint(*args, **kwargs):
        """No-op stand-in for ``pycompss.api.constraint.constraint``."""

        def decorator(func):
            return func

        return decorator

    def implement(*args, **kwargs):
        """No-op stand-in for ``pycompss.api.implement.implement``."""

        def decorator(func):
            return func

        return decorator

    COLLECTION_IN = COLLECTION_OUT = INOUT = None
