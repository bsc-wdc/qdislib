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

"""Central logging control for Qdislib.

One switch replaces the per-function ``verbose=`` flags that used to be scattered
across the cutting and find_cut routines: raise the package log level and the
diagnostic messages (candidate scores, chosen backend/method, cache hits, ...)
are emitted; leave it alone and the library stays silent (a ``NullHandler`` is
attached in :mod:`Qdislib`).

Usage::

    import Qdislib as qd
    qd.set_log_level("DEBUG")     # see everything
    qd.set_log_level("INFO")      # milestones only
    qd.set_log_level("WARNING")   # back to quiet (the default)

Or set the ``QDISLIB_LOG_LEVEL`` environment variable before import (e.g.
``QDISLIB_LOG_LEVEL=DEBUG``).
"""

import logging
import os

#: The package-root logger; every module logger (``Qdislib.<...>``) is a child of
#: it, so setting its level and handler governs the whole library.
_PACKAGE_LOGGER = "Qdislib"

__all__ = ["set_log_level"]


def set_log_level(level):
    """Enable Qdislib's diagnostic logging at ``level``.

    Sets the level of the ``Qdislib`` package logger and, on first use, attaches a
    stderr handler so messages are actually shown (the library otherwise stays
    silent behind a ``NullHandler``). This is the replacement for the old
    per-function ``verbose=True`` flags -- one call turns diagnostics on for the
    whole library.

    :param level: A level name (``"DEBUG"``, ``"INFO"``, ``"WARNING"``, ...) or a
        :mod:`logging` level constant.
    :return: The configured package logger.
    """
    logger = logging.getLogger(_PACKAGE_LOGGER)
    if isinstance(level, str):
        level = logging.getLevelName(level.upper())
    logger.setLevel(level)
    if not any(getattr(h, "_qdislib_handler", False) for h in logger.handlers):
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(name)s %(levelname)s: %(message)s"))
        handler._qdislib_handler = True  # marker so we attach it only once
        logger.addHandler(handler)
    return logger


def _init_from_env():
    """Apply ``QDISLIB_LOG_LEVEL`` if it is set (called once at package import)."""
    level = os.environ.get("QDISLIB_LOG_LEVEL")
    if level:
        set_log_level(level)
