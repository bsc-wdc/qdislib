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
Qdislib exception utils.

This file contains all auxiliary exceptions.

:class:`QdislibException` is the base class for every Qdislib-specific error, so
``except QdislibException`` catches them all. The concrete errors also subclass
the matching built-in (e.g. :class:`ValueError`), so existing
``except ValueError`` handlers keep working.
"""

__all__ = [
    "QdislibException",
    "UnsupportedBackendError",
    "InvalidCutError",
]


class QdislibException(Exception):
    """Base class for all Qdislib-specific exceptions."""

    def __init__(self, message: str) -> None:
        """Create a new QdislibException instance.

        :param message: Exception message.
        """
        msg = f"Qdislib Exception: {message}"
        super().__init__(msg)


class UnsupportedBackendError(QdislibException, ValueError):
    """Raised when an unknown quantum-software backend (``software``) is requested.

    Also a :class:`ValueError` for backward compatibility.
    """


class InvalidCutError(QdislibException, ValueError):
    """Raised when a cut specification (gates to cut or wire cuts) is malformed.

    Also a :class:`ValueError` for backward compatibility.
    """
