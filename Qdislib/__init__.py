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

"""Main Quantum Distributed Library package."""

import logging

__version__ = "1.0.0"

# Attach a NullHandler so the library never emits "no handlers" warnings and
# stays silent unless the host application configures logging itself. Diagnostic
# output is opt-in via qd.set_log_level(...) or the QDISLIB_LOG_LEVEL env var.
logging.getLogger(__name__).addHandler(logging.NullHandler())

from Qdislib.utils import logging_config as _logging_config

_logging_config._init_from_env()

from Qdislib.api import *
