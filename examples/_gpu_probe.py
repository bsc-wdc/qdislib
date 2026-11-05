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

"""Backends for :mod:`verify_gpu_pycompss`, recording from inside each task which GPU
COMPSs gave it (a module of its own: a worker has to import these classes)."""

import json
import os

from Qdislib.core.backend.aer_backend import AerBackend
from Qdislib.core.backend.cudaq_backend import CudaqBackend


def _record(kind, **extra):
    row = dict(kind=kind, gpus=os.environ.get("COMPSS_BINDED_GPUS"), **extra)
    with open(os.environ["QDISLIB_PROBE_FILE"], "a") as f:
        f.write(json.dumps(row) + "\n")


class ProbeCpu(AerBackend):
    def evaluate_joint_counts(self, *args, **kwargs):
        _record("cpu")
        return super().evaluate_joint_counts(*args, **kwargs)


class ProbeGpu(AerBackend):
    """Runs on the CPU, but asks for the tasks of a backend on a GPU."""

    uses_gpu = True

    def evaluate_joint_counts(self, *args, **kwargs):
        _record("gpu")
        return super().evaluate_joint_counts(*args, **kwargs)


class ProbeCudaq(CudaqBackend):
    def evaluate_joint_counts(self, *args, **kwargs):
        out = super().evaluate_joint_counts(*args, **kwargs)
        import cudaq

        _record("cudaq", target=cudaq.get_target().name)
        return out
