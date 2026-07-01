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

"""Run the Qdislib QPU listener on the BSC quantum cluster.

At BSC the classical (PyCOMPSs) cut job and the quantum-control job run on
separate clusters that can't share one scheduler allocation, so they communicate
through a shared filesystem directory. Launch this as its own job on the quantum
cluster: it watches ``$QDISLIB_SUBCIRCUITS_DIR`` for subcircuits submitted by a
Qdislib cut workflow on the classical cluster, executes each on the QPU (Blue or
Red) via qililab, and writes the results back.

    # both clusters must agree on the shared directory:
    export QDISLIB_SUBCIRCUITS_DIR=/gpfs/projects/.../subcircuits
    export RUNCARD_RED=/path/to/red_runcard.yml     # or RUNCARD_BLUE

    python examples/qpu_listener.py --qpu red --nshots 1000

Stop it by creating a ``STOP`` file in the shared directory, with Ctrl-C, or with
``--idle-timeout`` (exit after N idle seconds).
"""

import argparse

import Qdislib as qd


def main():
    parser = argparse.ArgumentParser(description="Qdislib QPU listener")
    parser.add_argument(
        "--qpu",
        default="red",
        choices=["red", "blue"],
        help="Which BSC QPU to use ($RUNCARD_RED / $RUNCARD_BLUE).",
    )
    parser.add_argument("--dir", default=None, help="Shared hand-off directory.")
    parser.add_argument("--runcard", default=None, help="Explicit runcard path.")
    parser.add_argument("--nshots", type=int, default=1000, help="Shots per circuit.")
    parser.add_argument(
        "--idle-timeout",
        type=float,
        default=None,
        help="Exit after this many idle seconds (default: run until STOP/Ctrl-C).",
    )
    parser.add_argument("--log-level", default="INFO", help="Qdislib log level.")
    args = parser.parse_args()

    qd.set_log_level(args.log_level)
    qd.run_qpu_listener(
        directory=args.dir,
        qpu=args.qpu,
        runcard=args.runcard,
        nshots=args.nshots,
        idle_timeout=args.idle_timeout,
    )


if __name__ == "__main__":
    main()
