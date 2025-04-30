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

"""Utils Graph Unit tests."""

from tests import BaseTimedTestCase


class GraphTest(BaseTimedTestCase):

    def test_circuit_to_dag(self):
        from Qdislib.utils.graph_qibo import circuit_qibo_to_dag
        # Do something to check that retrieves a valid dag from circuit
        self.assertTrue(True)

    def test_dag_to_circuit(self):
        from Qdislib.utils.graph_qibo import dag_to_circuit_qibo
        # Do something to check that the dag_to_circuit_qibo function works
        self.assertTrue(True)
