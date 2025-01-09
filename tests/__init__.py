#!/usr/bin/env python3
#
#  Copyright 2002-2024 Barcelona Supercomputing Center (www.bsc.es)
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

"""Main Quantum Distributed Library Tests."""

from time import time
import unittest
import numpy as np


class BaseTimedTestCase(unittest.TestCase):
    """Base test class that sets seed and measures the time."""

    def setUp(self):
        np.random.seed(521)
        self.start_time = time()
        self.end_time = self.start_time

    def tearDown(self):
        self.end_time = time()
        print(f"Test {self.id()} took: {self.end_time - self.start_time:.3f} seconds")
