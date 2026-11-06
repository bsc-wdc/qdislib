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


"""``import Qdislib`` hides CUDA-Q's two notices, and only those."""

import os
import subprocess
import sys
import tempfile
import textwrap
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

FAKE_CUDAQ = textwrap.dedent("""
    import warnings

    warnings.warn("cudaq-logical is in preview. Its APIs may change.", UserWarning)
    warnings.warn(
        "The CUDA-Q `sample` and `observe` algorithmic primitives will change in a "
        "future release.",
        FutureWarning,
    )
    warnings.warn("another notice from CUDA-Q", UserWarning)
    """)


class ImportWarningsTest(unittest.TestCase):
    def test_only_the_two_notices_are_hidden(self):
        with tempfile.TemporaryDirectory() as folder:
            os.mkdir(os.path.join(folder, "cudaq"))
            with open(os.path.join(folder, "cudaq", "__init__.py"), "w") as module:
                module.write(FAKE_CUDAQ)
            env = dict(os.environ, PYTHONPATH=os.pathsep.join([folder, ROOT]))
            run = subprocess.run(
                [sys.executable, "-W", "default", "-c", "import Qdislib"],
                capture_output=True,
                text=True,
                env=env,
                cwd=folder,
                timeout=120,
            )
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertNotIn("cudaq-logical", run.stderr)
        self.assertNotIn("algorithmic primitives", run.stderr)
        self.assertIn("another notice from CUDA-Q", run.stderr)


if __name__ == "__main__":
    unittest.main()
