#!/usr/bin/env python3
"""Run the Python code blocks of the documentation pages, in order, page by page.

    python docs/check_code_blocks.py                   # every page under docs/source
    python docs/check_code_blocks.py guide/shots.rst   # one page

Each page runs in a fresh process, as one script, so a page's blocks share a namespace
and have to work in the order shown. A block that cannot run here (it needs a device or
a GPU) has the comment ``.. not run: <why>`` just before it. Run this with Qdislib and
the SDKs of the pages (qiskit, qibo, and pyzx with lmdb for the cache) importable.
"""

import pathlib
import re
import subprocess
import sys
import tempfile
import time

SOURCE = pathlib.Path(__file__).resolve().parent / "source"
BLOCK = re.compile(r"^\.\. code-block:: python\s*$")


def python_blocks(path):
    """The page's Python blocks, as (line number, code), without the skipped ones."""
    lines = path.read_text().splitlines()
    blocks, skipped = [], 0
    for number, line in enumerate(lines):
        if not BLOCK.match(line):
            continue
        before = next((lines[i] for i in range(number - 1, -1, -1) if lines[i].strip()), "")
        if before.startswith(".. not run:"):
            skipped += 1
            continue
        code = []
        for body in lines[number + 1 :]:
            if body.strip() and not body.startswith("    "):
                break
            code.append(body[4:])
        blocks.append((number + 1, "\n".join(code)))
    return blocks, skipped


def main():
    pages = [SOURCE / name for name in sys.argv[1:]] or sorted(SOURCE.rglob("*.rst"))
    failed = 0
    for page in pages:
        blocks, skipped = python_blocks(page)
        if not blocks:
            continue
        script = "\n".join(
            f"print('--- line {number}')\n{code}" for number, code in blocks
        )
        start = time.time()
        with tempfile.TemporaryDirectory() as directory:
            result = subprocess.run(
                [sys.executable, "-c", "import matplotlib; matplotlib.use('Agg')\n" + script],
                capture_output=True, text=True, cwd=directory,
            )
        status = "ok" if result.returncode == 0 else "FAILED"
        print(f"{status:6s} {page.relative_to(SOURCE)}: {len(blocks)} blocks run, "
              f"{skipped} marked not run, {time.time() - start:.1f} s", flush=True)
        if result.returncode:
            failed += 1
            print(result.stdout[-400:], result.stderr[-1500:], sep="\n")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
