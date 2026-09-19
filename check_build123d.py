# /// script
# requires-python = ">=3.12"
# dependencies = ["build123d"]
# ///
import build123d, sys
from build123d import Box
b = Box(1, 2, 3)
print(f"build123d {build123d.__version__} on python {sys.version.split()[0]}")
print(f"box volume {b.volume:.1f}  faces {len(b.faces())}")
