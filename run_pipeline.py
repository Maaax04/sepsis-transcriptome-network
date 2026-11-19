"""Run the entire pipeline for sepsis transcriptome analysis."""

import subprocess
import sys
import time
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parent/"src"

scripts = [
        "01_fetch_data.py",
        "02_diff_expression.py",
        "03_visualize_de.py",
        "04_ppi_network.py",
        "05_functional_enrichment.py"
    ]

def main() -> None:
    for script in scripts:
        print(f"Running {script}...")
        start_time = time.perf_counter()
        subprocess.run([sys.executable, str(SRC_DIR/script)], check=True)
        end_time = time.perf_counter()
        print(f"Finished {script} in {end_time - start_time:.2f} seconds.")

if __name__ == "__main__":
    main()