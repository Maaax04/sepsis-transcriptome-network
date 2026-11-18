"""
Downloads and parses the GSE154918 bulk RNA-seq raw count matrix and metadata from GEO, filter lowly expressed genes and exports clean CSV files.
"""
import gzip
import io
import os
import pandas as pd
import requests
import time
from pathlib import Path
from requests.exceptions import ChunkedEncodingError, ConnectionError

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = (PROJECT_ROOT/"data")

COUNTS_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154918/"
    "suppl/GSE154918_Schughart_Sepsis_200320.txt.gz"
)
SERIES_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154918/"
    "matrix/GSE154918_series_matrix.txt.gz"
)


def download_gz_text(url: str, cache_filename: str) -> str | None:
    cache_path = os.path.join(DATA_DIR, cache_filename)
    if os.path.exists(cache_path):
        print(f"Loading cached {cache_path} ...")
        with gzip.open(cache_path, "rt", encoding="utf-8", errors="replace") as f:
            return f.read()

    print(f"Downloading {url} ...")
    headers = {
        "User-Agent": "Mozilla/5.0 (X11; Linux x86_64)",
        "Accept-Encoding": "identity",
    }

    for attempt in range(1, 6):
        try:
            response = requests.get(url, headers=headers, timeout=60)
            response.raise_for_status()
            with open(cache_path, "wb") as out:
                out.write(response.content)
            with gzip.GzipFile(fileobj=io.BytesIO(response.content)) as f:
                return f.read().decode("utf-8", errors="replace")
        except (ChunkedEncodingError, ConnectionError) as err:
            if attempt == 5:
                raise
            wait = 2 ** attempt
            print(f"Stream interrupted ({err}), retrying in {wait}s...")
            time.sleep(wait)

def parse_metadata(metadata_text: str) -> pd.DataFrame:
    titles, accessions, status, sex = [], [], [], []
    for line in metadata_text.splitlines():
        if line.startswith("!Sample_title"):
            titles = line.replace('"', '').split("\t")[1:]
        elif line.startswith("!Sample_geo_accession"):
            accessions = line.replace('"', '').split("\t")[1:]
        elif line.startswith("!Sample_characteristics_ch1") and "status" in line.lower():
            status = line.replace('"', '').replace("status: ", "").split("\t")[1:]
        elif line.startswith("!Sample_characteristics_ch1") and "sex" in line.lower():
            sex = line.replace('"', '').replace("Sex: ", "").split("\t")[1:]

    metadata_clean = pd.DataFrame({
        "sample_title": titles,
        "accession": accessions,
        "raw_condition": status,
        "sex": sex
    })

    def simplify_status(status: str) -> str:
        status = status.lower()
        if "hlty" in status:
            return "healthy"
        elif "seps_p" in status or "shock_p" in status:
            return "sepsis"
        else:
            return "unknown"

    metadata_clean["condition"] = metadata_clean["raw_condition"].apply(simplify_status)
    metadata_clean = metadata_clean.drop(columns=["raw_condition"])
    metadata_clean = metadata_clean[metadata_clean["condition"].isin(["healthy", "sepsis"])].copy()
    return metadata_clean


def main() -> None:
    os.makedirs(DATA_DIR, exist_ok=True)

    metadata = parse_metadata(download_gz_text(SERIES_MATRIX_URL, "GSE154918_series_matrix.txt.gz"))

    counts_text = download_gz_text(COUNTS_URL, "GSE154918_counts.txt.gz")
    expr_raw = pd.read_csv(io.StringIO(counts_text), sep=r"\s+")

    expr_raw = expr_raw.dropna(subset=["gene_symbol"])
    expr_raw = expr_raw.set_index("gene_symbol")

    expr = expr_raw[metadata["sample_title"]].copy()
    expr = expr.groupby(expr.index).mean()
    expr = expr.loc[expr.std(axis=1)>0]

    assert list(expr.columns) == list(metadata["sample_title"]), "Sample alignment mismatch!"

    expr_path = os.path.join(DATA_DIR, "GSE154918_counts_clean.csv.gz")
    metadata_path = os.path.join(DATA_DIR, "GSE154918_metadata_clean.csv.gz")

    expr.to_csv(expr_path)
    metadata.to_csv(metadata_path, index=False)

    print(f"Exported {expr_path} {expr.shape} and {metadata_path} {metadata.shape}")


if __name__ == "__main__":
    main()