"""
Downloads and pases the GSE154918 bulk RNA-seq raw count matrix and metadata from GEO, filter lowly expressed genes and exports clean CSV files.
"""
import gzip
import io

import pandas as pd
import requests

DATA_DIR ="data"
COUNTS_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154918/"
    "suppl/GSE154918_Schughart_Sepsis_200320.txt.gz"
)
SERIES_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154918/"
    "matrix/GSE154918_series_matrix.txt.gz"
)


def download_gz_text(url: str) -> str:
    print(f"Downloading {url} ...")
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    with gzip.GzipFile(fileobj=io.BytesIO(response.content)) as f:
        return f.read().decode("utf-8", errors="replace")

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


