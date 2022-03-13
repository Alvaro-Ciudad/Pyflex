import requests


def alphafold_downloader(pdb_id):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{pdb_id}-F1-model_v2.pdb"
    r = requests.get(url, allow_redirects=True).content.decode("utf-8")
    with open(f"{pdb_id}.pdb", "w") as f:
        f.write(r)
    return r
