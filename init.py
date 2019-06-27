import pandas as pd

__DATASET_PATH__ = "/data2/datasets.tsv"
__DATASETS__ = pd.read_csv(__DATASET_PATH__, sep="\t").set_index("species")

__SPECIES__ = [
    {"label": "H.sapiens", "value": "hg38"},
    {"label": "M.musculus", "value": "mm10"},
    {"label": "C.albicans", "value": "SC5314"},
    {"label": "D.melanogaster", "value": "BDGP6"},
    {"label": "C.elegans", "value": "WBcel235"},
    {"label": "M.mulatta", "value": "Mmul8"},
    {"label": "G.gallus", "value": "GRCg6"},
    {"label": "R.norvegicus", "value": "Rnor6.0"},
    {"label": "P.troglodytes", "value": "panTro3"},
    {"label": "D.rerio", "value": "GRCz11"},
]
