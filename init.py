import pandas as pd

__DATASET_PATH__ = "/data2/datasets.tsv"
__DATASETS__ = pd.read_csv(__DATASET_PATH__, sep="\t").set_index("species")
__SPECIES__ = sorted(
    [
        {"label": "H.sapiens (human)", "value": "hg38"},
        {"label": "M.musculus (mouse)", "value": "mm10"},
        {"label": "C.albicans (candida albicans)", "value": "SC5314"},
        {"label": "D.melanogaster (fruitfly)", "value": "BDGP6"},
        {"label": "C.elegans (caenorhabditis elegans)", "value": "WBcel235"},
        {"label": "M.mulatta (macaque)", "value": "Mmul8"},
        {"label": "G.gallus (chicken)", "value": "GRCg6"},
        {"label": "R.norvegicus (rat)", "value": "Rnor6.0"},
        {"label": "P.troglodytes (chimp)", "value": "panTro3"},
        {"label": "D.rerio (zebrafish)", "value": "GRCz11"},
        {"label": "S.cerevisiae (yeast)", "value": "R64-1-1"},
        {"label": "A.thaliana (arabidopsis)", "value": "TAIR10"},
        {"label": "S.pombe (pombe)", "value": "ASM294v2"},
        {"label": "E.coli (ecoli)", "value": "ASM584v2"},
        {"label": "P.falciparum (plasmodium falciparum)", "value": "Epr1"},
        # {"label": "C.griseus (hamster)", "value": "CriGri_1.0"},
        {"label": "T.brucei (african trypanosome)", "value": "TryBruApr2005chr11"},
        {"label": "X.laevis (frog)", "value": "JGI_4.2"},
    ],
    key=lambda species: species["label"],
)


def build_to_species(build):
    for element in __SPECIES__:
        if element["value"] == build:
            return element["label"]
