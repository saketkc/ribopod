#!/usr/bin/env python
# coding: utf-8

# In[12]:


from collections import defaultdict

from pysradb import SRAdb
import os
import glob
import pandas as pd
from riboraptor.helpers import path_leaf, parse_star_logs, millify, order_dataframe
from riboraptor.cutadapt_to_json import cutadapt_to_json
from riboraptor.utils import summary_starlogs_over_runs, mkdir_p


ROOT_DIRS = [
    "/data1/re-ribo-analysis",
    "/data2/re-ribo-analysis/",
    "/data4/re-ribo-analysis",
]  # , "/data3/re-ribo-analysis", "/data4/re-ribo-analysis"]
ROOT_DIRS_SUMMARY = ["/data2/re-ribo-analysis-summary-tables"]
# This directory stores the summarized ORFs and their counts
# inside each build directory
ORF_TABLES_DIRNAME = "re-ribo-analysis-orf-tables"


def check_ribotricer_output_exists(srp, srx, assembly):
    for rootdir in ROOT_DIRS:
        path = os.path.join(
            rootdir,
            assembly,
            srp,
            "ribotricer_results",
            "{}_translating_ORFs.tsv".format(srx),
        )
        if os.path.exists(path):
            return path


def summarise_ribotricer_output_exists(path):
    df = pd.read_csv(path, sep="\t", use_cols=["ORF_ID"])
    return df


def check_ribotricer_metagene_exists(srp, srx, assembly):
    for rootdir in ROOT_DIRS:
        path_5p = os.path.join(
            rootdir,
            assembly,
            srp,
            "ribotricer_results",
            "{}_metagene_profiles_5p.tsv".format(srx),
        )
        path_3p = os.path.join(
            rootdir,
            assembly,
            srp,
            "ribotricer_results",
            "{}_metagene_profiles_3p.tsv".format(srx),
        )
        path_5p_tsv = None
        path_3p_tsv = None
        if os.path.exists(path_5p):
            path_5p_tsv = path_5p
        if os.path.exists(path_3p):
            path_3p_tsv = path_3p
        if os.path.exists(path_5p) or os.path.exists(path_3p):
            return path_5p_tsv, path_3p_tsv
    return None, None


def check_ribotricer_metagene_plot_exists(srp, srx, assembly):
    for rootdir in ROOT_DIRS:
        path = os.path.join(
            rootdir,
            assembly,
            srp,
            "ribotricer_results",
            "{}_metagene_plots.pdf".format(srx),
        )
        if os.path.exists(path):
            return path


def check_ribotricer_protocol_exists(srp, srx, assembly):
    for rootdir in ROOT_DIRS:
        path = os.path.join(
            rootdir, assembly, srp, "ribotricer_results", "{}_protocol.txt".format(srx)
        )
        if os.path.exists(path):
            return path


def check_ribotricer_bam_summary_exists(srp, srx, assembly):
    for rootdir in ROOT_DIRS:
        path = os.path.join(
            rootdir,
            assembly,
            srp,
            "ribotricer_results",
            "{}_bam_summary.txt".format(srx),
        )
        if os.path.exists(path):
            return path


def check_summarized_orfs_exists(srp, assembly):
    for rootdir in ROOT_DIRS_SUMMARY:
        path = os.path.join(
            rootdir, assembly, ORF_TABLES_DIRNAME, "{}_summarized_orfs.tsv".format(srp)
        )
        if os.path.exists(path):
            return path


def check_summarized_phase_scores_exists(srp, assembly):
    for rootdir in ROOT_DIRS_SUMMARY:
        path = os.path.join(
            rootdir,
            assembly,
            ORF_TABLES_DIRNAME,
            "{}_summarized_phase_scores.tsv".format(srp),
        )
        if os.path.exists(path):
            return path


# In[14]:


def create_df_from_dir(rootdir):
    """Create a dataframe struture amenable fro ribotricer for samples with no metadata using their directory

    Parameters
    ----------
    path: string
          Directory location
    """
    srp = path_leaf(rootdir)
    samples = glob.glob("{}/ribotricer_results/*_translating_ORFs.tsv".format(rootdir))
    samples = list(
        sorted(
            [
                path_leaf(sample).replace("_translating_ORFs.tsv", "")
                for sample in samples
            ]
        )
    )
    df = []
    for sample in samples:
        df.append((srp, sample, sample))
    df = pd.DataFrame(df)
    # print(rootdir, df)
    df.columns = ["study_accession", "experiment_accession", "run_accession"]
    df["library_layout"] = "SINGLE"
    df["bases"] = ""
    df["spots"] = ""
    df["avg_read_length"] = ""
    df["library_source"] = ""
    df["library_selection"] = ""
    df["adapter_spec"] = ""
    df["library_strategy"] = ""
    df["library_name"] = ""
    df["experiment_title"] = ""
    df["taxon_id"] = ""
    return df


# In[15]:


def get_srp_table(srp, assembly, re_ribo_analysis_dir):
    sradb = SRAdb("/data2/SRAmetadb.sqlite")
    column_order = [
        "study_accession",
        "experiment_title",
        "experiment_accession",
        "run_accession",
        "taxon_id",
        "library_selection",
        "library_layout",
        "library_strategy",
        "library_source",
        "library_name",
        "adapter_spec",
        "bases",
        "spots",
        "avg_read_length",
        "pass1_adapter",
        "pass1_total_reads_processed",
        "pass1_reads_with_adapters",
        "pass2_adapter",
        "pass2_total_reads_processed",
        "pass2_reads_with_adapters",
        "mapping_total_reads_input",
        "uniquely_mapped",
        "uniquely_mapped_percent",
        "ribotricer_orfs",
    ]
    filepath = os.path.join(re_ribo_analysis_dir, assembly, srp)
    if os.path.exists(filepath):

        try:
            srp_df = sradb.sra_metadata(
                srp.split("_")[0], detailed=True
            )  # , expand_sample_attributes=True)
        except:
            if "Kadosh" in filepath and "Kadosh_30C_37C" not in filepath:
                srp_df = pd.read_csv(
                    "/data2/Kadosh_design_files/{}.tsv".format(srp), sep="\t"
                )
            else:
                srp_df = create_df_from_dir(filepath)

            # return pd.DataFrame()
        srp_df.library_layout = srp_df.library_layout.fillna("SINGLE")
        srp_df = srp_df[srp_df.library_layout.str.contains("SINGLE")]

        srp_df["pass1_reads_with_adapters"] = None
        srp_df["pass1_total_reads_processed"] = None
        srp_df["pass1_adapter"] = None
        srp_df["pass2_adapter"] = None
        srp_df["pass2_total_reads_processed"] = None
        srp_df["pass2_reads_with_adapters"] = None
        srp_df["mapping_total_reads_input"] = None
        srp_df["uniquely_mapped"] = None
        srp_df["uniquely_mapped_percent"] = None
        srp_df["ribotricer_orfs"] = None
        srp_df["ribotricer_metagene_5p"] = None
        srp_df["ribotricer_metagene_3p"] = None

        srp_df["ribotricer_metagene_plot"] = None
        srp_df["ribotricer_protocol"] = None
        srp_df["ribotricer_bam_summary"] = None
        # srp_df["summarized_orfs"] = None
        # srp_df["summarized_phase_scores"] = None

        srpdir = os.path.join(re_ribo_analysis_dir, assembly, srp)
        starlogsdir = os.path.join(srpdir, "starlogs")
        srp_srx_grouped = srp_df.groupby("experiment_accession")
        preprocess_step1_dir = os.path.join(srpdir, "preprocessed_step1")
        preprocess_step2_dir = os.path.join(srpdir, "preprocessed")
        for srx, srx_group in srp_srx_grouped:
            ribotricer_output = check_ribotricer_output_exists(srp, srx, assembly)
            ribotricer_metagene_5p, ribotricer_metagene_3p = check_ribotricer_metagene_exists(
                srp, srx, assembly
            )

            ribotricer_bam_summary = check_ribotricer_bam_summary_exists(
                srp, srx, assembly
            )
            ribotricer_protocol = check_ribotricer_protocol_exists(srp, srx, assembly)
            ribotricer_metagene_plot = check_ribotricer_metagene_plot_exists(
                srp, srx, assembly
            )

            # summarized_orfs = check_summarized_orfs_exists(srp, assembly)
            # summarized_phase_score = check_summarized_orfs_exists(srp, assembly)

            srrs = srx_group["run_accession"].tolist()
            if ribotricer_output:
                srp_df.loc[
                    srp_df.experiment_accession == srx, "ribotricer_orfs"
                ] = ribotricer_output

            srp_df.loc[
                srp_df.experiment_accession == srx, "ribotricer_metagene_5p"
            ] = ribotricer_metagene_5p
            srp_df.loc[
                srp_df.experiment_accession == srx, "ribotricer_metagene_3p"
            ] = ribotricer_metagene_3p

            srp_df.loc[
                srp_df.experiment_accession == srx, "ribotricer_bam_summary"
            ] = ribotricer_bam_summary
            srp_df.loc[
                srp_df.experiment_accession == srx, "ribotricer_protocol"
            ] = ribotricer_protocol
            srp_df.loc[
                srp_df.experiment_accession == srx, "ribotricer_metagene_plot"
            ] = ribotricer_metagene_plot
            # srp_df.loc[srp_df.experiment_accession == srx, "summarized_orfs"] = summarized_orfs
            # srp_df.loc[srp_df.experiment_accession == srx, "summarized_phase_scores"] = summarized_phase_score

            # starlogs_df = summary_starlogs_over_runs(starlogsdir, srrs)

            for srr in srrs:
                starlogs_df = None
                if os.path.isfile(os.path.join(starlogsdir, srr + "Log.final.out")):
                    starlogs_df = parse_star_logs(
                        os.path.join(starlogsdir, srr + "Log.final.out")
                    )
                # Preprocessed_step1 adapter info
                step1_txt = os.path.join(
                    preprocess_step1_dir, srr + ".fastq.gz_trimming_report.txt"
                )
                step2_txt = os.path.join(
                    preprocess_step2_dir, srr + "_trimmed.fq.gz_trimming_report.txt"
                )
                step1_cutadapt_json = None
                step2_cutadapt_json = None

                if os.path.isfile(step1_txt):
                    step1_cutadapt_json = cutadapt_to_json(step1_txt)

                if os.path.isfile(step2_txt):
                    step2_cutadapt_json = cutadapt_to_json(step2_txt)

                if step1_cutadapt_json:
                    adapters = step1_cutadapt_json["adapters"]
                    if len(step1_cutadapt_json["adapters"]) == 0:
                        srp_df.loc[
                            srp_df.run_accession == srr, "pass1_adapter"
                        ] = "Empty?"
                    elif isinstance(adapters, str):
                        srp_df.loc[
                            srp_df.run_accession == srr, "pass1_adapter"
                        ] = step1_cutadapt_json["adapters"]
                    else:
                        srp_df.loc[
                            srp_df.run_accession == srr, "pass1_adapter"
                        ] = step1_cutadapt_json["adapters"][
                            "{} - {}".format(srr, "Adapter 1")
                        ]
                        trim_info1 = step1_cutadapt_json["trim_info"][srr]
                        srp_df.loc[
                            srp_df.run_accession == srr, "pass1_total_reads_processed"
                        ] = trim_info1["r_processed"]
                        srp_df.loc[
                            srp_df.run_accession == srr, "pass1_reads_with_adapters"
                        ] = trim_info1["r_with_adapters"]
                if step2_cutadapt_json:
                    adapters = step2_cutadapt_json["adapters"]
                    if len(step2_cutadapt_json["adapters"]) == 0:
                        srp_df.loc[
                            srp_df.run_accession == srr, "pass2_adapter"
                        ] = "Empty?"
                    elif isinstance(adapters, str):
                        srp_df.loc[
                            srp_df.run_accession == srr, "pass2_adapter"
                        ] = step2_cutadapt_json["adapters"]
                    else:
                        srp_df.loc[
                            srp_df.run_accession == srr, "pass2_adapter"
                        ] = step2_cutadapt_json["adapters"][
                            "{} - {}".format(srr + "_trimmed", "Adapter 1")
                        ]
                        trim_info2 = step2_cutadapt_json["trim_info"][srr + "_trimmed"]
                        srp_df.loc[
                            srp_df.run_accession == srr, "pass2_reads_with_adapters"
                        ] = trim_info2["r_with_adapters"]
                        srp_df.loc[
                            srp_df.run_accession == srr, "pass2_total_reads_processed"
                        ] = trim_info2["r_processed"]

                if starlogs_df:
                    srp_df.loc[
                        srp_df.run_accession == srr, "mapping_total_reads_input"
                    ] = starlogs_df["total_reads"]
                    srp_df.loc[
                        srp_df.run_accession == srr, "uniquely_mapped"
                    ] = starlogs_df["uniquely_mapped"]
                    srp_df.loc[
                        srp_df.run_accession == srr, "uniquely_mapped_percent"
                    ] = starlogs_df["uniquely_mapped_percent"]

        cols = [
            "bases",
            "spots",
            "pass1_reads_with_adapters",
            "pass2_reads_with_adapters",
            "pass2_total_reads_processed",
            "pass1_total_reads_processed",
            "uniquely_mapped",
            "mapping_total_reads_input",
        ]
        for col in cols:
            try:
                srp_df[col] = srp_df[col].apply(lambda z: millify(z))
            except:
                pass
        sradb.close()
        return order_dataframe(srp_df, column_order)


# In[16]:


READ_LENGTH_DIRNAME = "read_lengths"
METAGENE_COVERAGE_DIRNAME = "metagene_coverages"
METAGENE_LENWISE_COVERAGE_DIRNAME = "metagene_coverage_lengthwise"

# Top level directory of the directories inside each of the ROOT_DIRS
__ASSEMBLIES__ = [os.listdir(dirname) for dirname in ROOT_DIRS]
__SPECIES__ = [
    {"label": "H.sapiens", "value": "hg38"},
    {"label": "M.musculus", "value": "mm10"},
    {"label": "C.albicans", "value": "SC5314"},
]
__ASSEMBLIES__ = list(
    sorted(set([item for sublist in __ASSEMBLIES__ for item in sublist]))
)
__ASSEMBLY_WISE_SRP__ = defaultdict(list)
__SRP_TO_ROOT_DIR_MAP__ = defaultdict(dict)

# DATASETS = {"hg38": pd.read_csv("/data1/hg_datasets.tsv", sep="\t"),
#            "mm10": pd.read_csv("/data1/mm_datasets.tsv", sep="\t")}

for root_dir in ROOT_DIRS:
    for assembly_build in os.listdir(root_dir):
        for srp_dir in filter(
            os.path.isdir, glob.glob(os.path.join(root_dir, assembly_build, "*"))
        ):
            srp = os.path.basename(srp_dir)
            __ASSEMBLY_WISE_SRP__[assembly_build].append(srp)
            __SRP_TO_ROOT_DIR_MAP__[srp][assembly_build] = os.path.join(
                root_dir, assembly_build, srp
            )


__ASSEMBLY_WISE_SRP__ = defaultdict(list)
__SRP_TO_ROOT_DIR_MAP__ = defaultdict(dict)
for root_dir in ROOT_DIRS:
    for assembly_build in os.listdir(root_dir):
        for srp_dir in filter(
            os.path.isdir, glob.glob(os.path.join(root_dir, assembly_build, "*"))
        ):
            srp = os.path.basename(srp_dir)
            __ASSEMBLY_WISE_SRP__[assembly_build].append(srp)
            __SRP_TO_ROOT_DIR_MAP__[srp][assembly_build] = os.path.join(
                root_dir, assembly_build, srp
            )


# In[18]:


__SRP_TO_ROOT_DIR_MAP__


# In[19]:


__ASSEMBLY_WISE_SRP__


# In[20]:


def get_fragment_lengths(file_path):
    try:
        return pd.read_csv(file_path, sep="\t").fragment_length.tolist()
    except:
        # Handle 3 headed files
        df = pd.read_csv(file_path, header=None, sep="\t")
        df.columns = ["fragment_length", "offset_5p", "profile"]
        return df.fragment_length.tolist()


# In[21]:


db = SRAdb("/data2/SRAmetadb.sqlite")
all_projects = []


for species, sample_list in __ASSEMBLY_WISE_SRP__.items():
    mkdir_p("/data2/re-ribo-analysis-metadata/{}".format(species))
    for srp in sample_list:
        basedir = os.path.dirname(
            os.path.dirname(__SRP_TO_ROOT_DIR_MAP__[srp][species])
        )
        if not os.listdir(__SRP_TO_ROOT_DIR_MAP__[srp][species]):
            continue
        print(srp, basedir)
        df = get_srp_table(srp, species, basedir)
        project_filepath = "{}/{}/{}".format(basedir, species, srp)
        metadata_filepath = "/data2/re-ribo-analysis-metadata/{}/{}.tsv".format(
            species, srp
        )
        df_subset = df[
            df.ribotricer_metagene_5p == df.ribotricer_metagene_5p
        ].ribotricer_metagene_5p.tolist()
        summarized_orfs = check_summarized_orfs_exists(srp, species)
        summarized_phase_score = check_summarized_phase_scores_exists(srp, species)
        fragment_lengths = []
        for f in df_subset:
            fragment_lengths += get_fragment_lengths(f)
        fragment_lengths = list(sorted(list(set(fragment_lengths))))
        all_projects.append(
            (
                species,
                srp,
                project_filepath,
                metadata_filepath,
                str(fragment_lengths),
                summarized_orfs,
                summarized_phase_score,
            )
        )
        df.to_csv(metadata_filepath, sep="\t", index=False, header=True)


# In[22]:


summary_df = pd.DataFrame(all_projects)
summary_df.columns = [
    "species",
    "srp",
    "project_output_path",
    "project_metadata_path",
    "fragment_lengths",
    "summarized_orfs",
    "summarized_phase_scores",
]
summary_df = summary_df.sort_values(by=["species", "srp"])
summary_df.to_csv("/data2/datasets.tsv", sep="\t", index=False, header=True)
