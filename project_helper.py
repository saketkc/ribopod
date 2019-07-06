import pandas as pd


def get_projects(datasets, species):
    """Get projects for a particular species.

    Parameters
    ----------
    datasets: pandas.DataFrame
              datasets.tsv
    species: string
             Assembly

    Returns
    -------
    projects: list
              List of projects available for given species
    """
    projects = datasets.loc[datasets.index == species].srp.tolist()
    projects = [{"label": project, "value": project} for project in projects]
    return projects


def get_srp_table(datasets, srp):
    """Get metadata table for SRP.

    Parameters
    ----------
    datasets: pandas.DataFrame
              datasets.tsv
    srp: string
         SRP ID

    Returns
    -------
    table: dash.html_table
           Metadata table for SRP
    """
    dataset = datasets[datasets.srp == srp].iloc[0]
    srp_metadata = pd.read_csv(dataset.project_metadata_path, sep="\t")
    return srp_metadata


def get_srp_read_lengths(datasets, srp):
    """Get fragment lengths for a SRP.

    Parameters
    ----------
    datasets: pandas.DataFrame
    srp: string
         SRP ID

    Returns
    -------
    fragment_lengths: list(dict)
    """
    dataset = datasets[datasets.srp == srp].iloc[0]
    fragment_lengths = eval(dataset.fragment_lengths)
    return [{"label": length, "value": length} for length in fragment_lengths]


def get_summarized_phase_scores(datasets, srp):
    """Read summarised_phase_scores as dataframe.

    Parameters
    ----------
    datasets: pandas.DataFrame
    srp: string
         SRP ID

    Returns
    -------
    phase_scores_df: pd.DataFrame
    """
    dataset = datasets[datasets.srp == srp].iloc[0]
    phase_scores_df = pd.read_csv(dataset.summarized_phase_scores, sep="\t")
    print(phase_scores_df.head())
    phase_scores_df = phase_scores_df.set_index("ORF_ID")
    return phase_scores_df


def get_summarized_orf_counts(datasets, srp):
    """Read summarized_orfs as dataframe.

    Parameters
    ----------
    datasets: pandas.DataFrame
    srp: string
         SRP ID

    Returns
    -------
    orf_counts_df: pd.DataFrame
    """
    dataset = datasets[datasets.srp == srp].iloc[0]
    orf_counts_df = pd.read_csv(dataset.summarized_orfs, sep="\t").set_index(
        "experiment_accession"
    )
    return orf_counts_df


def get_project_summary_file(datasets, srp):
    """Get location of project's summary file.

    Parameters
    ----------
    datasets: pandas.DataFrame
    srp: string
         SRP ID

    Returns
    -------
    file_path: string
    """
    dataset = datasets[datasets.srp == srp].iloc[0]
    return dataset.project_metadata_path
