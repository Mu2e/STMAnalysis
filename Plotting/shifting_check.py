"""
Plot particle distributions at shifting stage.
"""

# TODOs
# Add the stats boxes to the plots
# Fix the legend style

import array
import glob
import os

import ROOT
import pandas as pd
import uproot
import numpy as np

INPUT_DIRS = {
    "Ele": "InputDatasetBiasing/output_data/S1/Ele/StepPointMCs/",
    "Mu": "InputDatasetBiasing/output_data/S1/Mu/StepPointMCs/",
}

TBRANCH_NAMES = ["virtualdetectorId", "x", "y", "Ekin"]

# Plot settings
N_BINS = 100
INDIVIDUAL_PLOT_PDGIDS = [11, -11, 13, -13, 22, 2112]

# Data handling variables
TREE_NAME = "ROOTDump/ttree"
ROOT.gROOT.SetBatch(True)  # pylint: disable=no-member
ROOT.gStyle.SetOptStat(0)  # pylint: disable=no-member
cLow = ROOT.TCanvas("cLow", "", 800, 600)  # pylint: disable=no-member
cHigh = ROOT.TCanvas("cHigh", "", 1600, 1200)  # pylint: disable=no-member

# Boring stuff to make the code more legible
BOOL_VALUES = [True, False]

# The S1 simulation of EleBeamCat and MuBeamCat have effective POT counts, define the ratio of these:
MU_ELE_RATIO = 9.95e12 / 2.37e11

# Define the titles for the TH1D plots
TH1D_PLOT_TITLES = {
    "VD101": "Energy at VD101",
    "HPGe Absorber": "Energy at HPGe Absorber Region",
    "LaBr Absorber": "Energy at LaBr Absorber Region",
    "HPGe SSC": "Energy at HPGe SSC Aperture Region",
    "LaBr SSC": "Energy at LaBr SSC Aperture Region",
}

# Define the legend entries and titles for overlap plots
OVERLAP_PLOT_LEGEND_ENTRIES_AND_TITLE = {
    "HPGe SSC vs LaBr SSC": {
        "legend": {
            "Ele": "HPGe in SSC Aperture",
            "LaBr": "LaBr in SSC Aperture",
        },
        "title": "SSC Aperture energy distributions",
    },
    "HPGe Absorber vs LaBr Absorber": {
        "legend": {
            "Ele": "HPGe in Absorber Aperture",
            "LaBr": "LaBr in Absorber Aperture",
        },
        "title": "Absorber Aperture energy distributions",
    },
    "VD101 vs HPGe Absorber": {
        "legend": {
            "Ele": "VD101 Data",
            "LaBr": "HPGe in Absorber Aperture",
        },
        "title": "VD101 vs HPGe Absorber energy distributions",
    },
    "HPGe Absorber vs HPGe SSC": {
        "legend": {
            "Ele": "HPGe in Absorber Aperture",
            "LaBr": "HPGe in SSC Aperture",
        },
        "title": "HPGe Absorber vs SSC energy distributions",
    },
    "VD101 vs LaBr Absorber": {
        "legend": {
            "Ele": "VD101 Data",
            "LaBr": "LaBr in Absorber Aperture",
        },
        "title": "VD101 vs LaBr Absorber energy distributions",
    },
    "LaBr Absorber vs SSC": {
        "legend": {
            "Ele": "LaBr in Absorber Aperture",
            "LaBr": "LaBr in SSC Aperture",
        },
        "title": "LaBr Absorber vs SSC energy distributions",
    },
    "VD101 vs HPGe SSC": {
        "legend": {
            "Ele": "VD101 Data",
            "LaBr": "HPGe in SSC Aperture",
        },
        "title": "VD101 vs HPGe SSC energy distributions",
    },
    "VD101 vs LaBr SSC": {
        "legend": {
            "Ele": "VD101 Data",
            "LaBr": "LaBr in SSC Aperture",
        },
        "title": "VD101 vs LaBr SSC energy distributions",
    },
}


def parse_dir(data_dir: str) -> pd.DataFrame:
    """
    Parse all ROOT files in the specified directory and return a DataFrame
    of particles going through FILTER_VDID.
    """
    all_dfs = []
    files = sorted(glob.glob(f"{data_dir}/*.root"))
    num_files = len(files)
    if num_files == 0:
        raise ValueError(f"No ROOT files found in directory: {data_dir}")
    counter = 1
    for input_file in files:
        print(f"Reading file {counter}/{num_files}:", input_file)
        counter += 1
        with uproot.open(input_file) as f:
            arrs = f[TREE_NAME].arrays(TBRANCH_NAMES, library="np")
            if arrs is None or any(len(arrs[param]) == 0 for param in TBRANCH_NAMES):
                print(
                    "Empty arrs:",
                    [param for param in TBRANCH_NAMES if len(arrs[param]) == 0],
                )
                raise ValueError(
                    f"No entries found in table {TREE_NAME} in file: " f"{input_file}"
                )
            save_arrs = {}
            for param in TBRANCH_NAMES:
                save_arrs[param] = arrs[param]
            df_tmp = pd.DataFrame(save_arrs).copy()
            if df_tmp.empty:
                raise ValueError(
                    f"DataFrame is empty after loading table {TREE_NAME} in "
                    f"file: {input_file}"
                )
            df_tmp = df_tmp[df_tmp["virtualdetectorId"] == 101]
            all_dfs.append(df_tmp)
            # break  # For debugging, remove later
    return pd.concat(all_dfs, ignore_index=True)


def read_data_from_dir() -> pd.DataFrame:
    """Read data from all ROOT files in the specified directory."""

    ele_dir = INPUT_DIRS["Ele"]
    print("Reading ele files from directory:", ele_dir)
    read_ele_df = parse_dir(ele_dir)
    if read_ele_df.empty:
        raise ValueError(f"No data found in directory: {ele_dir}")

    mu_dir = INPUT_DIRS["Mu"]
    print("Reading mu files from directory:", mu_dir)
    read_mu_df = parse_dir(mu_dir)
    if read_mu_df.empty:
        raise ValueError(f"No data found in directory: {mu_dir}")

    return read_ele_df, read_mu_df


def convert_plot_title_to_file_name(title_key: str) -> str:
    """Convert plot title to a valid file name."""
    return title_key.replace(" ", "_").replace("-", "_").lower()


def plot_th1d(
    ele_data: pd.DataFrame,
    mu_data: pd.DataFrame,
    title_key: str,
    eMin: float = 0.0,
    eMax: float = 2.0,
) -> None:
    """Plot TH1D histograms for electron and muon data."""

    # Get energy values
    energies_ele = ele_data["Ekin"].values
    energies_mu = mu_data["Ekin"].values

    # Check if there are any energies to plot
    if len(energies_ele) + len(energies_mu) == 0:
        print(f"Warning: No energy data available for plot {title_key}, skipping.")
        return

    # Define the bin width
    bin_wid = (eMax - eMin) / N_BINS

    # Create histograms
    h_ele = ROOT.TH1D(  # pylint: disable=no-member
        "hEle", TH1D_PLOT_TITLES[title_key], N_BINS, eMin, eMax
    )
    for energy in energies_ele:
        h_ele.Fill(energy)
    h_mu = ROOT.TH1D(  # pylint: disable=no-member
        "hMu", TH1D_PLOT_TITLES[title_key], N_BINS, eMin, eMax
    )
    for energy in energies_mu:
        h_mu.Fill(energy)

    # Make these plots reflect the same number of POTs
    h_ele.Scale(MU_ELE_RATIO)

    # Check if there is data to merge
    if h_ele.GetEntries() + h_mu.GetEntries() == 0:
        print(f"Warning: Both histograms for plot {title_key} are empty, skipping.")
        return

    # Define the file name
    file_name = convert_plot_title_to_file_name(TH1D_PLOT_TITLES[title_key])

    # Create a merged histogram
    h_total = h_ele.Clone("hTotal")
    h_total.Add(h_mu)

    # Define the plot title and axis labels
    main_title = TH1D_PLOT_TITLES[title_key]
    x_title = "Kinetic Energy (MeV)"
    y_title = f"Counts / {bin_wid:.2f} MeV"
    h_total.SetTitle(f"{main_title};{x_title};{y_title}")

    # Plot and save the low res canvas
    cLow.cd()
    h_total.Draw("HIST")
    cLow.SaveAs(f"{file_name}_th1d_low.png")

    # Plot and save the high res canvas
    cHigh.cd()
    h_total.Draw("HIST")
    cHigh.SaveAs(f"{file_name}_th1d_high.png")


def plot_th1d_overlap(
    ele_hpge_data: pd.DataFrame,
    mu_hpge_data: pd.DataFrame,
    ele_labr_data: pd.DataFrame,
    mu_labr_data: pd.DataFrame,
    title_key: str,
    eMin: float = 0.0,
    eMax: float = 2.0,
) -> None:
    """Plot TH1D histograms for electron and muon data."""

    # Get energy values
    energies_ele_hpge = ele_hpge_data["Ekin"].values
    energies_mu_hpge = mu_hpge_data["Ekin"].values
    energies_ele_labr = ele_labr_data["Ekin"].values
    energies_mu_labr = mu_labr_data["Ekin"].values

    # Check if there are any energies to plot
    if (
        len(energies_ele_hpge)
        + len(energies_mu_hpge)
        + len(energies_ele_labr)
        + len(energies_mu_labr)
        == 0
    ):
        print(f"Warning: No energy data available for plot {title_key}, skipping.")
        return

    # Define the bin width
    bin_wid = (eMax - eMin) / N_BINS

    # Create histograms for HPGe data
    h_ele_hpge = ROOT.TH1D(  # pylint: disable=no-member
        "hEleHPGe",
        OVERLAP_PLOT_LEGEND_ENTRIES_AND_TITLE[title_key]["title"],
        N_BINS,
        eMin,
        eMax,
    )
    for energy in energies_ele_hpge:
        h_ele_hpge.Fill(energy)
    h_mu_hpge = ROOT.TH1D(  # pylint: disable=no-member
        "hMuHPGe",
        OVERLAP_PLOT_LEGEND_ENTRIES_AND_TITLE[title_key]["title"],
        N_BINS,
        eMin,
        eMax,
    )
    for energy in energies_mu_hpge:
        h_mu_hpge.Fill(energy)

    # Create histograms for LaBr data
    h_ele_labr = ROOT.TH1D(  # pylint: disable=no-member
        "hEleLabr",
        OVERLAP_PLOT_LEGEND_ENTRIES_AND_TITLE[title_key]["title"],
        N_BINS,
        eMin,
        eMax,
    )
    for energy in energies_ele_labr:
        h_ele_labr.Fill(energy)
    h_mu_labr = ROOT.TH1D(  # pylint: disable=no-member
        "hMuLabr",
        OVERLAP_PLOT_LEGEND_ENTRIES_AND_TITLE[title_key]["title"],
        N_BINS,
        eMin,
        eMax,
    )
    for energy in energies_mu_labr:
        h_mu_labr.Fill(energy)

    # Make these plots reflect the same number of POTs
    h_ele_hpge.Scale(MU_ELE_RATIO)
    h_ele_labr.Scale(MU_ELE_RATIO)

    # Check if there is data to merge
    if (
        h_ele_hpge.GetEntries()
        + h_mu_hpge.GetEntries()
        + h_ele_labr.GetEntries()
        + h_mu_labr.GetEntries()
        == 0
    ):
        print(f"Warning: Both histograms for plot {title_key} are empty, skipping.")
        return

    # Define the file name
    file_name = convert_plot_title_to_file_name(
        OVERLAP_PLOT_LEGEND_ENTRIES_AND_TITLE[title_key]["title"]
    )

    # Create a merged histograms
    h_total_hpge = h_ele_hpge.Clone("hTotalHPGe")
    h_total_hpge.Add(h_mu_hpge)
    h_total_labr = h_ele_labr.Clone("hTotalLabr")
    h_total_labr.Add(h_mu_labr)
    h_total_hpge.SetLineColor(ROOT.kRed)  # pylint: disable=no-member
    h_total_labr.SetLineColor(ROOT.kBlue)  # pylint: disable=no-member
    h_total_hpge.SetLineWidth(2)
    h_total_labr.SetLineWidth(2)

    # Define the THStack
    hs = ROOT.THStack("hs", "")
    hs.Add(h_total_hpge)
    hs.Add(h_total_labr)

    # Define the legend
    legend = ROOT.TLegend(0.65, 0.7, 0.85, 0.85)  # pylint: disable=no-member
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(
        h_total_hpge,
        OVERLAP_PLOT_LEGEND_ENTRIES_AND_TITLE[title_key]["legend"]["Ele"],
        "l",
    )
    legend.AddEntry(
        h_total_labr,
        OVERLAP_PLOT_LEGEND_ENTRIES_AND_TITLE[title_key]["legend"]["LaBr"],
        "l",
    )

    # Define the plot title and axis labels
    main_title = OVERLAP_PLOT_LEGEND_ENTRIES_AND_TITLE[title_key]["title"]
    x_title = "Kinetic Energy (MeV)"
    y_title = f"Counts / {bin_wid:.2f} MeV"
    hs.SetTitle(f"{main_title};{x_title};{y_title}")

    # Plot and save the low res canvas
    cLow.cd()
    hs.Draw("nostack HIST")
    legend.Draw()
    cLow.SaveAs(f"{file_name}_th1d_low.png")

    # Plot and save the high res canvas
    cHigh.cd()
    hs.Draw("nostack HIST")
    legend.Draw()
    cHigh.SaveAs(f"{file_name}_th1d_high.png")

    return


ele_df, mu_df = read_data_from_dir()

# Filter data within ±25mm of the aperture center
hpge_ssc_aperture_x = -3945.0
labr_ssc_aperture_x = -3865.0
aperture_tolerance = 25.0

ele_hpge_absorber_df = ele_df[
    (ele_df["x"] >= hpge_ssc_aperture_x - aperture_tolerance)
    & (ele_df["x"] <= hpge_ssc_aperture_x + aperture_tolerance)
    & (ele_df["y"] >= -aperture_tolerance)
    & (ele_df["y"] <= aperture_tolerance)
]
ele_labr_absorber_df = ele_df[
    (ele_df["x"] >= labr_ssc_aperture_x - aperture_tolerance)
    & (ele_df["x"] <= labr_ssc_aperture_x + aperture_tolerance)
    & (ele_df["y"] >= -aperture_tolerance)
    & (ele_df["y"] <= aperture_tolerance)
]

mu_hpge_absorber_df = mu_df[
    (mu_df["x"] >= hpge_ssc_aperture_x - aperture_tolerance)
    & (mu_df["x"] <= hpge_ssc_aperture_x + aperture_tolerance)
    & (mu_df["y"] >= -aperture_tolerance)
    & (mu_df["y"] <= aperture_tolerance)
]
mu_labr_absorber_df = mu_df[
    (mu_df["x"] >= labr_ssc_aperture_x - aperture_tolerance)
    & (mu_df["x"] <= labr_ssc_aperture_x + aperture_tolerance)
    & (mu_df["y"] >= -aperture_tolerance)
    & (mu_df["y"] <= aperture_tolerance)
]


# Filter data within radius of aperture centers
aperture_radius = np.sqrt(50) / np.pi

ele_hpge_ssc_df = ele_df[
    (ele_df["x"] - hpge_ssc_aperture_x) ** 2 + ele_df["y"] ** 2 <= aperture_radius**2
]
ele_labr_ssc_df = ele_df[
    (ele_df["x"] - labr_ssc_aperture_x) ** 2 + ele_df["y"] ** 2 <= aperture_radius**2
]

mu_hpge_ssc_df = mu_df[
    (mu_df["x"] - hpge_ssc_aperture_x) ** 2 + mu_df["y"] ** 2 <= aperture_radius**2
]
mu_labr_ssc_df = mu_df[
    (mu_df["x"] - labr_ssc_aperture_x) ** 2 + mu_df["y"] ** 2 <= aperture_radius**2
]

# Generate energy distribution plots at VD101
# plot_th1d(
#     ele_df,
#     mu_df,
#     "VD101",
# )

# # Generate energy distribution plots at HPGe absorber
# plot_th1d(
#     ele_hpge_absorber_df,
#     mu_hpge_absorber_df,
#     "HPGe Absorber",
# )

# # Generate energy distribution plots at LaBr absorber
# plot_th1d(
#     ele_labr_absorber_df,
#     mu_labr_absorber_df,
#     "LaBr Absorber",
# )

# # Generate energy distribution plots at HPGe SSC aperture
# plot_th1d(
#     ele_hpge_ssc_df,
#     mu_hpge_ssc_df,
#     "HPGe SSC",
# )

# # Generate energy distribution plots at LaBr SSC aperture
# plot_th1d(
#     ele_labr_ssc_df,
#     mu_labr_ssc_df,
#     "LaBr SSC",
# )


# Generate overlap plots comparing HPGe and LaBr SSC data
plot_th1d_overlap(
    ele_hpge_ssc_df,
    mu_hpge_ssc_df,
    ele_labr_ssc_df,
    mu_labr_ssc_df,
    "HPGe SSC vs LaBr SSC",
)

# plot_th1d_overlap(
#     ele_hpge_absorber_df,
#     mu_hpge_absorber_df,
#     ele_labr_absorber_df,
#     mu_labr_absorber_df,
#     "HPGe Absorber vs LaBr Absorber",
# )

# plot_th1d_overlap(
#     ele_hpge_absorber_df,
#     mu_hpge_absorber_df,
#     ele_hpge_ssc_df,
#     mu_hpge_ssc_df,
#     "HPGe Absorber vs HPGe SSC",
# )

# plot_th1d_overlap(
#     ele_labr_absorber_df,
#     mu_labr_absorber_df,
#     ele_labr_ssc_df,
#     mu_labr_ssc_df,
#     "LaBr Absorber vs LaBr SSC",
# )

# plot_th1d_overlap(
#     ele_df, mu_df, ele_hpge_absorber_df, mu_hpge_absorber_df, "VD101 vs HPGe Absorber"
# )

# plot_th1d_overlap(
#     ele_df, mu_df, ele_labr_absorber_df, mu_labr_absorber_df, "VD101 vs LaBr Absorber"
# )

# plot_th1d_overlap(
#     ele_df,
#     mu_df,
#     ele_hpge_ssc_df,
#     mu_hpge_ssc_df,
#     "VD101 vs HPGe SSC",
# )

# plot_th1d_overlap(
#     ele_df,
#     mu_df,
#     ele_labr_ssc_df,
#     mu_labr_ssc_df,
#     "VD101 vs LaBr SSC",
# )
