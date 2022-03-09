##################################################
## Author: Matt McGauley
## Credits: Emily Darin for ggplot code
## Email: mmcgau01@villanova.edu
##################################################

import os
import glob
import csv
import re
import shutil
import warnings
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotnine
from plotnine import *
from matplotlib import gridspec
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap
from matplotlib.backends.backend_pdf import PdfPages

warnings.filterwarnings("ignore")

plotPath = os.path.join(os.getcwd(), "Plots")
if not os.path.exists(plotPath):
    os.makedirs(plotPath)

dataPath = os.path.join(os.getcwd(), "Data")
if not os.path.exists(plotPath):
    os.makedirs(plotPath)

heatmapDataPath = os.path.join(dataPath, "HeatMapData")
if not os.path.exists(heatmapDataPath):
    os.makedirs(heatmapDataPath)

outputPath = os.path.join(dataPath, "OutputData")
if not os.path.exists(outputPath):
    os.makedirs(outputPath)

resampledPath = os.path.join(dataPath, "ResampledData")
if not os.path.exists(resampledPath):
    os.makedirs(resampledPath)

metadataPath = os.path.join(dataPath, "MetaData")
if not os.path.exists(metadataPath):
    os.makedirs(metadataPath)

heatmapPlotPath = os.path.join(plotPath, "HeatMapPlots")
if not os.path.exists(heatmapPlotPath):
    os.makedirs(heatmapPlotPath)

quartilePlotPath = os.path.join(plotPath, "QuartilePlots")
if not os.path.exists(quartilePlotPath):
    os.makedirs(quartilePlotPath)

scatterDataPath = os.path.join(dataPath, "ScatterData")
if not os.path.exists(scatterDataPath):
    os.makedirs(scatterDataPath)

scatterInflow = os.path.join(scatterDataPath, "Inflow")
if not os.path.exists(scatterInflow):
    os.makedirs(scatterInflow)

scatterOutflow = os.path.join(scatterDataPath, "Outflow")
if not os.path.exists(scatterOutflow):
    os.makedirs(scatterOutflow)

exportPath = os.path.join(os.getcwd(), "Export")
if not os.path.exists(exportPath):
    os.makedirs(exportPath)

all_merged_resample = []

def rpd(original, resampled):
    resampled = resampled.rename(columns = {"EMC" : "resampledEMC"})
    for_export = resampled.merge(original, how = "left")
    all_merged_resample.append(for_export[["Run", "N", "Location", "Analyte", "Quartile", "resampledEMC", "EMC"]])
    resampled = resampled.merge(original, how = "left", on = "Quartile").drop(columns = ["Location_x", "Analyte_x", "Flow_x", "Location_y", "Analyte_y", "Flow_y"])
    resampled["rpd"] = resampled.apply(lambda x: abs(x["resampledEMC"] - x["EMC"])/((x["resampledEMC"] + x["EMC"])/2), axis = 1)
    return resampled.groupby(["N", "Quartile"])["rpd"].agg("mean")

def create_output():
    all_original_files = []
    subDirs = glob.glob(os.path.join(outputPath, "*"))
    for dir in subDirs:
        csv_files = glob.glob(os.path.join(dir, "*.csv"))
        all_original_files.extend(csv_files)

    all_resampled_files = []
    resampledSubDirs = glob.glob(os.path.join(resampledPath, "*"))
    for dir in resampledSubDirs:
        csv_files = glob.glob(os.path.join(dir, "*.csv"))
        all_resampled_files.extend(csv_files)

    all_metadata_files = []
    metadataSubDirs = glob.glob(os.path.join(metadataPath, "*"))
    for dir in metadataSubDirs:
        csv_files = glob.glob(os.path.join(dir, "*.csv"))
        all_metadata_files.extend(csv_files)

    Copper_ScatterData = pd.DataFrame(columns = ["N", "Quartile", "median", "std", "rpd", "Location", "Flow", "Analyte", "NumObservations", "NearMedianDetectionLimit", "NumObservationsBelowDetect","changeRPD", "percentOriginal", "skewnessCategory", "skewnessValue", "coefficientOfVariation", "averageAnnualMonthsSampled", "inductionRatio", "distributionCategory"])
    TSS_ScatterData = pd.DataFrame(columns = ["N", "Quartile", "median", "std", "rpd", "Location", "Flow", "Analyte", "NumObservations", "NearMedianDetectionLimit", "NumObservationsBelowDetect","changeRPD", "percentOriginal", "skewnessCategory", "skewnessValue", "coefficientOfVariation", "averageAnnualMonthsSampled", "inductionRatio", "distributionCategory"])
    Phosphorus_ScatterData = pd.DataFrame(columns = ["N", "Quartile", "median", "std", "rpd", "Location", "Flow", "Analyte", "NumObservations", "NearMedianDetectionLimit", "NumObservationsBelowDetect", "changeRPD", "percentOriginal", "skewnessCategory", "skewnessValue", "coefficientOfVariation", "averageAnnualMonthsSampled", "inductionRatio",  "distributionCategory"])

    for dataFile in all_original_files:
        fileName = os.path.basename(dataFile).split(".")[0]
        resampledDataFiles = [resampledDataFile for resampledDataFile in all_resampled_files if fileName in resampledDataFile]
        metadataFiles = [metadataDataFile for metadataDataFile in all_metadata_files if fileName in metadataDataFile]

        originalDataframe = pd.read_csv(dataFile)
        resampledDataframes = [pd.read_csv(file) for file in resampledDataFiles]
        metadataDataframe = [pd.read_csv(file) for file in metadataFiles][0]

        for frame in resampledDataframes:
            frame["N"] =  "n=" + str(int(len(frame)/500))
            frame["N"] = frame["N"].astype("category")

        fullResampledDataframe = pd.concat(resampledDataframes, ignore_index = True)
        fullResampledDataframe["N"] = fullResampledDataframe["N"].astype("category")
        values = fullResampledDataframe["N"].unique()

        if "n=20" not in values and "n=25" not in values and "n=15" not in values:
            fullResampledDataframe["N"] = fullResampledDataframe["N"].cat.reorder_categories(["n=5", "n=10"])
            color_list = ["#2c925f", "#2e4057", "#000000"]
        if "n=20" not in values and "n=25" not in values and "n=15" in values:
            fullResampledDataframe["N"] = fullResampledDataframe["N"].cat.reorder_categories(["n=5", "n=10", "n=15"])
            color_list = ["#2c925f", "#2e4057", "#810f7c", "#000000"]
        if "n=20" in values and "n=25" not in values:
            fullResampledDataframe["N"] = fullResampledDataframe["N"].cat.reorder_categories(["n=5", "n=10", "n=15", "n=20"])
            color_list = ["#2c925f", "#2e4057", "#810f7c", "#de2d26", "#000000"]
        if "n=20" in values and "n=25" in values:
            fullResampledDataframe["N"] = fullResampledDataframe["N"].cat.reorder_categories(["n=5", "n=10", "n=15", "n=20", "n=25"])
            color_list = ["#2c925f", "#2e4057", "#810f7c", "#de2d26", "#41b6c4", "#000000"]

        Quartiles = [0.1, 0.25, 0.5, 0.75, 0.9]


        resampledQuartiles = fullResampledDataframe.groupby(["Run", "N", "Location", "Analyte", "Flow"])["EMC"].quantile(Quartiles).reset_index().rename(columns = {"level_5" : "Quartile"})
        originalQuartiles = originalDataframe.groupby(["Location", "Analyte", "Flow"])["EMC"].quantile(Quartiles).reset_index().rename(columns = {"level_3" : "Quartile"})
        originalQuartiles["resampled_vline_label"] = "Quartile\nEMC"
        originalQuartiles["original_vline_label"] = "Parent\ndataset\nmedian\ndetection\nlimit"

        heatmapData = resampledQuartiles.groupby(["N", "Quartile"])["EMC"].agg(["median", "std"])
        heatmapData["rpd"] = rpd(originalQuartiles, resampledQuartiles) * 100
        heatmapData["Location"] = originalQuartiles["Location"].iloc[0]
        heatmapData["Flow"] = originalQuartiles["Flow"].iloc[0]
        heatmapData["Analyte"] = originalQuartiles["Analyte"].iloc[0]
        heatmapData["NumObservations"] = metadataDataframe["Number_of_final_Observations"].iloc[0]
        heatmapData["NumObservationsBelowDetect"] = metadataDataframe["Number_of_below_detect_Observations"].iloc[0]
        heatmapData["NearMedianDetectionLimit"] = metadataDataframe["Number_of_Obs_near_dl"].iloc[0]
        heatmapData["skewnessCategory"] = metadataDataframe["skewnessCategory"].iloc[0]
        heatmapData["distributionCategory"] = metadataDataframe["distributionCategory"].iloc[0]
        heatmapData["skewnessValue"] = metadataDataframe["skewnessValue"].iloc[0]
        heatmapData["coefficientOfVariation"] = metadataDataframe["coefficientOfVariation"].iloc[0]
        heatmapData["averageAnnualMonthsSampled"] = metadataDataframe["averageAnnualMonthsSampled"].iloc[0]
        heatmapData["inductionRatio"] = metadataDataframe["inductionRatio"].iloc[0]
        heatmapData = heatmapData.reset_index()
        heatmapData["changeRPD"] = heatmapData.groupby(["Quartile"])["rpd"].apply(lambda row: abs(row - row.shift())).fillna(0)
        heatmapData["percentOriginal"] = heatmapData.apply(lambda row: (int(row["N"].split("=")[1])/row["NumObservations"]) * 100, axis = 1)

        originalQuartiles["mdl"] = metadataDataframe["mdl"].iloc[0]

        HeatMapDestination = os.path.join(heatmapDataPath, originalQuartiles["Location"].iloc[0])
        if not os.path.exists(HeatMapDestination):
            os.makedirs(HeatMapDestination)

        heatmapData.to_csv(os.path.join(HeatMapDestination, originalQuartiles["Location"].iloc[0] + "_" + originalQuartiles["Flow"].iloc[0] + "_" + originalQuartiles["Analyte"].iloc[0] + ".csv"))

        if originalQuartiles["Analyte"].iloc[0] == "Copper":
            Copper_ScatterData = Copper_ScatterData.append(heatmapData)
        if originalQuartiles["Analyte"].iloc[0] == "Phosphorus":
            Phosphorus_ScatterData = Phosphorus_ScatterData.append(heatmapData)
        if originalQuartiles["Analyte"].iloc[0] == "TSS":
            TSS_ScatterData = TSS_ScatterData.append(heatmapData)

        LocAnalyte = os.path.join(quartilePlotPath, originalQuartiles["Location"].iloc[0] + "_" + originalQuartiles["Analyte"].iloc[0])
        if not os.path.exists(LocAnalyte):
            os.makedirs(LocAnalyte)

        with PdfPages(os.path.join(LocAnalyte, originalQuartiles["Location"].iloc[0] + "_" + originalQuartiles["Flow"].iloc[0] + "_" + originalQuartiles["Analyte"].iloc[0] + ".pdf")) as pdf:
            for sample in resampledQuartiles["Quartile"].unique():
                if originalQuartiles["Analyte"].iloc[0] == "Copper":
                    xlabel = "EMC (µg/L)"
                else:
                    xlabel = "EMC (mg/L)"

                originalQuartiles["Analyte"].replace({
                        "Copper" : "Dissolved Copper",
                        "Phosphorus" : "Total Phosphorus",
                        "TSS" : "Total Suspended Solids"
                        }, inplace = True)

                finalPlot, ax1 = plt.subplots()
                left, bottom, width, height = [0.25, 0.6, 0.2, 0.2]
                ax2 = finalPlot.add_axes([left, bottom, width, height])

                thisPlot = (ggplot(resampledQuartiles[resampledQuartiles["Quartile"] == sample], plotnine.aes(x = "EMC", label = "Location", color = "N")) +
                                geom_line(stat = "density", size = 2) +
                                geom_vline(originalQuartiles[originalQuartiles["Quartile"] == sample], plotnine.aes(xintercept = "EMC", color = "resampled_vline_label"), linetype = "dashed", size = 2, show_legend = {"color" : True}) +
                                scale_color_manual(values = color_list) +
                                guides(color = guide_legend(title = "Number of\nMonitoring\nEvents\n\n\n")) +
                                ggtitle(originalQuartiles["Location"].iloc[0] + "\n" + originalQuartiles["Flow"].iloc[0] + "\n" + "Quartile " + str(int(sample * 100)) + "\n" + originalQuartiles["Analyte"].iloc[0])+
                                xlab(xlabel) +
                                expand_limits(x = 0) +
                                theme_light() +
                                theme(legend_position = (0.95, .75), axis_title_y = element_blank(), axis_text_y = element_blank(), text = element_text(size = 30, family = "serif"))
                            )


                scatter = (ggplot(originalDataframe, plotnine.aes(x = "EMC", y = originalDataframe.index)) +
                                geom_point() +
                                geom_vline(originalQuartiles[originalQuartiles["Quartile"] == sample], plotnine.aes(xintercept = "mdl", color = "original_vline_label"), linetype = "dashed", size = 2, show_legend = {"color" : True}) +
                                scale_color_manual(values = ["#FF0000"]) +
                                expand_limits(x = 0) +
                                guides(color = guide_legend(title = "")) +
                                theme_light() +
                                theme(axis_title_y = element_blank(), axis_text_y = element_blank(), text = element_text(size = 30, family = "serif"))
                            )

                fig = (ggplot() + theme_void() + theme(figure_size = (30, 20))).draw()

                gs = gridspec.GridSpec(40, 40)
                ax1 = fig.add_subplot(gs[:,:])

                _ = thisPlot._draw_using_figure(fig, [ax1])

                ax1.text(0.3 * (ax1.get_xlim()[1]), 1.01 * (ax1.get_ylim()[1]),
                "Number of raw observations: " + str(metadataDataframe["Number_of_raw_Observations"].iloc[0]) + "\n" +
                "Number of below-detection observations: " + str(metadataDataframe["Number_of_below_detect_Observations"].iloc[0]) + "\n" +
                "Number of excluded values: " + str(metadataDataframe["Number_of_excluded_Observations"].iloc[0]) + "\n" +
                "Number of final observations: " + str(metadataDataframe["Number_of_final_Observations"].iloc[0]),
                size = 30, fontfamily = "serif")

                ax1.set_title(originalQuartiles["Location"].iloc[0] + "\n" + originalQuartiles["Flow"].iloc[0] + "\n" + "Quartile " + str(int(sample * 100)) + "\n" + originalQuartiles["Analyte"].iloc[0], fontsize = 30, fontfamily = "serif", loc = "left")
                ax1.set_xlabel(xlabel, size = 30, fontfamily = "serif")

                if originalQuartiles[originalQuartiles["Quartile"] == sample]["EMC"].iloc[0]/ax1.get_xlim()[1] > 0.5:
                    ax2 = fig.add_subplot(gs[2:10:, 2:10])
                elif originalQuartiles[originalQuartiles["Quartile"] == sample]["EMC"].iloc[0]/ax1.get_xlim()[1] < 0.5:
                    ax2 = fig.add_subplot(gs[2:10, -10:-2])
                else:
                    ax2 = fig.add_subplot(gs[2:10:, 2:10])

                _ = scatter._draw_using_figure(fig, [ax2])


                ax2.set_title("Original Data Distribution", fontsize = 30, fontfamily = "serif")
                ax2.set_xlabel(xlabel, size = 30, fontfamily = "serif")

                pdf.savefig(fig)
            plt.close("all")

    Copper_ScatterData[Copper_ScatterData["Flow"] == "Inflow"].to_csv(os.path.join(scatterInflow, "Copper_inflow.csv"), index = False)
    TSS_ScatterData[TSS_ScatterData["Flow"] == "Inflow"].to_csv(os.path.join(scatterInflow, "TSS_inflow.csv"), index = False)
    Phosphorus_ScatterData[Phosphorus_ScatterData["Flow"] == "Inflow"].to_csv(os.path.join(scatterInflow, "Phosphorus_inflow.csv"), index = False)

    Copper_ScatterData[Copper_ScatterData["Flow"] == "Outflow"].to_csv(os.path.join(scatterOutflow, "Copper_outflow.csv"), index = False)
    TSS_ScatterData[TSS_ScatterData["Flow"] == "Outflow"].to_csv(os.path.join(scatterOutflow, "TSS_outflow.csv"), index = False)
    Phosphorus_ScatterData[Phosphorus_ScatterData["Flow"] == "Outflow"].to_csv(os.path.join(scatterOutflow, "Phosphorus_outflow.csv"), index = False)

    combined_ScatterData = pd.concat([Copper_ScatterData, TSS_ScatterData, Phosphorus_ScatterData])
    combined_ScatterData.to_csv(os.path.join(scatterDataPath, "all_scatter_data.csv"), index = False)

    merged_original_and_resampled_data = pd.concat(all_merged_resample)
    merged_original_and_resampled_data.to_csv(os.path.join(scatterDataPath, "merged_original_and_resampled_data.csv"), index = False)

def create_heatmap():
    subDirs = glob.glob(os.path.join(heatmapDataPath, "*"))

    for dir in subDirs:
        inflow = [file for file in glob.glob(os.path.join(dir, "*.csv")) if "inflow" in os.path.basename(file).lower()]
        inflow_dfs = [pd.read_csv(file) for file in inflow]
        if inflow_dfs:
            for df in inflow_dfs:
                values = df["N"].unique()
                if "n=20" not in values and "n=25" not in values and "n=15" not in values:
                    df["N"] = df["N"].astype("category").cat.reorder_categories(["n=10", "n=5"])
                if "n=20" not in values and "n=25" not in values and "n=15" in values:
                    df["N"] = df["N"].astype("category").cat.reorder_categories(["n=15", "n=10", "n=5"])
                if "n=20" in values and "n=25" not in values:
                    df["N"] = df["N"].astype("category").cat.reorder_categories(["n=20", "n=15", "n=10", "n=5"])
                if "n=20" in values and "n=25" in values:
                    df["N"] = df["N"].astype("category").cat.reorder_categories(["n=25", "n=20", "n=15", "n=10", "n=5"])

            inflow_combined = pd.concat(inflow_dfs)
            get_heatmap(inflow_combined)

        outflow = [file for file in glob.glob(os.path.join(dir, "*.csv")) if "outflow" in os.path.basename(file).lower()]
        outflow_dfs = [pd.read_csv(file) for file in outflow]
        if outflow_dfs:
            for df in outflow_dfs:
                values = df["N"].unique()
                if "n=20" not in values and "n=25" not in values and "n=15" not in values:
                    df["N"] = df["N"].astype("category").cat.reorder_categories(["n=10", "n=5"])
                if "n=20" not in values and "n=25" not in values and "n=15" in values:
                    df["N"] = df["N"].astype("category").cat.reorder_categories(["n=15", "n=10", "n=5"])
                if "n=20" in values and "n=25" not in values:
                    df["N"] = df["N"].astype("category").cat.reorder_categories(["n=20", "n=15", "n=10", "n=5"])
                if "n=20" in values and "n=25" in values:
                    df["N"] = df["N"].astype("category").cat.reorder_categories(["n=25", "n=20", "n=15", "n=10", "n=5"])

            outflow_combined = pd.concat(outflow_dfs)
            get_heatmap(outflow_combined)


def get_heatmap(df):
    phos_obs = 0
    phos_dl = 0

    tss_obs = 0
    tss_dl = 0

    copp_obs = 0
    copp_dl = 0

    if "Phosphorus" in df["Analyte"].values:
        phos_obs = df[df["Analyte"] == "Phosphorus"]["NumObservations"].iloc[0]
        phos_dl = df[df["Analyte"] == "Phosphorus"]["NumObservationsBelowDetect"].iloc[0]

    if "TSS" in df["Analyte"].values:
        tss_obs = df[df["Analyte"] == "TSS"]["NumObservations"].iloc[0]
        tss_dl = df[df["Analyte"] == "TSS"]["NumObservationsBelowDetect"].iloc[0]

    if "Copper" in df["Analyte"].values:
        copp_obs = df[df["Analyte"] == "Copper"]["NumObservations"].iloc[0]
        copp_dl = df[df["Analyte"] == "Copper"]["NumObservationsBelowDetect"].iloc[0]

    df["rpd"] = df["rpd"].astype(np.int)
    df["rpd_w_percent"] = df["rpd"].astype(str) + "%"
    df["Quartile"] = df["Quartile"].astype("category")

    df["Analyte"].replace({
            "Copper" : "Dissolved Copper\nBelow-detection = " + str(copp_dl) + "\nFinal Observations = " + str(copp_obs),
            "Phosphorus" : "Total Phosphorus\nBelow-detection = " + str(phos_dl) + "\nFinal Observations = " + str(phos_obs),
            "TSS" : "Total Suspended Solids\nBelow-detection = " + str(tss_dl) + "\nFinal Observations = " + str(tss_obs)
            }, inplace = True)

    heatmap_colors = ["#1cac78", "#b2ec5d", "#c5e384", "#ffae42", "#ffa343", "#ff7f49", "#ff5349", "#ff2b2b", "#fc2847", "#cb4154"]

    heatmap = (ggplot(df, aes("Quartile", "N", fill = "rpd", label = "rpd_w_percent")) +
                      geom_tile(aes(height = 1, width = 1)) +
                      ggtitle(df["Location"].iloc[0] + "\n" + df["Flow"].iloc[0] + "\n") +
                      scale_fill_gradientn(colors = heatmap_colors, limits = [0, 100], breaks = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]) +
                      geom_text(color = "black", size = 22, family = "serif") +
                      xlab("Quartile") +
                      ylab("Number of Random Events Monitored") +
                      facet_wrap("Analyte") +
                      labs(fill = "Relative Percent Difference (%)") +
                      theme(figure_size = (16, 6), legend_key_size = 35, legend_direction = "horizontal", legend_position = (0.78, 1.21), text = element_text(color = "black", size = 22, family = "serif"), legend_title_align = "center", plot_title = element_text(size = 24)) +
                      theme(panel_border = element_blank(), panel_grid_major = element_blank(), panel_grid_minor = element_blank(), panel_background = element_blank()) +
                      scale_x_discrete(breaks = [0.1, 0.25, 0.5, 0.75, 0.9], expand = (0.1, 0.1)) +
                      scale_y_discrete(breaks = ["n=25", "n=20", "n=15", "n=10", "n=5"], limits = df["N"].unique().tolist()[::-1], drop = False, expand = (0.1, 0.1))
    )

    LocAnalyte = os.path.join(heatmapPlotPath, df["Location"].iloc[0])
    if not os.path.exists(LocAnalyte):
        os.makedirs(LocAnalyte)

    heatmap.save(os.path.join(LocAnalyte, df["Location"].iloc[0] + "_" + df["Flow"].iloc[0] + ".png"))

create_output()
create_heatmap()

shutil.make_archive(os.path.join(os.getcwd(), exportPath, "HeatMapData"), "zip", heatmapDataPath)
shutil.make_archive(os.path.join(os.getcwd(), exportPath, "HeatMapPlots"), "zip", heatmapPlotPath)
shutil.make_archive(os.path.join(os.getcwd(), exportPath, "QuartilepPlots"), "zip", quartilePlotPath)
shutil.make_archive(os.path.join(os.getcwd(), exportPath, "ScatterData"), "zip", scatterDataPath)
