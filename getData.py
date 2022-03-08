##################################################
## Author: Matt McGauley
## Email: mmcgau01@villanova.edu
##################################################

import os
import glob
import csv
import re
import shutil
import numpy as np
import pandas as pd
import rpy2
from scipy.stats import variation
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

pandas2ri.activate()
numpy2ri.activate()
utils = importr("utils")
utils.chooseCRANmirror(ind=1)
utils.install_packages("gamlss")
gamlss = importr("gamlss")
base = importr("base")

rpy2.robjects.r["options"](warn=-1)


dataPath = os.path.join(os.getcwd(), "Data")
if not os.path.exists(dataPath):
    os.makedirs(dataPath)

locationsDataPath = os.path.join(dataPath, "LocationsData")
if not os.path.exists(locationsDataPath):
    os.makedirs(locationsDataPath)

zoneDataPath = os.path.join(dataPath, "ZoneData")
if not os.path.exists(zoneDataPath):
    os.makedirs(zoneDataPath)

outputPath = os.path.join(dataPath, "OutputData")
if not os.path.exists(outputPath):
    os.makedirs(outputPath)

resampledPath = os.path.join(dataPath, "ResampledData")
if not os.path.exists(resampledPath):
    os.makedirs(resampledPath)

metadataPath = os.path.join(dataPath, "MetaData")
if not os.path.exists(metadataPath):
    os.makedirs(metadataPath)

exportPath = os.path.join(os.getcwd(), "Export")
if not os.path.exists(exportPath):
    os.makedirs(exportPath)


def getSkew(output_df):
    skewValue = output_df["EMC"].skew()
    if abs(skewValue) > 1:
        if skewValue < 0:
            return "heavy_left_skew"
        if skewValue > 0:
            return "heavy_right_skew"
    if abs(skewValue) < 1 and abs(skewValue) >= 0.5:
        if skewValue < 0:
            return "moderate_left_skew"
        if skewValue > 0:
            return "moderate_right_skew"
    if abs(skewValue) < 0.5:
        return "no_skew"


def getDistribution(output_df):
    try:
        fit = base.invisible(gamlss.fitDist(
            output_df["EMC"].values, k=2, type="realplus", trace=False, try_gamlss=True))
        fit_type = dict(zip(fit.names, map(list, list(fit))))["family"][0]
    except TypeError:
        return None

    fit = base.invisible(gamlss.fitDist(
        output_df["EMC"].values, k=2, type="realplus", trace=False, try_gamlss=True))
    fit_type = dict(zip(fit.names, map(list, list(fit))))["family"][0]

    return fit_type


def getNumMonths(original):
    calc_average = original.copy()
    calc_average["datesample"] = pd.to_datetime(calc_average["datesample"])
    calc_average["month"] = calc_average["datesample"].dt.month
    calc_average["year"] = calc_average["datesample"].dt.year

    return(int(calc_average.groupby(calc_average["year"])["month"].nunique().mean()))

def getInduction(output_df):
    induction = max(output_df["EMC"])/min(output_df["EMC"])
    return induction

def getMetaData(output_df, original, fileName, dirName, analyte, flow):
    below_detect = original[original["wqqualifier"] == "U"]["value_subhalfdl"] * 2
    mdl = below_detect.median()

    if mdl == np.nan:
        metaData_df = pd.DataFrame({"Number_of_raw_Observations": len(original),
                                    "Number_of_below_detect_Observations": len(original[original["wqqualifier"] == "U"]),
                                    "Number_of_excluded_Observations": len(original[original["initialscreen_flag"] == "No"]),
                                    "Number_of_final_Observations": len(output_df),
                                    "Number_of_Obs_near_dl": np.nan,
                                    "mdl": mdl,
                                    "skewnessCategory": getSkew(output_df),
                                    "skewnessValue": output_df["EMC"].skew(),
                                    "coefficientOfVariation": variation(output_df["EMC"]),
                                    "averageAnnualMonthsSampled": getNumMonths(original),
                                    "inductionRatio": getInduction(output_df),
                                    "distributionCategory": getDistribution(output_df)},
                                   index=[0])
    else:
        metaData_df = pd.DataFrame({"Number_of_raw_Observations": len(original),
                                    "Number_of_below_detect_Observations": len(original[original["wqqualifier"] == "U"]),
                                    "Number_of_excluded_Observations": len(original[original["initialscreen_flag"] == "No"]),
                                    "Number_of_final_Observations": len(output_df),
                                    "Number_of_Obs_near_dl": len(output_df[(output_df["EMC"] >= mdl - (0.5 * mdl)) & (output_df["EMC"] <= mdl + (0.5 * mdl))]),
                                    "mdl": mdl,
                                    "skewnessCategory": getSkew(output_df),
                                    "skewnessValue": output_df["EMC"].skew(),
                                    "coefficientOfVariation": variation(output_df["EMC"]),
                                    "averageAnnualMonthsSampled": getNumMonths(original),
                                    "inductionRatio": getInduction(output_df),
                                    "distributionCategory": getDistribution(output_df)},
                                   index=[0])

    LocAnalyteMetaData = os.path.join(metadataPath, fileName + "_" + analyte)
    if not os.path.exists(LocAnalyteMetaData):
        os.makedirs(LocAnalyteMetaData)

    metaData_df.to_csv(os.path.join(LocAnalyteMetaData, dirName + "_"
                                    + fileName + "_" + flow + "_" + analyte.lower() + ".csv"), index=False)


def resampleData(df, type, dirName, fileName, analyte):
    numObs = len(df)

    if numObs == 15:
        sampleSizes = [5, 10]
    if numObs <= 20 and numObs > 15:
        sampleSizes = [5, 10, 15]
    if numObs <= 25 and numObs > 20:
        sampleSizes = [5, 10, 15, 20]
    if numObs > 25:
        sampleSizes = [5, 10, 15, 20, 25]

    ResampledLocAnalyte = os.path.join(resampledPath, fileName + "_" + analyte)
    if not os.path.exists(ResampledLocAnalyte):
        os.makedirs(ResampledLocAnalyte)

    for sampleSize in sampleSizes:
        resampled = pd.DataFrame()
        for i in range(1, 501):
            sample = df.sample(n=sampleSize)
            sample["Run"] = i
            resampled = resampled.append(sample, ignore_index=True)

        resampled.to_csv(os.path.join(ResampledLocAnalyte, dirName + "_" + fileName + "_"
                                      + type + "_" + analyte.lower() + "_resampled_" + str(sampleSize) + ".csv"), index=False)


def getInflow(df, dirName, fileName, analyte, names):
    original_inflow = df[df["msname"].isin(names)]
    inflow = df[(df["msname"].isin(names)) & (
        df["initialscreen_flag"] != "No")]

    if inflow.empty:
        return
    else:
        output_df = pd.DataFrame()
        output_df["EMC"] = inflow["value_subhalfdl"]
        output_df["Location"] = fileName
        output_df["Analyte"] = analyte
        output_df["Flow"] = "Inflow"

        output_df = output_df.dropna()
        output_df = output_df[output_df["EMC"] > 0]

        if len(output_df) < 15:
            return
        else:
            getMetaData(output_df, original_inflow,
                        fileName, dirName, analyte, "inflow")
            LocAnalyte = os.path.join(outputPath, fileName + "_" + analyte)
            if not os.path.exists(LocAnalyte):
                os.makedirs(LocAnalyte)

            output_df.to_csv(os.path.join(LocAnalyte, dirName + "_"
                                          + fileName + "_inflow_" + analyte.lower() + ".csv"), index=False)
            resampleData(output_df, "inflow", dirName, fileName, analyte)


def getOutflow(df, dirName, fileName, analyte, names):
    original_outflow = df[df["msname"].isin(names)]
    outflow = df[(df["msname"].isin(names)) & (
        df["initialscreen_flag"] != "No")]

    if outflow.empty:
        return
    else:
        output_df = pd.DataFrame()
        output_df["EMC"] = outflow["value_subhalfdl"]
        output_df["Location"] = fileName
        output_df["Analyte"] = analyte
        output_df["Flow"] = "Outflow"

        output_df = output_df.dropna()
        output_df = output_df[output_df["EMC"] > 0]

        if len(output_df) < 15:
            return
        else:
            getMetaData(output_df, original_outflow,
                        fileName, dirName, analyte, "outflow")
            LocAnalyte = os.path.join(outputPath, fileName + "_" + analyte)
            if not os.path.exists(LocAnalyte):
                os.makedirs(LocAnalyte)

            output_df.to_csv(os.path.join(LocAnalyte, dirName + "_"
                                          + fileName + "_outflow_" + analyte.lower() + ".csv"), index=False)
            resampleData(output_df, "outflow", dirName, fileName, analyte)


def createDataframes():
    desiredDirs = ["2", "3", "6", "1", "9", "NSWQ"]

    data = os.path.join(dataPath, "ZoneData")
    subdirs = glob.glob(os.path.join(data, "*"))

    in_names = pd.read_csv(os.path.join(locationsDataPath, "in_msnames.csv"), names=[
                           "InflowName"])["InflowName"].to_list()
    out_names = pd.read_csv(os.path.join(locationsDataPath, "out_msnames.csv"), names=[
                            "OutflowName"])["OutflowName"].to_list()

    for dir in subdirs:
        dirName = dir.split("/")[-1]

        if dirName in desiredDirs:
            files = glob.glob(os.path.join(dir, "*.xlsx"))
            for file in files:
                fileName = os.path.basename(file).split(".")[0]

                copperSheet = None
                tssSheet = None
                phosphorusSheet = None
                keySheet = "key"

                key = pd.read_excel(file, sheet_name=keySheet)

                key["analyte"] = key["analyte"].str.lower()

                copperSheet = np.where(key.loc[key["analyte"] == "copper, dissolved"]["sheet"].empty,
                                       copperSheet, key.loc[key["analyte"] == "copper, dissolved"]["sheet"])
                phosphorusSheet = np.where(key.loc[key["analyte"] == "phosphorus as p, total"]["sheet"].empty,
                                           phosphorusSheet, key.loc[key["analyte"] == "phosphorus as p, total"]["sheet"])
                tssSheet = np.where(key.loc[key["analyte"] == "total suspended solids, ns"]["sheet"].empty,
                                    tssSheet, key.loc[key["analyte"] == "total suspended solids, ns"]["sheet"])

                if phosphorusSheet.size != 0:
                    phosphorus = pd.read_excel(
                        file, sheet_name=str(phosphorusSheet[0]))
                    getInflow(phosphorus, dirName, fileName,
                              "Phosphorus", in_names)
                    getOutflow(phosphorus, dirName, fileName,
                               "Phosphorus", out_names)

                if tssSheet.size != 0:
                    tss = pd.read_excel(file, sheet_name=str(tssSheet[0]))
                    getInflow(tss, dirName, fileName, "TSS", in_names)
                    getOutflow(tss, dirName, fileName, "TSS", out_names)

                if copperSheet.size != 0:
                    copper = pd.read_excel(
                        file, sheet_name=str(copperSheet[0]))
                    getInflow(copper, dirName, fileName, "Copper", in_names)
                    getOutflow(copper, dirName, fileName, "Copper", out_names)

        else:
            continue


createDataframes()

shutil.make_archive(os.path.join(os.getcwd(), exportPath, "LocationsData"), "zip", locationsDataPath)
shutil.make_archive(os.path.join(os.getcwd(), exportPath,
                                 "ZoneData"), "zip", zoneDataPath)
shutil.make_archive(os.path.join(os.getcwd(), exportPath,
                                 "OutputData"), "zip", outputPath)
shutil.make_archive(os.path.join(os.getcwd(), exportPath,
                                 "ResampledData"), "zip", resampledPath)
shutil.make_archive(os.path.join(os.getcwd(), exportPath,
                                 "MetaDataData"), "zip", metadataPath)
