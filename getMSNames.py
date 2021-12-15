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

exportPath = os.path.join(os.getcwd(), "Export")
if not os.path.exists(exportPath):
    os.makedirs(exportPath)

dataPath = os.path.join(os.getcwd(), "Data")
if not os.path.exists(dataPath):
    os.makedirs(dataPath)

locationsDataPath = os.path.join(dataPath, "LocationsData")
if not os.path.exists(locationsDataPath):
    os.makedirs(locationsDataPath)

zoneDataPath = os.path.join(dataPath, "ZoneData")
if not os.path.exists(zoneDataPath):
    os.makedirs(zoneDataPath)

def get_msnames():
    desiredDirs = [2, 3, 6, 1, 9]

    subdirs = glob.glob(os.path.join(zoneDataPath, "*"))

    in_msnames = []
    out_msnames = []

    for dir in subdirs:
        dirName = dir[-1]

        if int(dirName) in desiredDirs:
            files = glob.glob(os.path.join(dir, "*.xlsx"))
            for file in files:
                fileName = os.path.basename(file).split(".")[0]

                copperSheet = None
                tssSheet = None
                phosphorusSheet = None
                keySheet = "key"

                key = pd.read_excel(file, sheet_name = keySheet)

                key["analyte"] = key["analyte"].str.lower()

                copperSheet = np.where(key.loc[key["analyte"] == "copper, dissolved"]["sheet"].empty, copperSheet, key.loc[key["analyte"] == "copper, dissolved"]["sheet"])
                phosphorusSheet = np.where(key.loc[key["analyte"] == "phosphorus as p, total"]["sheet"].empty, phosphorusSheet, key.loc[key["analyte"] == "phosphorus as p, total"]["sheet"])
                tssSheet = np.where(key.loc[key["analyte"] == "total suspended solids, ns"]["sheet"].empty, tssSheet, key.loc[key["analyte"] == "total suspended solids, ns"]["sheet"])


                if phosphorusSheet.size != 0:
                    phosphorus = pd.read_excel(file, sheet_name = str(phosphorusSheet[0]))
                    in_msnames.extend(phosphorus[phosphorus["msname"].str.contains("in", flags = re.IGNORECASE)]["msname"].unique())
                    out_msnames.extend(phosphorus[phosphorus["msname"].str.contains("out", flags = re.IGNORECASE)]["msname"].unique())

                if tssSheet.size != 0:
                    tss = pd.read_excel(file, sheet_name = str(tssSheet[0]))
                    in_msnames.extend(tss[tss["msname"].str.contains("in", flags = re.IGNORECASE)]["msname"].unique())
                    out_msnames.extend(tss[tss["msname"].str.contains("out", flags = re.IGNORECASE)]["msname"].unique())

                if copperSheet.size != 0:
                    copper = pd.read_excel(file, sheet_name = str(copperSheet[0]))
                    in_msnames.extend(copper[copper["msname"].str.contains("in", flags = re.IGNORECASE)]["msname"].unique())
                    out_msnames.extend(copper[copper["msname"].str.contains("out", flags = re.IGNORECASE)]["msname"].unique())

        else:
            continue

    in_df = pd.DataFrame(list(set(in_msnames)), columns = ["InflowName"])
    in_df.to_csv(os.path.join(locationsDataPath, "in_msnames.csv"), index = False)

    out_df = pd.DataFrame(list(set(out_msnames)), columns = ["OutflowName"])
    out_df.to_csv(os.path.join(locationsDataPath, "out_msnames.csv"), index = False)

get_msnames()

shutil.make_archive(os.path.join(os.getcwd(), exportPath, "LocationsData"), "zip", locationsDataPath)
shutil.make_archive(os.path.join(os.getcwd(), exportPath, "ZoneData"), "zip", zoneDataPath)
