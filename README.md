<h1>Instructions for Running the automated BMP Pre-Processing</h1>
<h2>Matt McGauley<sup>1</sup>, Emily Darin<sup>2</sup>, Bridget Wadzuk, Ph.D.<sup>1</sup>, Elizabeth Fassman-Beck, Ph.D.<sup>2</sup></h2>
<sup>1</sup>Villanova Center for Resilient Water Systems
<sup>2</sup>Southern California Coastal Water Research Project



<h2>Getting Started</h2>
1.	Install python if you don’t have it installed already

  1. Go to: https://www.python.org/downloads/release/python-397/
  2. Install the optional features for pip and Python launcher if using Windows

2. Make sure the following packages are installed with pip
  1.	os, glob, csv, re, shutil, warnings (these are part of the standard python library, so they should already be installed)
  2.	numpy
  3.	pandas
  4.	plotnine
  5.	matplotlib
  6.	openpyxl
  7.	xlrd
  8.	ggplot
  9.	rpy2

Install packages on the Command Prompt (Windows) or Terminal (MacOS) like so:
(Windows)   py -m pip install <<name of package>>
(MacOS)      python -m pip install <<name of package>>

See here for more information about installing pip here if you are having issues, it should come pre-installed when you download python. See here about a common issue with Windows where pip is not recognized as a command on the Command Prompt. This shouldn’t be an issue if you install pip when you installed Python.

3.	Create a folder where you’re going to run the scripts for the analysis (name does not matter)
  1.  Put the following in this folder:
    * Cloned github repo
 
You only need to run the following steps once and they must be run in the following order to get everything we need:

Step 1: getMSNames.py
1.	Open Command Prompt (on Windows) or Terminal (MacOS)
2.	Navigate to the path of the folder you created in Step 3 of Getting Started
  1.	Windows: cd /D <<path to your folder>>  (You can copy the path to your folder from the top bar of the file navigator, then paste it where <<path to your folder>> is in this command)
  2.	MacOS: cd <<path to your folder>>  (You can copy the path to your folder by right clicking on the folder in Finder while holding the option key, then clicking “Copy <<folder>> as Pathname”, then paste it where <<path to your folder>> is in this command)
3.	Run getMSNames.py
Now that you’ve navigated to the folder where everything is, type the following:
(Windows)   py getMSNames.py
(MacOS)      python getMSNames.py
**Don’t close the Command Prompt/Terminal Window or you’ll have to renavigate to the folder in each step!**
4.	This will give you two CSV files in the “LocationsData” folder within the “Data” folder in the folder you created in Step 3 of Getting Started, which we need to make a few small edits to.
  1.	Open “in_msnames.csv”
  2.	Delete rows with any names that are not actually inflow locations, getting rid of blank rows when you delete.
  3.	Save when you’re done.
  4.	Open “out_msnames.csv”
  5.	Delete rows with any names that are not actually outflow locations, getting rid of blank rows when you delete.
  6.	Save when you’re done.

This script went through every zones 1, 2, 3, 6, and 9 in ZoneData, found all possible Inflow names (msnames with “in” in them, regardless of case) and all possible Outflow names (msnames with “out” in them, regardless of case), then put them each in separate CSV files. You then went through and kept on the viable ones.
 
Step 2: getData.py
1.	From the same Command Prompt (Windows) or Terminal (MacOS) navigated to the folder where the data and scripts are (See Step 1.1 and 1.2 if you need to do this again), run getData.py
(Windows)   py getData.py
(MacOS)      python getData.py
**Don’t close the Command Prompt/Terminal Window or you’ll have to renavigate to the folder in each step!**

This script took the “in_msnames” and “out_msnames”, went through each of the files in ZoneData in zones 1, 2, 3, 6, and 9, extracted Inflow and Outflow data using the msnames and checked if there was 15 samples available, got rid of NaN observations and flagged observations from the BMP database (blank values, -99999, or “No” in the “initialscreen_flag” column of the input Zone Data files) then created a csv files with the columns EMC (substituting for the one-half detection limit where necessary), Location, Analyte, Flow and put them in subfolders of “OutputData”, organized by location and analyte. Then, another function resampled the data 500 times for sample sizes of 5, 10, and 15 samples and 20 and 25 samples, if there were enough observations. These files, with columns EMC, Location, Analyte, Flow, Run were put in subfolders of “ResampledData”, organized by location and analyte. At the same time, metadata to be displayed on the quartile plots was captured and stored in files in the “MetaData” folder in the “Data” folder. Archives (zip files) of the ZoneData, Output, and Resampled folders were then copied to a folder called “Export”.

 
Step 3: getData.py
1.	From the same Command Prompt (Windows) or Terminal (MacOS) navigated to the folder where the data and scripts are (See Step 1.1 and 1.2 if you need to do this again), run getData.py
(Windows)   py getPlots.py
(MacOS)      python getPlots.py

This script first compiled a list of all the names of the “Output” files and all of the names of the “Resampled” files. Then for each of the “Output” files, using the file name of the Output, the “Resampled” files matching that “Output” were identified. A dataframe was created for the “Output” file. For each of the “Resampled” files, a dataframe was created, with the sample size was added in a new column called “N” in the format “n=(sample size)”. They were then combined into a single “Resampled” dataframe, and N was assigned as factors ordered by sample size (like in R). Then quartiles (0.1, 0.25, 0.5, 0.75, 0.9) were assigned to the data in the “Output” dataframe. The same was done for the “Resampled” dataframe for each combination of (Run, N). Then, for each quartile, a plot was assembled where the “Output” quartile’s EMC was shown as a vertical line, and the distribution (density plot) of EMC’s at that quartile were plotted, labeled by sample size. This is the exact same ggplot code Emily Darin wrote, just reformatted to accept the data with Python syntax rather than R. A plot of each quartile was compiled into a single pdf and put in subfolders of “QuartilePlots” in the “Plots” folder, organized by location and analyte. The process is repeated for each of the “Output” files, matching “Resampled” files and each quartile. Using the metadata, heatmap data was constructed by including the relative percent difference between the original and resampled each quartile for each sample size, which is stored in the “HeatMapData” folder under “Data”. The heatmap data was used to create plots of the relative percent difference between the original and resampled data. They are in the “HeatMapPlots” folder under “Plots” organized by location and flow type. Scatter data is captured for the option to create scatter plots of the relationship between the relative percent difference and the percentage that the resampled data is of the original dataset. This data is in “ScatterData” under “Data”. An Archive (zip files) of the Plots folder was then copied to the folder called “Export” from the previous step.


What you now have:
1.	The raw data (ZoneData)
2.	The data we want extracted from the raw data (OutputData)
3.	The data we want, resampled (ResampledData)
4.  Data for the creation of heatmaps and quartile plots (HeatMapData, ResampledData, and MetaData)
5.	The plots (Plots) — (HeatMapPlots and QuartilePlots)
6.	Archives of the four items above (Export)
7.	Two CSV files with the MSNames for the inflow and outflow locations for the data we extracted (in_msnames.csv and out_msnames.csv in Data/LocationsData)
8.	The scripts (getMSNames.py, getData.py, getPlots.py)
