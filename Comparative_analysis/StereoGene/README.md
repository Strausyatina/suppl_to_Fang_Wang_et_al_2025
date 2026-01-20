Stereogene[V2.50] for calculating the positional correlation between and among the chromatin marks and RNA expression.

We calculated the correlation between and among the chromatin marks and rna as below cases:
- Dmod_v_Dmod: Changing chromatin modification between KO and WT of Mark A vs Mark B.
- Dmod_v_Drna: Changing chromatin modification between KO and WT of Mark A vs Changing RNA expression between KO and WT.
- mod_v_mod: Chromatin modification mark A from specific condition to the mark B from the same condition.
- mod_v_rna: Chromatin modification mark A from specific condition to the RNA expression` from the same condition.

In the respective directories, runstereogene.sh is the script processing the input, running the stereogene and delivering the output.

Input description:
  - As the input, stereogene expect the pairs of the bedGRaph or model files. (example bedgraph files are located in each case comparison)
  - model files: It is a text file describing the linear combination of the input tracks, example model files are located in each case comparison)
  - Update the path to the stereogene pacakge.
  - genome.chrom.sizes: contains the list of all chromosomes in a genome along with their lenght. The chrome.size file used during this analysis is provied in the stereoGene direcotry.

Parameter description:
  - wSize: stands for window size. It is one of the key params impacting resolution of correlation analysis, performance, and noise handling hence controls how the stereogene process and compares genomic signals.
  - bins: It is closely related to wSize. Bins are the sigment of the genome that stereogene used for summarizing data
  - kernels: It deteremines how much weight to give to positions around a central point when comparing two genomic tracks.Its a smoothing function.
  - trackPath: It is the location of all the tracks and/or model files
  - profPath: It contains all the profiles genereted from stereogene
  - chrom: requires the path to genome.chrom.sizes
  - report: location of the stereogene report
  - plotType: specify the format of the plot (ex. pdf)
  - writeDistr: specify if need the raw output of the distribution.

Output description:
 - It will output the directory for each comparison with corresponding comparison name.
 - Resulting directory contains the following file formats:
   -  .bkg: contains the correlation values for the background
   -  .fg: contaoins the correlation values for the foreground.
   -  .dist: contains the distance, its foreground, and background correlation value.
   -  .r: contain the R script to generate the plot.
   -  report: This directory contains the output of the .r script.

<ins>Customized Rscripts</ins>:<br>
Directory **Script** contains the DistancePlot.R, DensityPlot.R
- DistancePlot.R:
    - Input: It takes the case directory (ex. Dmod_v_Dmod) as input, please change this path in the script, also change the modality type of your interest. Then script looks for .dist file from all the comparisons from this case.
    - Output: distance plots as a pdf file.
    - Usage: Rscript DistancePlot.R
- DensityPlot.R
    - Input: It takes the case directory (ex. Dmod_v_Dmod) as input, please change this path in the script, also change the modality type of your interest. Then script looks for .fg and .bkg file from all the comparisons from this case.
    - Output: density plots as a pdf file.
    - Usages: Rscript DensityPlot.R



