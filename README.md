# README: M6 Enrichment

# General description

The m6 gene in _Drosophila_ is annotated to code for 6 different isoforms. 
Since the antibody staining and endogenously tagged M6 show a different 
distribution of the protein between tricellular vertices and the bicellular 
membrane, I wanted to assess whether the different isoforms locate 
differently in the cells. For this, I used UAS constructs with either the 
whole genomic region ("full-length", "pan-isoform") or the cDNA for one 
isoform. In addition, I wanted to test the effect of cysteine mutations 
which might prevent a palmitoylation of the protein.


# Project Organization

The project folder is divided into several folders. Those will be explained in
the following.

## data

The `./data` folder contains all the raw data. The folder contains one 
sub-folder for each analyzed construct, e.g. `./data/delta-palm` for the 
endogenous homozygous M6[DPalm] allele. Currently, there are the following 
sub-folders.

* `./data/trap_wt-homo` The endogenous wildtypic M6-GFP construct (homozygous).
* `./data/trap_wt-hetero` The endogenous wildtypic M6-GFP construct 
  (heterozygous).
* `./data/trap_delta-palm-homo` The endogenous Delta-Palm M6-GFP construct 
  (homozygous) generated using CRISPR.
* `./data/trap_delta-palm-hetero` The endogenous Delta-Palm M6-GFP construct 
  (heterozygous) generated using CRISPR.
* `./data/wt` The wildtype pan-isoform M6-GFP transgene under control of the 
  69B driver.
* `./data/delta-palm` The Delta-palm pan-isoform M6-GFP transgene under 
  control of the 69B driver.
* `./data/isoform-b` Isoform B cDNA transgene under control of the 69B driver.
* `./data/isoform-c` Isoform C cDNA transgene under control of the 69B driver.
* `./data/isoform-d` Isoform D cDNA transgene under control of the 69B driver.
* `./data/isoform-e` Isoform E cDNA transgene under control of the 69B driver.
* `./data/isoform-f` Isoform F cDNA transgene under control of the 69B driver.
* `./data/nrg` 
  
Each of these folders contains one folder for each analyzed embryo. The name 
of the folder contains the meta-information of that embryo. The naming 
convention is the following:

`date-<Date>_prot-<Protein>_driver-<Driver>_cons-<Construct>_exp-<Expression>_emb-<Embryo>`

Each embryo folder contains the following:

* An overview image of the embryo with the same name as the folder with a 
  leading underscore.

* The image stack for analysis with the same name as the folder.

* A CSV-file for each analyzed cell. The file contains the coordinates for 
  all vertices of a single cell in the image stack. The X and Y coordinates 
  are given in µm. The files have the following naming convention<br>
  `Cell-<Cell>.csv`
    
* An image file for each analyzed cell. The naming convention is
  `Cell-<Cell>.tif`
  

## output

The `./output` folder contains all intermediate and final output of the 
project. This includes

* `./outout/_enrichment.csv` The raw mean intensities for each analyzed cell 
  along with additional data. This is the output from the Python script.
* `./output/enrichment.xlsx` An Excel file containing all the processed data 
  for each cell and embryo as well as some summary statistics. The same 
  information is also available as separate CSV files.
* Plots in PDF format.


## src

The `./src` folder contains R and Python source code that will be called 
from the actual analysis scripts in the root folder. These files are 
documented separately.
  

# Data Acquisition

Embryos were collected overnight on 25°C. After dechorionation embryos were 
aligned under the fluorescence bino, selected by gut morphology (roughly 
stage 15), and aligned laterally. Images were acquired at the new Leica SP8
(2) with the following settings for the overview:

* *Objective:* HC PL APO CS2 40x/1.30 OIL
* *Resolution:* 2048x1024px
* *Zoom:* 0.75x
* *Scan speed*: 100 lines/sec
* *Laser:* 488 nm with roughly 0.6%&ndash;0.8% laser power.

The image stack for analysis was acquired with the following settings:

* *Objective:* HC PL APO CS2 40x/1.30 OIL
* *Resolution:* 768x768px
* *Zoom:* 5.00x
* *Scan speed*: 100 lines/sec
* *Laser:* 488 nm with roughly 0.3%&ndash;0.6% laser power.
* *Line accumulation:* 2x
* *Z-stack:* 25 slices with a spacing of 0.35 µm (8.31 µm in total)

> Isoform E had a pretty weak signal. The laser power was increased to 
> roughly 1.2% and line accumulation was increased to 3 for the image stack 
> and 2 for the overview image.

This results in a resolution of $13.1915 px/µm$.


# Data Analysis


## Image Analysis

For each embryo roughly 5&ndash;10 cells were analyzed. For each cell, the 
vertex coordinates were selected in ImageJ and exported as a CSV file. The 
actual analysis was done in Python with a self-written script. From the 
vertex coordinates masks for the vertices, the remaining membrane and the 
background (the cytosol of the cell) were determined. The enrichment for 
each cell was then calculated as

$$ E_{cell} = \frac{\bar{I}_{TCJ} - \bar{I}_{BG}}{\bar{I}_{BCJ} - \bar{I}_
{BG}} $$

The enrichment of an embryo was calculated as the mean of the analyzed cells:

$$ E_{embryo} = \frac{1}{n_i} \cdot \sum_{j=1}^{n_i} E_j $$

where $E_j$ is the $j$-th cell of the $i$-embryo.


## Statistical Analysis

The data for the cysteine mutants was first tested for normality using the 
Shapiro-Wilk test. A normal distribution was only rejected for the 
transgenes ($\alpha = 0.05$). For the heterozygous and homozygous trap embryos 
equality of variances was tested using an F-test. Equality was not rejected 
for both ($\alpha = 0.05$). Accordingly, the transgenes were tested using 
Wilcoxon test and the other two compairsons were conducted using two-sample 
t-test with equal variances.