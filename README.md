
<!-- README.md is generated from README.Rmd. Please edit that file -->

# phobos

<!-- badges: start -->
<!-- badges: end -->

phobos is an R package of custom functions that accompany the Ingalls
Lab MARS Project.

The MARS project helps identify and rank unknown mass features (MFs) in
untargeted metabolomics data. It consists of two major sections:

-   A central “homebase” database containing compiled mass features
    detected by the Ingalls lab with *m/z* (mass/charge), retention time
    (in seconds), chromatographic column, ion charge, and fragmentation
    data, as well as all previous MARS missions(?) and metadata (What
    metadata might be included?).
-   This package, *phobos*: A package of functions for annotating,
    ranking, and scoring the unknown mass features, as well as other
    common data transformations associated with the MARS project.

**phobos aims to be: Flexible, updatable, searchable, rankable,
exportable.**

![MARS](example_data/MARS_Workflow.png)

## Installation

You can install the development version of phobos from
[Github](https://github.com/IngallsLabUW/phobos) with:

``` r
devtools::install_github("IngallsLabUW/phobos")
#> Downloading GitHub repo IngallsLabUW/phobos@HEAD
#> 
#>      checking for file ‘/private/var/folders/7j/gxsnrqzs7qg_g_7hp143b8xr0000gn/T/RtmpvR4rXp/remotes5bee166a79f3/IngallsLabUW-phobos-0eb15cc/DESCRIPTION’ ...  ✓  checking for file ‘/private/var/folders/7j/gxsnrqzs7qg_g_7hp143b8xr0000gn/T/RtmpvR4rXp/remotes5bee166a79f3/IngallsLabUW-phobos-0eb15cc/DESCRIPTION’
#>   ─  preparing ‘phobos’:
#>      checking DESCRIPTION meta-information ...  ✓  checking DESCRIPTION meta-information
#>   ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#>   ─  building ‘phobos_0.1.0.tar.gz’
#>      
#> 
```

Once installed, phobos can be loaded like any other package:

``` r
library(phobos)
```

## Description

A common problem in metabolomics is how to handle unknown but
determinable mass features. The output of non-targeted metabolite runs
produces an abundance of spectra that can be discretely separated and
isolated, but often remain unidentified, or identified with non-standard
confidence.

The phobos package, designed as part of the Ingalls MARS (Metabolite
Annotation, Rank and Sort) Project, prioritizes simplicity, efficiency,
and confidence in metabolite identification. To encourage confidence
consistency, phobos adheres to the proposed minimum reporting standards
for chemical analysis as laid out in the Metabolomics Standards
Initiative (MSI) [2007
paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3772505/pdf/nihms504189.pdf).
In addition to utilizing the [Ingalls Lab
Standards](https://github.com/IngallsLabUW/Ingalls_Standards) for
foremost confidence annotations (known as “Confidence Level 1”), phobos
provides functionality for incorporating spectral comparison to
third-party sources such as [MassBank of North America
(MoNA)](https://mona.fiehnlab.ucdavis.edu/),,
[KEGG](https://www.genome.jp/kegg/), and
[Metlin](https://metlin.scripps.edu/landing_page.php?pgcontent=mainPage),
**TODO: METLIN TBD**.

For a complete description of the MARS annotation process, as well as
all references for the equations used in the phobos package, please
refer to the MARS Design Doc.

## Usage

phobos simplifies the identification process by transforming each
individual MSI rank into a step for annotation. Users begin with
experimental values, and systematically compare those to theoretical
values. With each successive function, a new dataframe is created with
an additional level of annotation.

Before running the annotation functions, The experimental data is
transformed *by the user* into a standardized format. It is important to
pay attention to both the column names and the class; the code will not
run if they are not in the correct configuration!

Below is a list of the columns required to use MARS. The dataframe must
contain only these columns, in this order.

-   “MassFeature”: Any identifier for a unique mass feature, character.
-   “mz”: The m/z value, numeric.
-   “rt”: The retention time, in seconds, numeric.
-   “column”: The column the mass feature was run on, character. There
    are two options for this variable: “HILIC” or “RP” (short for
    Reverse Phase). In the Ingalls lab, cyano compounds are only run in
    positive mode, therefore the combination of an “RP” column
    observation and a “1” polarity observation will result in the
    correct analysis.
-   “z”: The polarity, numeric.
-   “MS2”: MS2 data for those compounds that have it, character.
    -   **Important**: The MS2 data must be in the filtered, scaled and
        concatenated format of “mz, intensity;”, as below. For example
        purposes, the MS2 column below has been cut off at two instances
        of “mz, intensity;” but real data will likely have many more
        entries. There are several functions included in this package to
        help modify the data to various forms: ConcatToScan() for
        transforming concatenated data to scantable format,
        ScanToConcat() for doing the opposite, and ScaleMS2() for
        scaling the MS2 data to a maximum intensity of 100 and minimum
        intensity of 0.5

For clarity, phobos also adds an additional primary column to the
dataframe titled “primary\_key”: a simple key column consisting of 1-n
rows. (What’s 1-n here?) Throughout the project, the search parameter
column names in the dataframes (mz, rt, z, column, MS2) will be modified
according to their source. For example, all “theoretical” mz values, aka
values that are derived from a standards sheet, are identified as
“mz\_theoretical”, where as experimental values are identified as
“mz\_experimental”, MoNA values are identified as “mz\_massbank”, etc.
This also applies to confidence level rankings and similarity scores: A
column titled mz\_similarity\_score1 refers to the similarity score
calculated during the Confidence Level 1 annotation.

Below is an example of an experimental dataframe in the correct format.

| MassFeature     |  mz |  rt | column |   z | MS2                             |
|:----------------|----:|----:|:-------|----:|:--------------------------------|
| I102.0549R12.14 | 102 | 729 | HILIC  |  -1 | 59.01263, 100; 74.99261, 41.1   |
| Uracil          | 111 | 240 | HILIC  |  -1 | 110.97466, 100; 110.03485, 39.1 |
| I124.0064R10.88 | 124 | 653 | HILIC  |  -1 | 124.01125, 100; 123.90107, 58.1 |
| Isethionic Acid | 125 | 432 | HILIC  |  -1 | 124.99041, 100; 94.97973, 27.2  |
| I126.904R4.8    | 127 | 288 | HILIC  |  -1 | 126.90401, 100; 85.03849, 1.4   |

Experimental Data

------------------------------------------------------------------------

### Confidence Level 1

In order to receive a confidence level 1 annotation, a mass feature must
be compared to an authentic chemical standard using two independent and
orthogonal data. Within the MARS project, a mass feature must match to a
chemical standard by column, charge, mz, retention time, and fall below
a user-defined ppm.

The experimental data is compared to the theoretical data (the Ingalls
Standards csv) for confidence level 1. A version of the theoretical
data, modified for the phobos package, is included in the phobos
example\_data/ subdirectory and can be imported using the
read.csv(“example\_data/Example\_Theoretical\_Data.csv”) command. The
complete Ingalls Standards csv is located on the [Ingalls Standards
Github](https://github.com/IngallsLabUW/Ingalls_Standards/blob/master/Ingalls_Lab_Standards_NEW.csv)
and can be downloaded for a more updated version. Like the experimental
data, the theoretical data must be in a specific format before
comparison. **TODO: Right now the github standards doesn’t include the
MS2 data, it will once consensus is reached!**

|     | compound   |  mz |   rt | column |   z | MS2                             |
|:----|:-----------|----:|-----:|:-------|----:|:--------------------------------|
| 1   | Methionine | 150 |  478 | HILIC  |   1 | 104.05333, 100; 133.03199, 69.3 |
| 9   | Alanine    |  90 |  678 | HILIC  |   1 | 72.08152, 100; 90.09258, 78.5   |
| 18  | Arginine   | 175 | 1076 | HILIC  |   1 | 175.11914, 100; 116.071, 27.1   |
| 25  | Asparagine | 133 |  712 | HILIC  |   1 | NA                              |

Theoretical Data

To annotate the experimental data for confidence level 1, use the
AnnotateConfidenceLevel1() function. The values for mz.flexibility and
rt.flexibility are in Daltons and seconds, respectively.

``` r
Confidence.Level.1 <- AnnotateConfidenceLevel1(experimental.values = experimental.values, theoretical.values = theoretical.values, mz.flexibility = 0.02, rt.flexibility = 30)
#> [1] "Columns are correctly named and ordered."
#> [1] "Below are your column classes. Please check them to make sure they are correct!"
#> $MassFeature
#> [1] "character"
#> 
#> $mz
#> [1] "numeric"
#> 
#> $rt
#> [1] "numeric"
#> 
#> $column
#> [1] "character"
#> 
#> $z
#> [1] "integer"
#> 
#> $MS2
#> [1] "character"
#> Joining, by = c("MassFeature", "primary_key")
#> Joining, by = c("MassFeature", "primary_key")
```

| primary\_key | compound\_theoretical | mz\_similarity\_score1 | rt\_similarity\_score1 | ppm\_mass\_error1 | MS2\_cosine\_similarity1 | total\_similarity\_score1 | confidence\_rank | confidence\_source |
|-------------:|:----------------------|-----------------------:|-----------------------:|------------------:|-------------------------:|--------------------------:|-----------------:|:-------------------|
|            1 | NA                    |                     NA |                     NA |                NA |                       NA |                        NA |               NA | NA                 |
|            2 | Uracil                |                      1 |                   1.00 |               0.0 |                     0.00 |                        67 |                1 | Ingalls\_Standards |
|            3 | Taurine               |                      1 |                   0.88 |               3.6 |                     0.75 |                        88 |                1 | Ingalls\_Standards |
|            4 | Isethionic Acid       |                      1 |                   1.00 |               0.0 |                     1.00 |                       100 |                1 | Ingalls\_Standards |
|            5 | NA                    |                     NA |                     NA |                NA |                       NA |                        NA |               NA | NA                 |

Confidence Level 1

In the above table, the results have been given an m/z, rt, MS2, and ppm
scores, which are used to create a total similarity score. Depending on
the results, a mass feature is given a confidence rank and source (for
the first step, either NA or Ingalls\_Standards). For more information
about interpreting the MARS output, please jump to the bottom of this
README to the “Interpreting Output” section.

------------------------------------------------------------------------

### Confidence Level 2

Compounds that receive a confidence level of 2 are considered putatively
annotated compounds; there are no available chemical reference
standards, but they have good MS1 and MS2 spectral similarity with
public/commercial spectral libraries. External libraries used for
confidence level 2 are MoNA and Metlin **TODO: incude Metlin**. These
libraries are updated every quarter **TODO: include details on
updates**.

In order to compare the experimental data to MoNA, it **must** have
already been annotated for Confidence Level 1 using the
AnnotateConfidenceLevel1() function. The output of that function is then
run through the AnnotateMoNAConfidenceLevel2() function.

The ready-to-compare MoNA spectra can be found on the shared Ingalls Lab
Google Drive under Collaborative\_Projects –&gt; MARS\_Project –&gt;
MoNA\_RelationalSpreadsheets, but the JSON data can also be downloaded
directly from the [public
website](https://mona.fiehnlab.ucdavis.edu/downloads) for the most
recent version. They are downloaded in JSON form and so require
transformation to csv for input to phobos. The code to transform the
JSONs to csv is included in the phobos package, in the **TODO MAKE THIS
FUNCTION**.

``` r
Confidence.Level.2 <- read.csv("example_data/Example_ConfidenceLevel2.csv")
```

|     | MassFeature      | primary\_key | compound\_theoretical | massbank\_match2 | massbank\_ppm | MS2\_cosine\_similarity2 | confidence\_rank | confidence\_source |
|:----|:-----------------|-------------:|:----------------------|:-----------------|--------------:|:-------------------------|:-----------------|:-------------------|
| 1   | I102.0549R12.14  |            1 | NA                    | NA               |            NA | NA                       | NA               | NA                 |
| 2   | Uracil           |            2 | Uracil                | NA               |            NA | NA                       | 1                | Ingalls\_Standards |
| 95  | Sulfolactic Acid |           15 | Sulfolactic Acid      | NA               |            NA | NA                       | 1                | Ingalls\_Standards |
| 119 | I217.93R10.87    |           24 | NA                    | NA               |            NA | NA                       | NA               | NA                 |

Confidence Level 2

The output of AnnotateMoNAConfidenceLevel2() produces a dataframe that
includes both Confidence Level 1 & 2 annotations. Often, a mass feature
will be annotated on both confidence levels, which results in a
concatenated column separated by semicolon (e.g, “1; 2” means a mass
feature has been matched on both confidence levels). This format is
echoed in the “confidence\_source” column. For example, if a mass
feature has a rank of “1; 2”, then its accompanying confidence\_source
would be “Ingalls\_Standards; MassBank”.

------------------------------------------------------------------------

### Confidence Level 3

Compounds that receive a confidence level of 3 are also considered
putatively annotated based on fewer restrictions. The MARS project
defines a confidence level 3 as a good MS1 match from MoNA, Metlin, or
KEGG within the column and z specifications.

For each database used, there is a separate, ordered set of functions.
For a full explanation of the order, please refer to the image at the
top of this README. Because databases are used, the user needs to ensure
that they have all of the required information accessible. The data from
MoNA is the same as in the previous step. The data from KEGG is located
in the same place as the MoNA: the shared Ingalls Lab Google Drive under
Collaborative\_Projects –&gt; MARS\_Project –&gt; KEGG\_Spreadsheets.
**TODO: When ready, add more info here about updating/creating the KEGG
data.**

MoNA:

``` r
Confidence.Level.3_MoNA <- read.csv("example_data/Example_ConfidenceLevel3_MoNA.csv")
```

|     | MassFeature         | primary\_key | compound\_theoretical | massbank\_match2                                                                                                                                                                                                                                                                     | massbank\_match3                                                                                                                                                                                                                                                     | MH\_mass\_experimental | MH\_mass\_MoNA | mz\_similarity\_score3 | confidence\_rank | confidence\_source |
|:----|:--------------------|-------------:|:----------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------:|---------------:|-----------------------:|:-----------------|:-------------------|
| 1   | I102.0549R12.14     |            1 | NA                    | NA                                                                                                                                                                                                                                                                                   | NA                                                                                                                                                                                                                                                                   |                     NA |             NA |                     NA | NA               | NA                 |
| 2   | Uracil              |            2 | Uracil                | NA                                                                                                                                                                                                                                                                                   | NA                                                                                                                                                                                                                                                                   |                     NA |             NA |                     NA | 1                | Ingalls\_Standards |
| 26  | I126.904R4.8        |            5 | NA                    | NA                                                                                                                                                                                                                                                                                   | NA                                                                                                                                                                                                                                                                   |                     NA |             NA |                     NA | NA               | NA                 |
| 27  | L-Pyroglutamic acid |            6 | Taurine, D4           | ID: PR100586; L-Pyroglutamic acid; pGlu; pyroGlu; Pyroglutamate; Pidolic acid; L-Glutimic acid; L-5-Oxoproline; (S)-(?)-2-Pyrrolidone-5-carboxylic acid; (S)-5-Oxo-2-pyrrolidinecarboxylic acid; L-a-Aminoglutaric Acid Lactam; L-5-Oxo-2-pyrrolidinecarboxylic acid ID:ID: PR100586 | ID: PR100586; L-Pyroglutamic acid; pGlu; pyroGlu; Pyroglutamate; Pidolic acid; L-Glutimic acid; L-5-Oxoproline; (S)-(?)-2-Pyrrolidone-5-carboxylic acid; (S)-5-Oxo-2-pyrrolidinecarboxylic acid; L-a-Aminoglutaric Acid Lactam; L-5-Oxo-2-pyrrolidinecarboxylic acid |                    128 |            128 |                      1 | 2; 3             | MassBank; MassBank |

Confidence Level 3: MoNA

KEGG:

``` r
Confidence.Level.3_KEGG <- read.csv("example_data/Example_ConfidenceLevel3_KEGG.csv")
```

|      | MassFeature     | primary\_key | compound\_theoretical | massbank\_match2 | massbank\_match3 | Compound\_KEGG                                                                                                                                                                                 | mz\_similarity\_scoreKEGG | confidence\_rank | confidence\_source |
|:-----|:----------------|-------------:|:----------------------|:-----------------|:-----------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------:|:-----------------|:-------------------|
| 1    | I102.0549R12.14 |            1 | NA                    | NA               | NA               | cpd:C00155; cpd:C00334; cpd:C01026; cpd:C01205; cpd:C02261; cpd:C02356; cpd:C02721; cpd:C03284; cpd:C03665; cpd:C05145; cpd:C05330; cpd:C05689; cpd:C05698; cpd:C11109; cpd:C11735; cpd:C19514 |                         1 | 3                | NA; KEGG           |
| NA   | NA              |           NA | NA                    | NA               | NA               | NA                                                                                                                                                                                             |                        NA | NA               | NA                 |
| NA.1 | NA              |           NA | NA                    | NA               | NA               | NA                                                                                                                                                                                             |                        NA | NA               | NA                 |
| NA.2 | NA              |           NA | NA                    | NA               | NA               | NA                                                                                                                                                                                             |                        NA | NA               | NA                 |

Confidence Level 3: KEGG

------------------------------------------------------------------------

### Confidence Level 4

Confidence level 4 is everything else!

Any mass feature that has not yet been assigned a confidence level in a
previous step will receive a confidence level of 4.

``` r
# Beautiful! Love to see the function used here, very clear
Confidence.Level.4 <- AnnotateConfidenceLevel4(Confidence.Level.3_KEGG)
```

|      | MassFeature     | primary\_key | compound\_theoretical | massbank\_match2 | massbank\_match3 | Compound\_KEGG                                                                                                                                                                                 | mz\_similarity\_scoreKEGG | confidence\_rank | confidence\_source       |
|:-----|:----------------|-------------:|:----------------------|:-----------------|:-----------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------:|:-----------------|:-------------------------|
| 1    | I102.0549R12.14 |            1 | NA                    | NA               | NA               | cpd:C00155; cpd:C00334; cpd:C01026; cpd:C01205; cpd:C02261; cpd:C02356; cpd:C02721; cpd:C03284; cpd:C03665; cpd:C05145; cpd:C05330; cpd:C05689; cpd:C05698; cpd:C11109; cpd:C11735; cpd:C19514 |                         1 | 3                | NA; KEGG                 |
| 117  | Gluconic Acid   |           21 | Gluconic Acid         | NA               | NA               | cpd:C16356; cpd:C16360                                                                                                                                                                         |                         1 | 1; 3             | Ingalls\_Standards; KEGG |
| NA   | NA              |           NA | NA                    | NA               | NA               | NA                                                                                                                                                                                             |                        NA | NA               | NA                       |
| NA.1 | NA              |           NA | NA                    | NA               | NA               | NA                                                                                                                                                                                             |                        NA | NA               | NA                       |
| NA.2 | NA              |           NA | NA                    | NA               | NA               | NA                                                                                                                                                                                             |                        NA | NA               | NA                       |
| NA.3 | NA              |           NA | NA                    | NA               | NA               | NA                                                                                                                                                                                             |                        NA | NA               | NA                       |

Confidence Level 4

------------------------------------------------------------------------

### Common Contaminants

In order to check against common contaminants, download the
“Common\_Contaminants.csv” data from the shared Ingalls Google Drive in
the MARS\_Project folder.

Run the checkContaminants() function using the final annotated data and
the downloaded Common Contaminants csv.

------------------------------------------------------------------------

### Interpreting Output

Once the experimental data has been completely through the MARS project,
the final output is a set of data frames with additional columns
containing various matching parameters: matches, similarity scores,
ranks, sources, etc.

In general, if a column as a suffix of “1”, e.g. mz\_similarity\_score1,
that column is associated with Confidence Level 1. The same logic
applies to columns with suffixes of 2 or 3.

When deciding between potential matches, start by prioritizing
confidence level, and then sort within level by similarity scores. If a
massfeature has received matches on confidence levels 1, 2 and 3, assume
that the confidence level 1 output is the best one to choose, and then
start comparing similarity scores within that confidence level to narrow
down matches to the best option.

When looking at similarity scores, there are several to choose from
depending on the confidence level. mz\_similarity\_score,
rt\_similarity\_score, ppm\_mass\_error, MS2\_similarity\_score, and
total\_similarity\_score are all potential filtering options for the
user to consider.

“Good” similarity scores are the following: - mz\_similarity\_score1
&gt; 0.9 - rt\_similarity\_score1 &gt; 0.75 - ppm\_mass\_error1 &lt; 7

MS2 cosine similarity does not act as a filtering score in the code due
to collision energies and NAs, and are instead included a column for
user review. The same is true for the total similarity score, which
includes all available similarity scores (see Design Doc for more
information).
