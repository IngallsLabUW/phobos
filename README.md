
<!-- README.md is generated from README.Rmd. Please edit that file -->

# phobos

<!-- badges: start -->
<!-- badges: end -->

phobos is an R package of custom functions that accompany the Ingalls
Lab MARS Project.

The MARS project is an easy way to identify and rank unknown mass
features (MFs) in your data. It consists of three major sections:

-   A central database containing compiled mass features detected by the
    Ingalls lab with mz (mass/charge), rt (retention time, seconds),
    column, charge and ms2 data
-   *phobos*: A package of functions for annotating, ranking, and
    scoring the unknown mass features, as well as other common data
    transformations.
-   A central database containing annotated mass features from previous
    MARS missions. **TODO: This may end up being combined with the first
    database as one complete “homebase”.**

**Flexible, updatable, searchable, rankable, exportable.**

## Installation

You can install the development version of phobos from
[Github](https://github.com) with:

``` r
devtools::install_github("IngallsLabUW/phobos")
```

Once installed, phobos can be loaded like any other package:

``` r
library("phobos")
```

## Description

A common problem in metabolomics is how to handle unknown but
determinable mass features. The output of non-targeted metabolite runs
produces an abundance of spectra that can be discretely separated and
isolated, but often remain unidentified, or identified with non-standard
confidence.

The phobos package, designed as part of the Ingalls MARS (Metabolite
Annotation, Rank and Sort) Project **TODO: LINK TO FUTURE GITHUB WEBSITE
HERE**, prioritizes simplicity, efficiency, and confidence in metabolite
identification. To encourage confidence consistency, phobos adheres to
the proposed minimum reporting standards for chemical analysis as laid
out in the Metabolomics Standards Initiative (MSI) [2007
paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3772505/pdf/nihms504189.pdf).
In addition to utilizing the [Ingalls Lab
Standards](https://github.com/IngallsLabUW/Ingalls_Standards) for
foremost confidence annotations (known as “Confidence Level 1”), phobos
provides functionality for incorporating spectral comparison to
third-party sources such as [MassBank of North America
(MoNA)](https://mona.fiehnlab.ucdavis.edu/),
[Metlin](https://metlin.scripps.edu/landing_page.php?pgcontent=mainPage),
**TODO: METLIN TBD** and [KEGG](https://www.genome.jp/kegg/) **TODO:
KEGG TBD**.

For a complete description of the MARS annotation process, as well as
all references for the equations used in the phobos package, please
refer to the **TODO link website or github page here**.

## Usage

**TODO: Put visual description here or somewhere else in the
usage/description? An overview of some kind needs to be here**

phobos simplifies the identification process by transforming each
individual MSI rank into a step for annotation. Each step corresponds to
a custom phobos function. Users begin with experimental (unknown)
values, and systematically compare those to theoretical (known) values.

The user begins with a dataframe of experimental data, transformed *by
the user* into a standardized format. It is important to pay attention
to both the column names and the class; the code will not run if they
are not in the correct configuration!

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
        entries. For more information on how to get your MS2 data into
        the correct format, **TODO PUT INFORMATION HERE**.

For clarity, phobos also adds an additional primary column to the
dataframe titled “primary\_key”: a simple key column consisting of 1-n
rows. Throughout the project, the search parameter column names in the
dataframes (mz, rt, z, column, MS2) will be modified according to their
source. For example, all “theoretical” mz values, aka values that are
derived from a standards sheet, are identified as “mz\_theoretical”,
where as experimental values are identified as “mz\_experimental”, MoNA
values are identified as “mz\_massbank”, etc. This also applies to
confidence level rankings and similarity scores: A column titled
mz\_similarity\_score1 refers to the similarity score calculated during
the Confidence Level 1 annotation.

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
read.csv(“example\_data/Theoretical\_Data.csv”) command. The complete
Ingalls Standards csv is located on the [Ingalls Standards
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
the first step, either NA or Ingalls\_Standards). **TODO: add more
specifics here about interpreting the output**

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

|     | MassFeature         | primary\_key | compound\_theoretical | massbank\_match                                                                                                                                                                                                                                                    | massbank\_ppm | massbank\_cosine\_similarity | confidence\_rank | confidence\_source |
|:----|:--------------------|-------------:|:----------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------:|-----------------------------:|:-----------------|:-------------------|
| 1   | I102.0549R12.14     |            1 | NA                    | NA                                                                                                                                                                                                                                                                 |            NA |                           NA | NA               | NA                 |
| 2   | Uracil              |            2 | Uracil                | NA                                                                                                                                                                                                                                                                 |            NA |                           NA | 1                | Ingalls\_Standards |
| 45  | L-Pyroglutamic acid |            6 | Taurine, D4           | L-Pyroglutamic acid; pGlu; pyroGlu; Pyroglutamate; Pidolic acid; L-Glutimic acid; L-5-Oxoproline; (S)-(?)-2-Pyrrolidone-5-carboxylic acid; (S)-5-Oxo-2-pyrrolidinecarboxylic acid; L-a-Aminoglutaric Acid Lactam; L-5-Oxo-2-pyrrolidinecarboxylic acid ID:PR100586 |           4.3 |                         0.78 | 2                | MassBank           |
| 88  | N-acetyltaurine     |           12 | N-acetyltaurine       | NA                                                                                                                                                                                                                                                                 |            NA |                           NA | 1                | Ingalls\_Standards |
| 124 | I192.9131R3.72      |           18 | NA                    | NA                                                                                                                                                                                                                                                                 |            NA |                           NA | NA               | NA                 |
| 125 | I192.9808R2.93      |           19 | NA                    | NA                                                                                                                                                                                                                                                                 |            NA |                           NA | NA               | NA                 |
| 201 | I277.9097R1.79      |           44 | NA                    | NA                                                                                                                                                                                                                                                                 |            NA |                           NA | NA               | NA                 |
| 209 | I325.0188R11.15     |           52 | NA                    | NA                                                                                                                                                                                                                                                                 |            NA |                           NA | NA               | NA                 |

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
defines a confidence level 3 as a good MS1 match from MoNA, Metlin or
KEGG within the column and z specifications.

Confidence level 3 uses the same external data as confidence level 2, so
no additional downloading is required.

<!-- # ```{r Confidence Level 3} -->
<!-- # Confidence.Level.3 <- read.csv("example_data/Example_ConfidenceLevel3.csv") -->
<!-- # ``` -->
<!-- #  -->
<!-- # ```{r Display CL3, echo=FALSE} -->
<!-- # CL3 <- Confidence.Level.3 %>% -->
<!-- #   arrange(primary_key) %>% -->
<!-- #   select(MassFeature, massbank_match3, mz_experimental, MH_mass_MoNA, confidence_rank, confidence_source) %>% -->
<!-- #   unique() -->
<!-- # knitr::kable(CL3[c(1, 26, 66, 314), ], caption = "Confidence Level 3") -->
<!-- # ``` -->

------------------------------------------------------------------------

### Confidence Level 4

Confidence level 4 is everything else! **TODO: It should be filtered for
possible contaminants and made sure that the potential mass feature is
quality.**

\*\*TODO: This is commented out while the join explosions are dealt
with, to see if all the data can be kept in the package/github. Will be
brought back in later.

<!-- # ```{r load CL3 data, echo=FALSE} -->
<!-- # Confidence.Level.3 <- read.csv("example_data/Example_ConfidenceLevel3.csv") -->
<!-- # ``` -->
<!-- #  -->
<!-- # ```{r Confidence Level 4} -->
<!-- # Confidence.Level.4 <- AnnotateConfidenceLevel4(Confidence.Level.3 = Confidence.Level.3) -->
<!-- # ``` -->
<!-- #  -->
<!-- # **TODO: Display Confidence Level 4 here, need to deal with some join explosions first. -->
<!-- # ```{r Display CL4, echo=FALSE} -->
<!-- # CL4 <- Confidence.Level.4 %>% -->
<!-- #   arrange(primary_key) %>% -->
<!-- #   select(MassFeature, compound_experimental:MoNA_Names, MH_mass_experimental, MH_mass_MoNA, mz_similarity_score3:confidence_source) %>% -->
<!-- #   unique() -->
<!-- # knitr::kable(CL4[c(1, 26, 73), ], caption = "Confidence Level 4") -->
<!-- # ``` -->
