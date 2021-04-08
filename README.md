
<!-- README.md is generated from README.Rmd. Please edit that file -->

# phobos

<!-- badges: start -->
<!-- badges: end -->

phobos is a package of custom functions that accompany the Ingalls Lab
MARS Project.

The MARS project is an easy way to identify and rank unknown mass
features (MFs) in your data. It consists of three major sections:

-   A central database containing compiled unknown mass features with mz
    (mass/charge), rt (retention time, seconds), column, charge and ms2
    data detected by the Ingalls lab
-   *phobos*: A package of functions for annotating, ranking, and
    scoring the unknown mass features
-   A central database containing annotated mass features from previous
    MARS missions

**Flexible, updatable, searchable, rankable, exportable.**

## Installation

You can install the development version of phobos from
[Github](https://github.com) with:

``` r
devtools::install_github("IngallsLabUW/phobos")
```

Once installed, phobos should load like any other package:

``` r
install.packages("phobos")
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
The MARS project In addition to utilizing the [Ingalls Lab
Standards](https://github.com/IngallsLabUW/Ingalls_Standards) for
foremost confidence annotations, phobos provides functionality for
incorporating spectral comparison to third-party sources such as
[MassBank of North America (MoNA)](https://mona.fiehnlab.ucdavis.edu/),
[Metlin](https://metlin.scripps.edu/landing_page.php?pgcontent=mainPage),
**TODO: METLIN TBD** and [KEGG](https://www.genome.jp/kegg/) **TODO:
KEGG TBD**.

**TODO: better to have the full description all in the README or link to
a github website?**

## Usage

**TODO: Put visual description here or somewhere else in the
usage/description? An overview of some kind needs to be here**

phobos simplifies the identification process by transforming the each
individual MSI standard into a step for annotation. Each step
corresponds to a custom phobos function. Users begin with experimental
(unknown) values, and systematically compare those to theoretical
(known) values.

The user begins with a dataframe of experimental data, transformed *by
the user* into a standardized format. It is important to pay attention
to both the column names and the class; the code will not run if they
are not in the correct configuration!

Below is a list of the columns required to use MARS. The dataframe must
contain only these columns, in this order.

-   “MassFeature”: Any identifier for a unique mass feature, character.
-   “mz”: The mz value, numeric.
-   “rt”: The retention time, in seconds, numeric.
-   “column”: Column the mass feature was run on, character.
-   “z”: The ion mode, numeric.
-   “MS2”: MS2 data for those compounds that have it, character.
    -   **Important**: The MS2 data must be in the filtered, scaled and
        concatenated format of “mz, intensity;”, as below. For example
        purposes, the MS2 column below has been cut off at two instances
        of “mz, intensity;” but real data should have many more entries.
        For more information on how to get your MS2 data into the
        correct format, **TODO PUT INFORMATION HERE**.

For clarity, phobos also adds an additional primary column to the
dataframe titled “compound\_experimental”: a simple keys column
consisting of 1-n rows.

Below is an example of an experimental dataframe in the correct format.

| compound\_experimental | MassFeature     |  mz |  rt | column |   z | MS2                             |
|-----------------------:|:----------------|----:|----:|:-------|----:|:--------------------------------|
|                      1 | I102.0549R12.14 | 102 | 729 | HILIC  |  -1 | 59.01263, 100; 74.99261, 41.1   |
|                      2 | Uracil          | 111 | 240 | HILIC  |  -1 | 110.97466, 100; 110.03485, 39.1 |
|                      3 | I124.0064R10.88 | 124 | 653 | HILIC  |  -1 | 124.01125, 100; 123.90107, 58.1 |
|                      4 | Isethionic Acid | 125 | 432 | HILIC  |  -1 | 124.99041, 100; 94.97973, 27.2  |
|                      5 | I126.904R4.8    | 127 | 288 | HILIC  |  -1 | 126.90401, 100; 85.03849, 1.4   |

Experimental Data

|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| \#\#\# Confidence Level 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| In order to receive a confidence level 1 annotation, a mass feature must be compared to an authentic chemical standard using two independent and orthogonal data. Within the MARS project, a mass feature must match to a chemical standard by column, charge, mz, retention time, and ppm.                                                                                                                                                                                                                                                                                                                                             |
| The experimental data is compared to the theoretical data (the Ingalls Standards csv) for confidence level 1. This csv is located on the [Ingalls Standards Github](https://github.com/IngallsLabUW/Ingalls_Standards/blob/master/Ingalls_Lab_Standards_NEW.csv) and can be downloaded. Like the experimental data, the theoretical data must be in a specific format before comparison. **TODO: The modified csv should probably just live on the MARS github, or wherever the external data like the MoNA sheets live. Users won’t need to edit it every time.**                                                                      |
| Table: Theoretical Data                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| \| \|compound \| mz\| rt\|column \| z\|MS2 \| \|:–\|:———-\|—:\|—-:\|:——\|–:\|:——————————-\| \|1 \|Methionine \| 150\| 478\|HILIC \| 1\|104.05333, 100; 133.03199, 69.3 \| \|9 \|Alanine \| 90\| 678\|HILIC \| 1\|72.08152, 100; 90.09258, 78.5 \| \|18 \|Arginine \| 175\| 1076\|HILIC \| 1\|175.11914, 100; 116.071, 27.1 \| \|25 \|Asparagine \| 133\| 712\|HILIC \| 1\|NA \|                                                                                                                                                                                                                                                         |
| To annotate the experimental data for confidence level 1, use the AnnotateConfidenceLevel1() function.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| `r Confidence.Level.1 <- AnnotateConfidenceLevel1(experimental.values = experimental.values, theoretical.values = theoretical.values, mz.flexibility = 0.02, rt.flexibility = 30)`                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| Table: Confidence Level 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| \| compound\_experimental\|compound\_theoretical \| mz\_similarity\_score\| rt\_similarity\_score\| ppm\_mass\_error\| MS2\_cosine\_similarity\| total\_similarity\_score\| confidence\_rank\|confidence\_source \| \|———————:\|:——————–\|——————-:\|——————-:\|————–:\|———————:\|———————-:\|—————:\|:—————–\| \| 1\|NA \| NA\| NA\| NA\| NA\| NA\| NA\|NA \| \| 2\|Uracil \| 1\| 1.00\| 0.0\| 0.00\| 67\| 1\|Ingalls\_Standards \| \| 3\|Taurine \| 1\| 0.88\| 3.6\| 0.75\| 88\| 1\|Ingalls\_Standards \| \| 4\|Isethionic Acid \| 1\| 1.00\| 0.0\| 1.00\| 100\| 1\|Ingalls\_Standards \| \| 5\|NA \| NA\| NA\| NA\| NA\| NA\| NA\|NA \| |
| In the above table, the results have been given an m/z, rt, MS2, and ppm scores, which are used to create a total similarity score. Depending on the results, a mass feature is given a confidence rank and source (for the first step, either NA or Ingalls\_Standards).                                                                                                                                                                                                                                                                                                                                                               |

### Confidence Level 2

Compounds that receive a confidence level of 2 are considered putatively
annotated compounds: there are no available chemical reference
standards, but they have good MS1 and MS2 spectral similarity with
public/commercial spectral libraries. External libraries used for
confidence level 2 are MoNA and Metlin **TODO: incude Metlin**. These
libraries are updated every quarter **TODO: include details on
updates**.

In order to compare the experimental data to MoNA, it **must** have
first been annotated for Confidence Level 1 using the
AnnotateConfidenceLevel1() function. The output of that function can
then be run through the AnnotateMoNAConfidenceLevel2() function.

**TODO: Download and define the MoNA spectra. For now they are included
in the data but won’t be in the final package**

The ready-to-compare MoNA spectra can be found on the *TODO: CONFIRM
Ingalls Google Drive/GitHub website*, but they can also be downloaded
directly from the [public
website](https://mona.fiehnlab.ucdavis.edu/downloads) for the most
recent version. They are downloaded in JSON form and so require
transformation to csv. The code to transform the JSONs to csv is
included in the phobos package, in the **TODO MAKE THIS FUNCTION**.

``` r
Confidence.Level.2 <- AnnotateMoNAConfidenceLevel2(Confidence.Level.1 = Confidence.Level.1, mz.flexibility = 0.02, rt.flexibility = 30)
#> [1] "Making potential candidates"
#> Joining, by = c("compound_experimental", "z_experimental", "mz_experimental", "MS2_experimental")
```

|     | compound\_experimental | compound\_theoretical | massbank\_match                                                                                                                                                                                                                                                    | massbank\_ppm | massbank\_cosine\_similarity | confidence\_rank | confidence\_source |
|:----|-----------------------:|:----------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------:|-----------------------------:|:-----------------|:-------------------|
| 1   |                      1 | NA                    | NA                                                                                                                                                                                                                                                                 |            NA |                           NA | NA               | NA                 |
| 2   |                      2 | Uracil                | NA                                                                                                                                                                                                                                                                 |            NA |                           NA | 1                | Ingalls\_Standards |
| 45  |                      6 | Taurine, D4           | L-Pyroglutamic acid; pGlu; pyroGlu; Pyroglutamate; Pidolic acid; L-Glutimic acid; L-5-Oxoproline; (S)-(?)-2-Pyrrolidone-5-carboxylic acid; (S)-5-Oxo-2-pyrrolidinecarboxylic acid; L-a-Aminoglutaric Acid Lactam; L-5-Oxo-2-pyrrolidinecarboxylic acid ID:PR100586 |           4.3 |                         0.91 | 2                | MoNA               |
| 88  |                     12 | N-acetyltaurine       | NA                                                                                                                                                                                                                                                                 |            NA |                           NA | 1                | Ingalls\_Standards |
| 124 |                     18 | NA                    | NA                                                                                                                                                                                                                                                                 |            NA |                           NA | NA               | NA                 |
| 125 |                     19 | NA                    | NA                                                                                                                                                                                                                                                                 |            NA |                           NA | NA               | NA                 |
| 187 |                     46 | NA                    | NA                                                                                                                                                                                                                                                                 |            NA |                           NA | NA               | NA                 |
| 195 |                     54 | NA                    | NA                                                                                                                                                                                                                                                                 |            NA |                           NA | NA               | NA                 |

Confidence Level 2

The output of AnnotateMoNAConfidenceLevel2() produces a dataframe that
includes both Confidence Level 1 & 2 annotations. Often, a compound will
be annotated on both confidence levels, which results in a concatenated
column separated by semicolon (e.g, “1; 2” means a compound has been
matched on both confidence levels).

|                                                                                                                                                                                                                                                                        |
|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| \#\#\# Confidence Level 3                                                                                                                                                                                                                                              |
| Compounds that receive a confidence level of 3 are also considered putatively annotated with fewer restrictions. The MARS project and the phobos package define a confidence level 3 as an MS1 match from MoNA, Metlin or KEGG within the column and z specifications. |
| Confidence level 3 uses the same external data as confidence level 2, so no additional downloading is required.                                                                                                                                                        |

### Confidence Level 4

Confidence level 4 is everything else. It should be filtered for
possible contaminants and made sure that the potential mass feature is
quality.
