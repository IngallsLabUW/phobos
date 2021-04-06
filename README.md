
<!-- README.md is generated from README.Rmd. Please edit that file -->

# phobos

<!-- badges: start -->
<!-- badges: end -->

phobos is a package comprising of custom functions that accompany the
Ingalls Lab MARS Project.

The MARS project is an easy way to identify and rank unknown mass
features (MFs) in your data. It consists of three major sections:

-   A central database containing compiled unknown mass features with mz
    (mass/charge), rt (retention time, seconds) and ms2 data detected by
    the Ingalls lab
-   *phobos*: A series of processes and functions for annotating,
    ranking, and scoring the unknown mass features
-   A central database containing annotated mass features from previous
    MARS missions

Flexible, updatable, searchable, rankable, exportable.

## Installation

You can install the development version of phobos from
[Github](https://github.com) with:

``` r
devtools::install_github("IngallsLabUW/phobos")
```

Load phobos:

``` r
install.packages("phobos")
```

## Description

A common problem in metabolomics is how to handle unknown but
determinable mass features. The output of non-targeted metabolite runs
produces an abundance of spectra that can be discretely separated and
isolated, but often remain unidentified.

The phobos package, designed as part of the Ingalls MARS Project *TODO
LINK TO FUTURE GITHUB WEBSITE HERE?*, prioritizes simplicity,
efficiency, and confidence in metabolite identification. phobos adheres
to the proposed minimum reporting standards for chemical analysis as
laid out in the Metabolomics Standards Initiative (MSI) [2007
paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3772505/pdf/nihms504189.pdf).
In addition to utilizing the [Ingalls Lab
Standards](https://github.com/IngallsLabUW/Ingalls_Standards) for
foremost confidence annotations, phobos provides functionality for
incorporating spectral comparison to third-party sources such as
\[MassBank of North America
(MoNA)\](<https://mona.fiehnlab.ucdavis.edu/>,
[Metlin](https://metlin.scripps.edu/landing_page.php?pgcontent=mainPage),
*TODO METLIN TBD* and [KEGG](https://www.genome.jp/kegg/) *TODO KEGG
TBD*.

*Maybe can remove this section?* There are currently several available
packages in R that handle metabolite identification, including the
well-established
[xcms](https://bioconductor.org/packages/release/bioc/html/xcms.html),
as well as newly-released packages such as
[RHermes](https://rogerginber.github.io/RHermes/articles/RHermes_UserGuide.html).

## Usage

phobos simplifies the identification process by transforming the MSI
standards
[paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3772505/pdf/nihms504189.pdf)
into steps for annotation. Users begin with experimental (unknown)
values, and systematically compare those values to various dataframes of
theoretical (known) values. Each step corresponds to a custom phobos
function.

    #> ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──
    #> ✓ ggplot2 3.3.3     ✓ purrr   0.3.4
    #> ✓ tibble  3.0.6     ✓ dplyr   1.0.4
    #> ✓ tidyr   1.1.2     ✓ stringr 1.4.0
    #> ✓ readr   1.4.0     ✓ forcats 0.5.1
    #> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    #> x dplyr::filter() masks stats::filter()
    #> x dplyr::lag()    masks stats::lag()

The user begins with a dataframe of experimental data, transformed *by
the user* into a standardized format. It is important to pay attention
to both the column names and the class; the code will not run if they
are not in the correct configuration!

-   “MassFeature”: An identifier for a unique mass feature, character.
-   “mz”: The mz value, numeric.
-   “rt”: The retention time, in seconds, numeric.
-   “column”: Column the mass feature was run on, character.
-   “z”: The ion mode, numeric.
-   “MS2”: MS2 data for those compounds that have it, character.
    -   **Important**: The MS2 data must be in the filtered, scaled and
        concatenated format of “mz, intensity;”, as below. For example
        purposes, the MS2 column has been cut off at two instances of
        “mz, intensity;” but real data should have many more entries.
        For more information on how to get your MS2 data into the
        correct format, *TODO PUT INFORMATION HERE*.

For clarity, phobos also adds an additional identifier column to the
dataframe titled “compound\_experimental”: a simple keys column
consisting of 1-n rows.

``` r
library(phobos)

knitr::kable(head(experimental.values, 5), caption = "Experimental Data")
```

| compound\_experimental | MassFeature     |  mz |  rt | column |   z | MS2                             |
|-----------------------:|:----------------|----:|----:|:-------|----:|:--------------------------------|
|                      1 | I102.0549R12.14 | 102 | 729 | HILIC  |  -1 | 59.01263, 100; 74.99261, 41.1   |
|                      2 | Uracil          | 111 | 240 | HILIC  |  -1 | 110.97466, 100; 110.03485, 39.1 |
|                      3 | I124.0064R10.88 | 124 | 653 | HILIC  |  -1 | 124.01125, 100; 123.90107, 58.1 |
|                      4 | Isethionic Acid | 125 | 432 | HILIC  |  -1 | 124.99041, 100; 94.97973, 27.2  |
|                      5 | I126.904R4.8    | 127 | 288 | HILIC  |  -1 | 126.90401, 100; 85.03849, 1.4   |

Experimental Data

``` r
knitr::kable(theoretical.values[c(1, 9, 18, 25), ], caption = "Theoretical Data")
```

|     | compound   |  mz |   rt | column |   z | MS2                             |
|:----|:-----------|----:|-----:|:-------|----:|:--------------------------------|
| 1   | Methionine | 150 |  478 | HILIC  |   1 | 104.05333, 100; 133.03199, 69.3 |
| 9   | Alanine    |  90 |  678 | HILIC  |   1 | 72.08152, 100; 90.09258, 78.5   |
| 18  | Arginine   | 175 | 1076 | HILIC  |   1 | 175.11914, 100; 116.071, 27.1   |
| 25  | Asparagine | 133 |  712 | HILIC  |   1 | NA                              |

Theoretical Data

|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| \#\#\# Confidence Level 1 Create a dataframe annotated for Confidence Level 1 by comparing experimental and theoretical data, using the AnnotateConfidenceLevel1() function.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `r Confidence.Level.1 <- AnnotateConfidenceLevel1(experimental.values = experimental.values, theoretical.values = theoretical.values, mz.flexibility = 0.02, rt.flexibility = 30)`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| *TODO: fix display of this table*                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| `r knitr::kable(head(Confidence.Level.1, 10))`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| \|MassFeature \| compound\_experimental\|compound\_theoretical \| mz\_experimental\| mz\_theoretical\| rt\_sec\_experimental\| rt\_sec\_theoretical\|column\_experimental \|column\_theoretical \| z\_experimental\| z\_theoretical\|MS2\_experimental \|MS2\_theoretical \| mz\_similarity\_score\| rt\_similarity\_score\| ppm\_mass\_error\| MS2\_cosine\_similarity\| total\_similarity\_score\| confidence\_rank\|confidence\_source \| \|:—————\|———————:\|:——————–\|—————:\|————–:\|——————-:\|——————:\|:——————-\|:——————\|————–:\|————-:\|:——————————-\|:——————————\|——————-:\|——————-:\|————–:\|———————:\|———————-:\|—————:\|:—————–\| \|I102.0549R12.14 \| 1\|NA \| 102\| NA\| 729\| NA\|HILIC \|NA \| -1\| NA\|59.01263, 100; 74.99261, 41.1 \|NA \| NA\| NA\| NA\| NA\| NA\| NA\|NA \| \|Uracil \| 2\|Uracil \| 111\| 111\| 240\| 240\|HILIC \|HILIC \| -1\| -1\|110.97466, 100; 110.03485, 39.1 \|111.01879, 100; 111.04389, 2.4 \| 1\| 1.00\| 0.0\| 0.00\| 67\| 1\|Ingalls\_Standards \| \|Uracil \| 2\|Uracil \| 111\| 111\| 240\| 240\|HILIC \|HILIC \| -1\| -1\|110.97466, 100; 110.03485, 39.1 \|111.01884, 100; 111.05526, 1.9 \| 1\| 1.00\| 0.0\| 0.00\| 67\| 1\|Ingalls\_Standards \| \|Uracil \| 2\|Uracil \| 111\| 111\| 240\| 240\|HILIC \|HILIC \| -1\| -1\|110.97466, 100; 110.03485, 39.1 \|111.01889, 100; 111.04259, 2.1 \| 1\| 1.00\| 0.0\| 0.00\| 67\| 1\|Ingalls\_Standards \| \|Uracil \| 2\|Uracil \| 111\| 111\| 240\| 240\|HILIC \|HILIC \| -1\| -1\|110.97466, 100; 110.03485, 39.1 \|111.01886, 100; 111.04398, 2 \| 1\| 1.00\| 0.0\| 0.00\| 67\| 1\|Ingalls\_Standards \| \|Uracil \| 2\|Uracil \| 111\| 111\| 240\| 240\|HILIC \|HILIC \| -1\| -1\|110.97466, 100; 110.03485, 39.1 \|111.01885, 100; 111.04394, 1.8 \| 1\| 1.00\| 0.0\| 0.00\| 67\| 1\|Ingalls\_Standards \| \|Uracil \| 2\|Uracil \| 111\| 111\| 240\| 240\|HILIC \|HILIC \| -1\| -1\|110.97466, 100; 110.03485, 39.1 \|111.01881, 100; 111.04394, 1.9 \| 1\| 1.00\| 0.0\| 0.00\| 67\| 1\|Ingalls\_Standards \| \|Uracil \| 2\|Uracil \| 111\| 111\| 240\| 240\|HILIC \|HILIC \| -1\| -1\|110.97466, 100; 110.03485, 39.1 \|68.01299, 100; 52.17011, 58.5 \| 1\| 1.00\| 0.0\| 0.00\| 67\| 1\|Ingalls\_Standards \| \|Uracil \| 2\|Uracil \| 111\| 111\| 240\| 240\|HILIC \|HILIC \| -1\| -1\|110.97466, 100; 110.03485, 39.1 \|111.01888, 100; 67.02911, 88.4 \| 1\| 1.00\| 0.0\| 0.00\| 67\| 1\|Ingalls\_Standards \| \|I124.0064R10.88 \| 3\|Taurine \| 124\| 124\| 653\| 668\|HILIC \|HILIC \| -1\| -1\|124.01125, 100; 123.90107, 58.1 \|124.00623, 100; 79.95613, 76.5 \| 1\| 0.88\| 3.6\| 0.75\| 88\| 1\|Ingalls\_Standards \| |

### Confidence Level 2

Compounds that receive a confidence level of 2 are considered putatively
annotated compounds: there are no available chemical reference
standards, but they have good MS1 and MS2 spectral similarity with
public/commercial spectral libraries. External libraries used for
confidence level 2 are MoNA and Metlin.

The comparison MoNA spectra can be downloaded directly from the [public
website](https://mona.fiehnlab.ucdavis.edu/downloads). Unlike the
Ingalls Standards sheets, the downloaded MoNA sheets require
transformation from JSON to csv, and are large enough to need external
storage. The code to transform the JSONs to csvs is included in the
phobos package. The pre-downloaded and cleaned MoNA csvs live on the
shared Ingalls Google Drive [MARS
project](https://drive.google.com/drive/folders/1lzIsDJJ7EyDpJrCTLdspHc0KS6zxUIDQ)
folder.

``` r
Confidence.Level.2 <- AnnotateMoNAConfidenceLevel2(Confidence.Level.1 = Confidence.Level.1, mz.flexibility = 0.02, rt.flexibility = 30)
#> [1] "Making potential candidates"
#> Joining, by = c("compound_experimental", "z_experimental", "mz_experimental", "MS2_experimental")
```

|                                                                                                                                                                                                                                 |
|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| \#\#\# Confidence Level 3                                                                                                                                                                                                       |
| Compounds that receive a confidence level of 3 are also considered putatively annotated with fewer restrictions. The MARS project and the phobos package define a confidence level 3 as an MS1 match from MoNA, Metlin or KEGG. |

### Confidence Level 4

Confidence level 4 is everything else. It should be filtered for
possible contaminants and made sure that the potential mass feature is
quality.
