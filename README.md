
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

------------------------------------------------------------------------

Create a dataframe annotated for Confidence Level 1 by comparing
exerimental and theoretical data, using the AnnotateConfidenceLevel1()
function.

``` r
ConfidenceLevel1 <- AnnotateConfidenceLevel1(experimental.values = experimental.values, theoretical.values = theoretical.values, 
                                             mz.flexibility = 0.02, rt.flexibility = 30)
```

``` r
knitr::kable(head(ConfidenceLevel1, 10))
```

| compound\_experimental | mz\_similarity\_score | rt\_similarity\_score | ppm\_mass\_error | MS2\_cosine\_similarity | total\_similarity\_score | MassFeature     | confidence\_rank | confidence\_source |
|-----------------------:|----------------------:|----------------------:|-----------------:|------------------------:|-------------------------:|:----------------|-----------------:|:-------------------|
|                      1 |                    NA |                    NA |               NA |                      NA |                       NA | I102.0549R12.14 |               NA | NA                 |
|                      2 |                     1 |                  1.00 |              0.0 |                    0.00 |                       67 | NA              |                1 | Ingalls\_Standards |
|                      2 |                     1 |                  1.00 |              0.0 |                    0.00 |                       67 | NA              |                1 | Ingalls\_Standards |
|                      2 |                     1 |                  1.00 |              0.0 |                    0.00 |                       67 | NA              |                1 | Ingalls\_Standards |
|                      2 |                     1 |                  1.00 |              0.0 |                    0.00 |                       67 | NA              |                1 | Ingalls\_Standards |
|                      2 |                     1 |                  1.00 |              0.0 |                    0.00 |                       67 | NA              |                1 | Ingalls\_Standards |
|                      2 |                     1 |                  1.00 |              0.0 |                    0.00 |                       67 | NA              |                1 | Ingalls\_Standards |
|                      2 |                     1 |                  1.00 |              0.0 |                    0.00 |                       67 | NA              |                1 | Ingalls\_Standards |
|                      2 |                     1 |                  1.00 |              0.0 |                    0.00 |                       67 | NA              |                1 | Ingalls\_Standards |
|                      3 |                     1 |                  0.88 |              3.6 |                    0.75 |                       88 | NA              |               NA | NA                 |

2.  Run through Annotate\_Confidence\_Level1.R script. All required
    extra data (Ingalls standards and MS2) are included in this
    repository. When you have produced your confidence level 1
    datasheet, write the csv and save it for the next step.

3.  Move to Annotate\_Confidence\_Level2\_MoNA.R. As in the previous
    step, all the required extra data (scraped MoNA sheets) are
    available in the data\_extra/MoNA\_RelationalSpreadsheets directory.

-   Be aware that some of the column selections/renaming are to avoid
    explosion joins. If you want to see more data included in your
    joins, adjust the selections as you need. This will cause your
    datasheet to balloon!

Over time, phobos increases in efficiency as the user applies the
package to more and more datasets.

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/master/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />
