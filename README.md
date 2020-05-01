# Spawning Potential Ratio for black rockfish near Kodiak

General structure is:

scott reproduces the excel spreadsheet - basically a familiarization and training ground

R has 3 files
 - helper
 - vonb
 - brf_spr
 
Run the vonb.R file first - it generates von Bertalanffy parameters for male and female fish and stores the output in the output folder
The brf_spr.R file will estimate the current selectivity, M, F, and SPR for a pooled dataset (e.g., all years are bunched together, weighted by catch by area and year). 
It will also estimate fishing mortality at target SPR level.
 
The TMB folder is all C++ code - play in there at your own peril.

# Reproducibility

Please contact carrie.worton@alaska.gov to obtain the necessary data files stored on the Kodiak wiki:
data/dockside_AWSLM.csv

**Session info:**

```
> devtools::session_info()
- Session info --------------------------------------------------------------
 setting  value                       
 version  R version 3.6.3 (2020-02-29)
 os       Windows >= 8 x64            
 system   x86_64, mingw32             
 ui       RStudio                     
 language (EN)                        
 collate  English_United States.1252  
 ctype    English_United States.1252  
 tz       America/Anchorage           
 date     2020-05-01                  

- Packages ------------------------------------------------------------------
 ! package     * version date       lib source                             
   assertthat    0.2.1   2019-03-21 [1] CRAN (R 3.6.3)                     
   backports     1.1.6   2020-04-05 [1] CRAN (R 3.6.3)                     
   broom         0.5.6   2020-04-20 [1] CRAN (R 3.6.3)                     
   callr         3.4.3   2020-03-28 [1] CRAN (R 3.6.3)                     
   cellranger    1.1.0   2016-07-27 [1] CRAN (R 3.6.3)                     
   cli           2.0.2   2020-02-28 [1] CRAN (R 3.6.3)                     
   colorspace    1.4-1   2019-03-18 [1] CRAN (R 3.6.1)                     
   crayon        1.3.4   2017-09-16 [1] CRAN (R 3.6.3)                     
   DBI           1.1.0   2019-12-15 [1] CRAN (R 3.6.3)                     
   dbplyr        1.4.3   2020-04-19 [1] CRAN (R 3.6.3)                     
   desc          1.2.0   2018-05-01 [1] CRAN (R 3.6.3)                     
   devtools    * 2.3.0   2020-04-10 [1] CRAN (R 3.6.3)                     
   digest        0.6.25  2020-02-23 [1] CRAN (R 3.6.3)                     
   dplyr       * 0.8.5   2020-03-07 [1] CRAN (R 3.6.3)                     
   ellipsis      0.3.0   2019-09-20 [1] CRAN (R 3.6.3)                     
   fansi         0.4.1   2020-01-08 [1] CRAN (R 3.6.3)                     
   farver        2.0.3   2020-01-16 [1] CRAN (R 3.6.3)                     
   forcats     * 0.5.0   2020-03-01 [1] CRAN (R 3.6.3)                     
   fs            1.4.1   2020-04-04 [1] CRAN (R 3.6.3)                     
   funcr       * 0.1.0   2020-04-30 [1] Github (ben-williams/funcr@e0668cd)
   generics      0.0.2   2018-11-29 [1] CRAN (R 3.6.3)                     
   ggplot2     * 3.3.0   2020-03-05 [1] CRAN (R 3.6.3)                     
   glue          1.4.0   2020-04-03 [1] CRAN (R 3.6.3)                     
   gtable        0.3.0   2019-03-25 [1] CRAN (R 3.6.3)                     
   haven         2.2.0   2019-11-08 [1] CRAN (R 3.6.3)                     
   here        * 0.1     2017-05-28 [1] CRAN (R 3.6.3)                     
   hms           0.5.3   2020-01-08 [1] CRAN (R 3.6.3)                     
   httr          1.4.1   2019-08-05 [1] CRAN (R 3.6.3)                     
   jsonlite      1.6.1   2020-02-02 [1] CRAN (R 3.6.3)                     
   labeling      0.3     2014-08-23 [1] CRAN (R 3.6.0)                     
   lattice       0.20-38 2018-11-04 [2] CRAN (R 3.6.3)                     
   lifecycle     0.2.0   2020-03-06 [1] CRAN (R 3.6.3)                     
   lubridate     1.7.8   2020-04-06 [1] CRAN (R 3.6.3)                     
   magrittr      1.5     2014-11-22 [1] CRAN (R 3.6.3)                     
   Matrix        1.2-18  2019-11-27 [2] CRAN (R 3.6.3)                     
   memoise       1.1.0   2017-04-21 [1] CRAN (R 3.6.3)                     
   modelr        0.1.7   2020-04-30 [1] CRAN (R 3.6.3)                     
   munsell       0.5.0   2018-06-12 [1] CRAN (R 3.6.3)                     
   nlme          3.1-144 2020-02-06 [2] CRAN (R 3.6.3)                     
   pillar        1.4.3   2019-12-20 [1] CRAN (R 3.6.3)                     
   pkgbuild      1.0.7   2020-04-25 [1] CRAN (R 3.6.3)                     
   pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 3.6.3)                     
   pkgload       1.0.2   2018-10-29 [1] CRAN (R 3.6.3)                     
   prettyunits   1.1.1   2020-01-24 [1] CRAN (R 3.6.3)                     
   processx      3.4.2   2020-02-09 [1] CRAN (R 3.6.3)                     
   ps            1.3.2   2020-02-13 [1] CRAN (R 3.6.3)                     
   purrr       * 0.3.4   2020-04-17 [1] CRAN (R 3.6.3)                     
   R6            2.4.1   2019-11-12 [1] CRAN (R 3.6.3)                     
   Rcpp          1.0.4.6 2020-04-09 [1] CRAN (R 3.6.3)                     
   readr       * 1.3.1   2018-12-21 [1] CRAN (R 3.6.3)                     
   readxl        1.3.1   2019-03-13 [1] CRAN (R 3.6.3)                     
   remotes       2.1.1   2020-02-15 [1] CRAN (R 3.6.3)                     
   reprex        0.3.0   2019-05-16 [1] CRAN (R 3.6.3)                     
   rlang         0.4.5   2020-03-01 [1] CRAN (R 3.6.3)                     
   rprojroot     1.3-2   2018-01-03 [1] CRAN (R 3.6.3)                     
   rstudioapi    0.11    2020-02-07 [1] CRAN (R 3.6.3)                     
   rvest         0.3.5   2019-11-08 [1] CRAN (R 3.6.3)                     
   scales        1.1.0   2019-11-18 [1] CRAN (R 3.6.3)                     
   sessioninfo   1.1.1   2018-11-05 [1] CRAN (R 3.6.3)                     
   stringi       1.4.6   2020-02-17 [1] CRAN (R 3.6.2)                     
   stringr     * 1.4.0   2019-02-10 [1] CRAN (R 3.6.3)                     
   testthat      2.3.2   2020-03-02 [1] CRAN (R 3.6.3)                     
   tibble      * 3.0.1   2020-04-20 [1] CRAN (R 3.6.3)                     
   tidyr       * 1.0.2   2020-01-24 [1] CRAN (R 3.6.3)                     
   tidyselect    1.0.0   2020-01-27 [1] CRAN (R 3.6.3)                     
   tidyverse   * 1.3.0   2019-11-21 [1] CRAN (R 3.6.3)                     
 D TMB         * 1.7.16  2020-01-15 [1] CRAN (R 3.6.3)                     
   usethis     * 1.6.1   2020-04-29 [1] CRAN (R 3.6.3)                     
   utf8          1.1.4   2018-05-24 [1] CRAN (R 3.6.3)                     
   vctrs         0.2.4   2020-03-10 [1] CRAN (R 3.6.3)                     
   viridisLite   0.3.0   2018-02-01 [1] CRAN (R 3.6.3)                     
   withr         2.2.0   2020-04-20 [1] CRAN (R 3.6.3)                     
   xml2          1.3.2   2020-04-23 [1] CRAN (R 3.6.3)                     
   yaml          2.2.1   2020-02-01 [1] CRAN (R 3.6.3)                     

[1] C:/Users/jysullivan/Documents/R/win-library/3.6
[2] C:/Program Files/R/R-3.6.3/library
```
