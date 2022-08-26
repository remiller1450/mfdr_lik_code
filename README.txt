Manuscript title: "Marginal false discovery rate control for likelihood-based penalized regression models"
Authors: Ryan Miller and Patrick Breheny 
Code development: Ryan Miller (contact: millerry@grinnell.edu)

R version 3.5.1 (2018-07-02)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] GEOquery_2.48.0          Biobase_2.40.0           BiocGenerics_0.26.0      ggplot2_3.0.0            gridExtra_2.3           
 [6] selectiveInference_1.2.4 intervals_0.15.1         covTest_1.02             MASS_7.3-50              glmpath_0.98            
[11] glmnet_2.0-16            foreach_1.4.4            lars_1.2                 survival_2.42-3          Matrix_1.2-14           
[16] ncvreg_3.11-0         

Directions:
Once all R packages are installed the .R code file may be ran as is.  The code will save .pdf files of each figure and table
in the manuscript.  The case study data are included as .RData files in the same directory as the code.

The second simulation takes many many hours to run (over 240 on my pc), the bottleneck is the lars.glm function which needs to be
run to perform covariance testing (a competitor to our method), it is unavoidable.  The first simulation and both case studies 
take only a few hours to run.