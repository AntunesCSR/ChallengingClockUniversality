# Installs needed to run the R code in this folder

install.packages('biocmanager')

BiocManager::install('clusterProfiler')
BiocManager::install('RDAVIDWebService', lib ="C:/Users/antun/R/win-library/4.3")  # not available for Bioconductor 3.18
BiocManager::install('BiocVersion',lib ="C:/Users/antun/R/win-library/4.3") 
BiocManager::install('topGO',lib ="C:/Users/antun/R/win-library/4.3")
BiocManager::install('org.Hs.eg.db',lib ="C:/Users/antun/R/win-library/4.3")
BiocManager::install('GOstats',lib ="C:/Users/antun/R/win-library/4.3")
BiocManager::install('AnnotationDbi',lib ="C:/Users/antun/R/win-library/4.3")
BiocManager::install('DOSE',lib ="C:/Users/antun/R/win-library/4.3") 
