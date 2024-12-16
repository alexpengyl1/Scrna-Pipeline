# scrna
2024-12-16
for function findallmarker， we can replace with presto:https://github.com/immunogenomics/presto


Presto performs a fast Wilcoxon rank sum test and auROC analysis. Latest benchmark ran 1 million observations, 1K features, and 10 groups in 16 seconds (sparse input) and 85 seconds (dense input).


for doulets-remove， we can use scDblFinder for better results:https://biostatsquid.com/scdblfinder-tutorial/
