# [配置bioconda](https://bioconda.github.io/index.html#set-up-channels)

## 1.Install conda

## 2.Set up channels
The conda-forge channel contains many general-purpose packages not already found in the defaults channel.

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    
## 3.Install packages
Browse the packages to see what’s available.
Bioconda is now enabled, so any packages on the bioconda channel can be installed into the current conda environment:

    conda install bwa
