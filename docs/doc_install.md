<!-- 这个文档用来说明如何安装流程, 能让它运行起来.  -->

# Installation 

This pipeline does not need to be compiled, but it will only run successfully with the help of depended packages

## Fresh Installation

If [git](https://git-scm.com) is installed, you can download this pipeline using following command.

<!-- todo 修改为最终路径 -->
``` bash
git clone 
```

<!-- todo 下载的路径 -->
You can also download the package as a zip file and unzip it manually.

Also [Python 3](https://www.python.org) and [R](https://www.r-project.org) need to be installed before you can proceed with the following steps.

## install dependencies packages: python

A versatile helper script `info.py` is provided by this pipeline . 

When this script is invoked in command line, it first checks whether the software-related python packages are installed properly. 

It will report an error once any python package is not installed. users can check the error message to find the missing package.

## install dependencies: R

These R package may be required for analysis: `tidyverse`, `Biostrings` , `rtracklayer`, `GenomicFeatures`, the latter 3 packages are part of  [bioconductor](https://www.bioconductor.org).

use the following code to install them in the R console:

To install `tidyverse`:

```R
install.packages("tidyverse")
```

To install other 3 bioconductor packages

```R
install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("rtracklayer")
BiocManager::install("GenomicFeatures")

```

If you need to install  `Biostrings` `rtracklayer` `GenomicFeatures` in an older version of R that does not work with `BiocManager`, please check bioconductor page of these packages.

## modify the default configuration template

The default configure template lies in the `pysrc/body/default.cpysrc/body/default.cfg`. Each time you invoke `new_config.py`, you're actually making a copy of this file.

This template is almost empty after a fresh installation. 
You can modify this file to save your time in future use.

More details can be found in `conf.md`
