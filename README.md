# R package: "misclassGLM" #
## Description ##
Estimates models that extend the standard GLM to take misclassification into account. The models require side information from a secondary data set on the misclassification process, i.e. some sort of misclassification probabilities conditional on some common covariates. A detailed description of the algorithm can be found in Dlugosz, Mammen and Wilke (2017) Computational Statistics & Data Analysis 110:145-159 [http://dx.doi.org/10.1016/j.csda.2017.01.003](http://dx.doi.org/10.1016/j.csda.2017.01.003).

## Installation ##
### From CRAN ###
The easiest way to use any of the functions in the misclassGLM package is to install the CRAN version. It can be installed from within R using the command:

```
#!R

install.packages("misclassGLM")
```


### From bitbucket ###
The devtools package contains functions that allow you to install R packages directly from bitbucket or github. If you've installed and loaded the devtools package, the installation command is

```
#!R

install_bitbucket("sdlugosz/misclassGLM")
```
