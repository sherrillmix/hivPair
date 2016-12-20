# Code for "Resistance to type 1 interferons is a major determinant of HIV-1 transmission fitness"

## Necessary packages
The code uses several packages. If you'd like to install them all in one shot, you can do:
    install.packages(c('cluster','rstan','vioplot','png','vipor','ROCR','pROC'))
The code also uses that `parallel` package but that is included in base R and should not require installation. Versions used in the paper were:

Package|Version
-------|---------
cluster|2.0.4
rstan|2.13.2
vioplot|0.2
png|0.1-7
vipor|0.4.4
ROCR|1.0-7
pROC|1.8
R|3.3.1


## Generating all plots
To generate all plots, start R and use the command:
    source('runAll.R')
This will generate PCA, ROC, box and whisker and Bayesian plots in the `out` directory. 

