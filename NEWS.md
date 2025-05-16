## version 2.4.3
- Adds conf.level argument to the summary function to change the confidence level of the confidence intervals.
- Updates documentation and references.

## version 2.4.2
- Fixes https://github.com/FBartos/zcurve/issues/22

## version 2.4.1
- Fixes ODR computation with censored data

## version 2.4.0
- Implementation of ggploting function.

## version 2.3.0
- Implementation of parallel bootstrap.

## version 2.2.0
- Implementation of z-curve version that deals with clustered estimates. 
- Fixing bootstrap of censored data -- now both the censored and uncensored data are bootstrapped.

## version 2.1.2
- Fixing number of statistically significant studies displayed in the annotation of z-curve plot.

## version 2.1.1
- Allowing to pass the xlim and ylim arguments to the z-curve plot function.

## version 2.1.0
- Adding 'zcurve_data' function to transform different tests statistics into a z-curve data object.

## version 2.0.0
- Adding the possibility to fit z-curve to censored/rounded input.

## version 1.0.9
- Changing the ERR estimates to be take the directionality of the expected replication into account.

## version 1.0.8
- Fixing an input test of z-scores (checking whether p-values weren't used by a mistake) being triggered when all z-scores were negative.

## version 1.0.7
- Fixing incorrect passing of the significance level when the lower fitting range is not equal to alpha.

## version 1.0.6
- Fixing the formatting of the NEWS.md file.

## version 1.0.5
- Adding convenience functions for extracting the z-curve estimates directly from the fitted object. (issue #3)

## version 1.0.4
- Adding observed discovery rate, number of included studies, and number of significant studies to the summary output. (issue #3)

## version 1.0.3
- Fixing a help file entry and making sure of proper control specification '[[]]' vs '$' list accessing.

## version 1.0.2
- Fixing a bug when producing figure with no non-significant z-scores.

## version 1.0.1
- Added check for supplying p-values as z-scores. (issue #2)

## version 1.0.
- The release version of the package.