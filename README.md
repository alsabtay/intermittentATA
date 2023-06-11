# intermittentATA
Intermittent demand time series analysis using the Ata Method based on the 'ATAforecasting package with the 'fable' Framework. 

This package provides a tidy R interface to Intermittent demand time series analysis using the Ata Method
procedure using [fabletools](https://github.com/tidyverts/fabletools). This
package makes use of the [ATAforecasting
package](https://cran.r-project.org/package=ATAforecasting) for R.

## Installation

You can install the **development** version from
[Github](https://github.com/alsabtay/intermittentATA) with:

``` r
# install.packages("remotes")
remotes::install_github("alsabtay/intermittentATA")
```

## Example

fmcgData: Intermittent sales data from a SKU form a store in Türkiye 2012--2019

``` r
library(intermittentATA)
as_tsibble(fmcgData) %>% model(crostonata = intermittentATA(value ~ d_trend(type = "M", parQ = 1) + i_trend("A") + intermittent("croston"))) %>% forecast(h=6)
``` 

## Links

[Github page](https://github.com/alsabtay/intermittentATA)

[Github.io page](https://alsabtay.github.io/intermittentATA/index.html)

[Github - Fable Modelling Wrappers for ATAforecasting Package](https://github.com/alsabtay/fable.ata)

[Github.io - Fable Modelling Wrappers for ATAforecasting Package](https://alsabtay.github.io/fable.ata/index.html)

[Github - ATAforecasting](https://github.com/alsabtay/ATAforecasting)

[Github.io - ATAforecasting](https://alsabtay.github.io/ATAforecasting/)

[Project team website](https://atamethod.wordpress.com/)


## License
This package is free and open source software, licensed under GPL-3.
