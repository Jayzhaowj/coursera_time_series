#### download google trends from R
# install
install.packages("gtrendsR")
library(gtrendsR)

## nba and nfl
res <- gtrends(c("nfl", "nba"), geo = c("CA", "US"))
plot(res)

##
crypto_data <- gtrends("bitcoin")
plot(crypto_data)
