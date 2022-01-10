usethis::use_package("data.table")

usethis::use_testthat()

library(tidyverse)
library(data.table)

# Edit one or more files below R/.
# document() (if you’ve made any changes that impact help files or NAMESPACE)
# load_all()
# Run some examples interactively.
# test() (or test_file())
# check()

devtools::document()
devtools::load_all()
devtools::test()

donors.uk <- read.csv2("D:/2.PhD/HEADS/kars/files/donors.uk.csv")
#donors.uk$urgent <- 0
dim(donors.uk)

cabs <- abs

usethis::use_data(cabs)

# https://rich-iannone.github.io/pointblank/articles/VALID-III.html


library(tidyverse)
library(data.table)
lima1()

pts_PRA(cPRA = 79)

pts_HLA(
  dA = c('1','2'), dB = c('3','16'), dDR = c('4','7'),
  cA = c('1','2'), cB = c('3','15'), cDR = c('4','4')
)

mmHLA(dA = c('1','2'), dB = c('3','15'), dDR = c('7','7'),
      cA = c('1','2'), cB = c('3','15'), cDR = c('4','4'))

library(tictoc)
tic()
pt1() #data.table::merge.data.table --0.44
toc()


pt1() %>% view()

merge(candidates,
      xmatch(),
      all.x=TRUE) %>%
  #rowwise() %>%
  mutate(ptsdial = 0.1 * dialysis)


pts_HLA(mm.A = 0
        , mm.B = 1
        , mm.DR = 0)

a <- list('mmA' = 0,
     'mmB' = 0,
     'mmDR' = 0,
     'mmHLA' = 0)
a[['mmA']]
