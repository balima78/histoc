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


library(tictoc)
tic()
lima1() #data.table::merge.data.table --0.44
toc()


lima1()
