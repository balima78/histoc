library(tidyverse)

# Edit one or more files below R/.
# document() (if you’ve made any changes that impact help files or NAMESPACE)
# load_all()
# Run some examples interactively.
# test() (or test_file())
# check()

devtools::document()
devtools::load_all()
devtools::test()

donors <- read.csv2("D:/2.PhD/HEADS/kars/files/donors.csv")
candidates$urgent <- 0
dim(donors)

hlaDRet <- read.csv2("D:/2.PhD/HEADS/kars/files/hlaDRet.csv")
dim(hlaDRet)

usethis::use_data(donors)
