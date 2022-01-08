context("Histocompatibility")
library(histoc)

test_that("compatibility ABO", {
  expect_false(abo(cABO = 'A', dABO = 'O', iso = TRUE))
  expect_true(abo(cABO = 'A', dABO = 'O', iso = FALSE))
  expect_false(abo(cABO = 'AB', dABO = 'A', iso = TRUE))
  expect_true(abo(cABO = 'AB', dABO = 'A', iso = FALSE))
  expect_true(abo(cABO = 'AB', dABO = 'AB', iso = TRUE))
})

test_that("computes HLA mismatchs", {
  expect_equal(mmHLA(dA = c('01','02'),
                     dB = c('05','07'),
                     dDR = c('01','04'),
                     cA = c('01','02'),
                     cB = c('03','15'),
                     cDR = c('04','07'))[['mmA']], 0)
  expect_equal(mmHLA(dA = c('01','02'),
                     dB = c('05','07'),
                     dDR = c('01','04'),
                     cA = c('01','02'),
                     cB = c('03','15'),
                     cDR = c('04','07'))[['mmB']], 2)
  expect_equal(mmHLA(dA = c('01','02'),
                     dB = c('05','07'),
                     dDR = c('01','04'),
                     cA = c('01','02'),
                     cB = c('03','15'),
                     cDR = c('04','07'))[['mmDR']], 1)
  expect_equal(mmHLA(dA = c('01','02'),
                     dB = c('05','07'),
                     dDR = c('01','04'),
                     cA = c('01','02'),
                     cB = c('03','15'),
                     cDR = c('04','07'))[['mmHLA']], 3)
})

test_that("virtual crossmatch", {
  xmatch(dA = c('1','2'),
         dB = c('5','7'),
         dDR = c('1','4'),
         df.abs = cabs)$xm %>% .[1] %>%
    expect_equal("NEG")
  xmatch(dA = c('1','2'),
         dB = c('5','7'),
         dDR = c('1','4'),
         df.abs = cabs)$xm %>% .[6] %>%
    expect_equal("POS")
  xmatch(dA = c('1','2'),
         dB = c('5','7'),
         dDR = c('1','4'),
         df.abs = cabs)$xm %>% .[10] %>%
    expect_equal("POS")
})


test_that("Hiperimunized patients", {
  expect_true(hiper(cPRA = 99, cutoff = 85))
  expect_false(hiper(cPRA = 80, cutoff = 85))
  expect_true(hiper(cPRA = 90, cutoff = 80))
  expect_true(hiper(cPRA = 80, cutoff = 80))
  expect_false(hiper(cPRA = 1, cutoff = 85))
  expect_false(hiper(cPRA = 50, cutoff = 60))
})

test_that("Senior program classification", {
  expect_equal(sp(dage = 66, cage = 70), 1)
  expect_equal(sp(dage = 50, cage = 64), 2)
  expect_equal(sp(dage = 66, cage = 64), 3)
  expect_equal(sp(dage = 50, cage = 70), 3)
})

test_that("Tx Score (5 year survival probability)", {
  expect_equal(txscore(ageR = 40
                       , race = "White"
                       #, insurance = 0
                       , causeESRD = "Other"
                       , timeD = 60
                       , diabetesR = F
                       , coronary = F
                       , albumin = 1.5
                       , hemoglobin = 10
                       , ageD = 40
                       , diabetesD = "Absence"
                       , ECD = F
                       #, mmHLA = "0"
                       , mmHLA_A = 0
                       , mmHLA_B = 0
                       , mmHLA_DR = 0)$prob5y, 40)
  expect_equal(txscore(ageR = 72
                       , race = "Black"
                       #, insurance = 0
                       , causeESRD = "Hypertension"
                       , timeD = 21
                       , diabetesR = F
                       , coronary = F
                       , albumin = 2.7
                       , hemoglobin = 6.4
                       , ageD = 63
                       , diabetesD = "Unknown"
                       , ECD = T
                       #, mmHLA = "0"
                       , mmHLA_A = 0
                       , mmHLA_B = 1
                       , mmHLA_DR = 2)$prob5y, 57.25)
  expect_gte(txscore(ageR = 52
                       , race = "Black"
                       #, insurance = 0
                       , causeESRD = "Hypertension"
                       , timeD = 0
                       , diabetesR = F
                       , coronary = F
                       , albumin = 1
                       , hemoglobin = 10
                       , ageD = 70
                       , diabetesD = "Unknown"
                       , ECD = F
                       #, mmHLA = "0"
                       , mmHLA_A = 2
                       , mmHLA_B = 1
                       , mmHLA_DR = 0)$prob5y, 0)
  expect_lte(txscore(ageR = 52
                     , race = "Black"
                     #, insurance = 0
                     , causeESRD = "Hypertension"
                     , timeD = 0
                     , diabetesR = F
                     , coronary = F
                     , albumin = 1
                     , hemoglobin = 10
                     , ageD = 70
                     , diabetesD = "Unknown"
                     , ECD = F
                     #, mmHLA = "0"
                     , mmHLA_A = 2
                     , mmHLA_B = 1
                     , mmHLA_DR = 0)$prob5y, 100)
  })





