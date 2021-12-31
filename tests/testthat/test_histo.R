context("Histocompatibility")
library(histoc)

test_that("compatibility ABO", {
  expect_equal(abo(cABO = 'A', dABO = 'O', iso = TRUE), FALSE)
  expect_equal(abo(cABO = 'A', dABO = 'O', iso = FALSE), TRUE)
  expect_equal(abo(cABO = 'AB', dABO = 'A', iso = TRUE), FALSE)
  expect_equal(abo(cABO = 'AB', dABO = 'A', iso = FALSE), TRUE)
  expect_equal(abo(cABO = 'AB', dABO = 'AB', iso = TRUE), TRUE)
})

test_that("computes HLA mismatchs", {
  expect_equal(mmHLA(dA = c('01','02'),
                     dB = c('05','07'),
                     dDR = c('01','04'),
                     cA = c('01','02'),
                     cB = c('03','15'),
                     cDR = c('04','07'))$mmA, 0)
  expect_equal(mmHLA(dA = c('01','02'),
                     dB = c('05','07'),
                     dDR = c('01','04'),
                     cA = c('01','02'),
                     cB = c('03','15'),
                     cDR = c('04','07'))$mmB, 2)
  expect_equal(mmHLA(dA = c('01','02'),
                     dB = c('05','07'),
                     dDR = c('01','04'),
                     cA = c('01','02'),
                     cB = c('03','15'),
                     cDR = c('04','07'))$mmDR, 1)
  expect_equal(mmHLA(dA = c('01','02'),
                     dB = c('05','07'),
                     dDR = c('01','04'),
                     cA = c('01','02'),
                     cB = c('03','15'),
                     cDR = c('04','07'))$mmHLA, 3)
})

test_that("virtual crossmatch", {
  expect_equal(xmatch(dA = c('1','2'),
                      dB = c('5','7'),
                      dDR = c('1','4'),
                      df.abs = abs)$xm %>% .[1], "NEG")
  expect_equal(xmatch(dA = c('1','2'),
                      dB = c('5','7'),
                      dDR = c('1','4'),
                      df.abs = abs)$xm %>% .[6], "POS")
  expect_equal(xmatch(dA = c('1','2'),
                      dB = c('5','7'),
                      dDR = c('1','4'),
                      df.abs = abs)$xm %>% .[10], "POS")
})


test_that("Hiperimunized patients", {
  expect_equal(hiper(cPRA = 99, cutoff = 85), TRUE)
  expect_equal(hiper(cPRA = 80, cutoff = 85), FALSE)
  expect_equal(hiper(cPRA = 90, cutoff = 80), TRUE)
  expect_equal(hiper(cPRA = 80, cutoff = 80), TRUE)
  expect_equal(hiper(cPRA = 1, cutoff = 85), FALSE)
  expect_equal(hiper(cPRA = 50, cutoff = 60), FALSE)
})

test_that("Senior program classification", {
  expect_equal(sp(dage = 66, cage = 70), 1)
  expect_equal(sp(dage = 50, cage = 64), 2)
  expect_equal(sp(dage = 66, cage = 64), 3)
  expect_equal(sp(dage = 50, cage = 70), 3)
})



