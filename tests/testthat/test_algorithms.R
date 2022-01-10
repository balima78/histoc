context("Algorithms")
library(histoc)

test_that("Points cPRA PT's algorithm", {
  expect_equal(pts_PRA(cPRA = 0, pts.80 = 8, pts.50 = 4), 0)
  expect_equal(pts_PRA(cPRA = 81, pts.80 = 8, pts.50 = 4), 8)
  expect_equal(pts_PRA(cPRA = 51, pts.80 = 8, pts.50 = 4), 4)
  expect_equal(pts_PRA(cPRA = 90, pts.80 = 10, pts.50 = 4), 10)
  expect_equal(pts_PRA(cPRA = 90, pts.80 = 12, pts.50 = 4), 12)
  expect_equal(pts_PRA(cPRA = 90, pts.80 = 14, pts.50 = 4), 14)
  expect_equal(pts_PRA(cPRA = 40, pts.80 = 8, pts.50 = 4), 0)
  expect_equal(pts_PRA(cPRA = 60, pts.80 = 8, pts.50 = 6), 6)
  expect_equal(pts_PRA(cPRA = 60, pts.80 = 8, pts.50 = 8), 8)
  expect_equal(pts_PRA(cPRA = 60, pts.80 = 8, pts.50 = 10), 10)
  })

test_that("Points mmHLA PT's algorithm", {
  expect_equal(pts_HLA(itemA = 12, itemB = 8, itemC = 4, itemD = 2, itemE = 1,
                       mm.A = 0, mm.B = 1, mm.DR = 0), 4)
  expect_equal(pts_HLA(itemA = 12, itemB = 8, itemC = 4, itemD = 2, itemE = 1,
                       mm.A = 0, mm.B = 0, mm.DR = 0), 12)
  expect_equal(pts_HLA(itemA = 15, itemB = 8, itemC = 4, itemD = 2, itemE = 1,
                       mm.A = 0, mm.B = 0, mm.DR = 0), 15)
  expect_equal(pts_HLA(itemA = 12, itemB = 8, itemC = 4, itemD = 2, itemE = 1,
                       mm.A = 0, mm.B = 1, mm.DR = 1), 2)
  expect_equal(pts_HLA(itemA = 12, itemB = 8, itemC = 4, itemD = 2, itemE = 1,
                       mm.A = 1, mm.B = 0, mm.DR = 0), 8)
  expect_equal(pts_HLA(itemA = 12, itemB = 8, itemC = 4, itemD = 2, itemE = 1,
                       mm.A = 1, mm.B = 1, mm.DR = 1), 2)
  expect_equal(pts_HLA(itemA = 12, itemB = 8, itemC = 4, itemD = 2, itemE = 1,
                       mm.A = 0, mm.B = 2, mm.DR = 0), 1)

  })

test_that("Points age differences PT's algorithm", {
  expect_equal(pts_age(dage = 60, cage = 40, pts = 4), 4)
  expect_equal(pts_age(dage = 60, cage = 40, pts = 6), 6)
  expect_equal(pts_age(dage = 20, cage = 60, pts = 4), 0)
  expect_equal(pts_age(dage = 20, cage = 20, pts = 4), 4)
  expect_equal(pts_age(dage = 70, cage = 70, pts = 4), 4)
})
