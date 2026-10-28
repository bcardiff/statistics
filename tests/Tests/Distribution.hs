{-# LANGUAGE FlexibleInstances, ScopedTypeVariables,
    ViewPatterns #-}
module Tests.Distribution (tests) where

import qualified Control.Exception as E
import Data.List (find)
import Data.Typeable (Typeable)
import Data.Word
import Numeric.MathFunctions.Constants (m_tiny,m_huge,m_epsilon)
import Numeric.MathFunctions.Comparison
import Statistics.Distribution
import Statistics.Distribution.Beta           (BetaDistribution)
import Statistics.Distribution.Binomial       (BinomialDistribution)
import Statistics.Distribution.CauchyLorentz
import Statistics.Distribution.ChiSquared     (ChiSquared)
import Statistics.Distribution.Exponential    (ExponentialDistribution,exponential)
import Statistics.Distribution.FDistribution  (FDistribution,fDistribution)
import Statistics.Distribution.Gamma          (GammaDistribution,gammaDistr)
import Statistics.Distribution.Geometric
import Statistics.Distribution.Hypergeometric
import Statistics.Distribution.Laplace        (LaplaceDistribution)
import Statistics.Distribution.Lognormal      (LognormalDistribution)
import Statistics.Distribution.NegativeBinomial (NegativeBinomialDistribution)
import Statistics.Distribution.Normal         (NormalDistribution)
import Statistics.Distribution.Poisson        (PoissonDistribution)
import Statistics.Distribution.StudentT
import Statistics.Distribution.Transform      (LinearTransform(..))
import Statistics.Distribution.Uniform        (UniformDistribution)
import Statistics.Distribution.Weibull        (WeibullDistribution)
import Statistics.Distribution.DiscreteUniform (DiscreteUniform)
import Test.Tasty                 (TestTree, testGroup)
import Test.Tasty.QuickCheck      (testProperty)
import Test.Tasty.ExpectedFailure (ignoreTest)
import Test.QuickCheck as QC
import Test.QuickCheck.Monadic as QC
import Text.Printf (printf)

import Tests.ApproxEq  (ApproxEq(..))
import Tests.ExactDistribution (exactDistributionTests)
import Tests.Helpers   (T(..), Double01(..), testAssertion, typeName)
import Tests.Helpers   (monotonicallyIncreasesIEEE,isDenorm)
import Tests.Orphanage ()

-- | Tests for all distributions
tests :: TestTree
tests = testGroup "Tests for all distributions"
  [ contDistrTests (T :: T BetaDistribution        )
  , contDistrTests (T :: T CauchyDistribution      )
  , contDistrTests (T :: T ChiSquared              )
  , contDistrTests (T :: T ExponentialDistribution )
  , contDistrTests (T :: T GammaDistribution       )
  , contDistrTests (T :: T LaplaceDistribution     )
  , contDistrTests (T :: T LognormalDistribution   )
  , contDistrTests (T :: T NormalDistribution      )
  , contDistrTests (T :: T UniformDistribution     )
  , contDistrTests (T :: T WeibullDistribution     )
  , contDistrTests (T :: T StudentT                )
  , contDistrTests (T :: T (LinearTransform NormalDistribution))
  , contDistrTests (T :: T FDistribution           )

  , discreteDistrTests (T :: T BinomialDistribution       )
  , discreteDistrTests (T :: T GeometricDistribution      )
  , discreteDistrTests (T :: T GeometricDistribution0     )
  , discreteDistrTests (T :: T HypergeometricDistribution )
  , discreteDistrTests (T :: T NegativeBinomialDistribution )
  , discreteDistrTests (T :: T PoissonDistribution        )
  , discreteDistrTests (T :: T DiscreteUniform            )

  , exactDistributionTests
  , linearTransformTests
  , unitTests
  ]

----------------------------------------------------------------
-- Tests
----------------------------------------------------------------

-- Tests for continuous distribution
contDistrTests :: (Param d, ContDistr d, QC.Arbitrary d, Typeable d, Show d) => T d -> TestTree
contDistrTests t = testGroup ("Tests for: " ++ typeName t) $
  cdfTests t ++
  [ testProperty "PDF sanity"              $ pdfSanityCheck     t
  , (if quantileIsInvCDF_enabled t then id else ignoreTest)
  $ testProperty "Quantile is CDF inverse" $ quantileIsInvCDF t
  , testProperty "quantile fails p<0||p>1" $ quantileShouldFail t
  , testProperty "log density check"       $ logDensityCheck    t
  , testProperty "complQuantile"           $ complQuantileCheck t
  ]

-- Tests for discrete distribution
discreteDistrTests :: (Param d, DiscreteDistr d, QC.Arbitrary d, Typeable d, Show d) => T d -> TestTree
discreteDistrTests t = testGroup ("Tests for: " ++ typeName t) $
  cdfTests t ++
  [ testProperty "Prob. sanity"         $ probSanityCheck       t
  , testProperty "CDF is sum of prob."  $ discreteCDFcorrect    t
  , testProperty "Discrete CDF is OK"   $ cdfDiscreteIsCorrect  t
  , testProperty "log probability check" $ logProbabilityCheck   t
  ]

-- Tests for distributions which have CDF
cdfTests :: (Param d, Distribution d, QC.Arbitrary d, Show d) => T d -> [TestTree]
cdfTests t =
  [ testProperty "C.D.F. sanity"        $ cdfSanityCheck         t
  , testProperty "CDF limit at +inf"    $ cdfLimitAtPosInfinity  t
  , (if cdfLimitAtNegInfinity_enabled t then id else ignoreTest)
  $ testProperty "CDF limit at -inf"    $ cdfLimitAtNegInfinity  t
  , testProperty "CDF at +inf = 1"      $ cdfAtPosInfinity       t
  , testProperty "CDF at -inf = 1"      $ cdfAtNegInfinity       t
  , testProperty "CDF is nondecreasing" $ cdfIsNondecreasing     t
  , testProperty "1-CDF is correct"     $ cdfComplementIsCorrect t
  ]


----------------------------------------------------------------

-- CDF is in [0,1] range
cdfSanityCheck :: (Distribution d) => T d -> d -> Double -> Bool
cdfSanityCheck _ d x = c >= 0 && c <= 1
  where c = cumulative d x

-- CDF never decreases
cdfIsNondecreasing :: (Distribution d) => T d -> d -> Double -> Double -> Bool
cdfIsNondecreasing _ d = monotonicallyIncreasesIEEE $ cumulative d

-- cumulative d +∞ = 1
cdfAtPosInfinity :: (Distribution d) => T d -> d -> Bool
cdfAtPosInfinity _ d
  = cumulative d (1/0) == 1

-- cumulative d - ∞ = 0
cdfAtNegInfinity :: (Distribution d) => T d -> d -> Bool
cdfAtNegInfinity _ d
  = cumulative d (-1/0) == 0

-- CDF limit at +∞ is 1
cdfLimitAtPosInfinity :: (Param d, Distribution d) => T d -> d -> Bool
cdfLimitAtPosInfinity _ d
  = Just 1.0 == find (>=1) probs
  where
    probs = map (cumulative d)
          $ takeWhile (< (m_huge/2))
          $ iterate (*1.4) 1

-- CDF limit at -∞ is 0
cdfLimitAtNegInfinity :: (Param d, Distribution d) => T d -> d -> Bool
cdfLimitAtNegInfinity _ d
  = Just 0 == find (<=0) probs
  where
    probs = map (cumulative d)
          $ takeWhile (> (-m_huge/2))
          $ iterate (*1.4) (-1)


-- CDF's complement is implemented correctly
cdfComplementIsCorrect :: (Distribution d, Param d) => T d -> d -> Double -> Property
cdfComplementIsCorrect _ d x
  = counterexample ("err. tolerance = " ++ show tol)
  $ counterexample ("difference     = " ++ show delta)
  $ delta <= tol
  where
    tol   = prec_complementCDF d
    delta = 1 - (cumulative d x + complCumulative d x)

-- CDF for discrete distribution uses <= for comparison
cdfDiscreteIsCorrect :: (Param d, DiscreteDistr d) => T d -> d -> Property
cdfDiscreteIsCorrect _ d
  = counterexample (unlines badN)
  $ null badN
  where
    -- We are checking that:
    --
    -- > CDF(i) - CDF(i-e) = P(i)
    --
    -- Approximate equality is tricky here. Scale is set by maximum
    -- value of CDF and probability. Case when all probabilities are
    -- zero should be treated specially.
    badN = [ printf "N=%3i    p[i]=%g\tp[i+1]=%g\tdP=%g\trelerr=%g" i p p1 dp ((p1-p-dp) / max p1 dp)
           | i <- [0 .. 100]
           , let p      = cumulative d $ fromIntegral i - 1e-6
                 p1     = cumulative d $ fromIntegral i
                 dp     = probability d i
                 relerr = ((p1 - p) - dp) / max p1 dp
           , p  > m_tiny || p == 0
           , p1 > m_tiny
           , dp > m_tiny
           , relerr > tol
           ]
    tol = prec_discreteCDF d

logDensityCheck :: (Param d, ContDistr d) => T d -> d -> Double -> Property
logDensityCheck _ d x
  = not (isDenorm x)
  ==> ( counterexample (printf "density    = %g" p)
      $ counterexample (printf "logDensity = %g" logP)
      $ counterexample (printf "log p      = %g" (log p))
      $ counterexample (printf "ulps[log]  = %i" ulpsLog)
      $ counterexample (printf "ulps[lin]  = %i" ulpsLin)
      $ or [ p == 0      && logP == (-1/0)
           , p <= m_tiny && logP < log m_tiny
             -- To avoid problems with roundtripping error in case
             -- when density is computed as exponent of logDensity we
             -- accept either inequality
           ,  (ulpsLog <= n) || (ulpsLin <= n)
           ])
  where
    p       = density d x
    logP    = logDensity d x
    n       = prec_logDensity d
    ulpsLog = ulpDistance (log p) logP
    ulpsLin = ulpDistance p       (exp logP)

-- PDF is positive
pdfSanityCheck :: (ContDistr d) => T d -> d -> Double -> Bool
pdfSanityCheck _ d x = p >= 0
  where p = density d x

complQuantileCheck :: (ContDistr d) => T d -> d -> Double01 -> Property
complQuantileCheck _ d (Double01 p)
  = counterexample (printf "x0 = %g" x0)
  $ counterexample (printf "x1 = %g" x1)
  $ counterexample (printf "abs err = %g" $ abs (x1 - x0))
  $ counterexample (printf "rel err = %g" $ relativeError x1 x0)
  -- We avoid extreme tails of distributions
  --
  -- FIXME: all parameters are arbitrary at the moment
  $ and [ p > 0.01
        , p < 0.99
        , not $ isInfinite x0
        , not $ isInfinite x1
        ] ==> (if x0 < 1e6 then abs (x1 - x0) < 1e-6 else relativeError x1 x0 < 1e-12)
  where
    x0 = quantile      d (1 - p)
    x1 = complQuantile d p

-- Quantile is inverse of CDF
quantileIsInvCDF :: (Param d, ContDistr d) => T d -> d -> Double01 -> Property
quantileIsInvCDF _ d (Double01 p) =
  and [ p > m_tiny
      , p < 1
      , x > m_tiny
      , dens > 0
      ] ==>
    ( counterexample (printf "Quantile      = %g" x )
    $ counterexample (printf "Probability   = %g" p )
    $ counterexample (printf "Probability'  = %g" p')
    $ counterexample (printf "Rel. error    = %g" (relativeError p p'))
    $ counterexample (printf "Abs. error    = %e" (abs $ p - p'))
    $ counterexample (printf "Expected err. = %g" err)
    $ counterexample (printf "Distance      = %i" (ulpDistance p p'))
    $ counterexample (printf "Err/est       = %g" (fromIntegral (ulpDistance p p') / err))
    $ ulpDistance p p' <= round err
    )
  where
    -- Algorithm for error estimation is taken from here
    --
    -- http://sepulcarium.org/posts/2012-07-19-rounding_effect_on_inverse.html
    dens = density    d x
    err  = eps + eps' * abs (x / p) * dens
    --
    x    = quantile   d p
    p'   = cumulative d x
    (eps,eps') = prec_quantile_CDF d

-- Test that quantile fails if p<0 or p>1
quantileShouldFail :: (ContDistr d) => T d -> d -> Double -> Property
quantileShouldFail _ d p =
  p < 0 || p > 1 ==> QC.monadicIO $ do r <- QC.run $ E.catch
                                              (False <$ (return $! quantile d p))
                                              (\(_ :: E.SomeException) -> return True)
                                       QC.assert r


-- Probability is in [0,1] range
probSanityCheck :: (DiscreteDistr d) => T d -> d -> Int -> Bool
probSanityCheck _ d x = p >= 0 && p <= 1
  where p = probability d x

-- Check that discrete CDF is correct
discreteCDFcorrect :: (DiscreteDistr d) => T d -> d -> Int -> Int -> Property
discreteCDFcorrect _ d a b
  = counterexample (printf "CDF   = %g" p1)
  $ counterexample (printf "Sum   = %g" p2)
  $ counterexample (printf "Delta = %g" (abs (p1 - p2)))
  $ abs (p1 - p2) < 3e-10
  -- Avoid too large differences. Otherwise there is to much to sum
  --
  -- Absolute difference is used guard against precision loss when
  -- close values of CDF are subtracted
  where
    n  = min a b
    m  = n + (abs (a - b) `mod` 100)
    p1 = cumulative d (fromIntegral m + 0.5) - cumulative d (fromIntegral n - 0.5)
    p2 = sum $ map (probability d) [n .. m]

logProbabilityCheck :: (Param d, DiscreteDistr d) => T d -> d -> Int -> Property
logProbabilityCheck _ d x
  = counterexample (printf "probability    = %g" p)
  $ counterexample (printf "logProbability = %g" logP)
  $ counterexample (printf "log p          = %g" (log p))
  $ counterexample (printf "ulps[log]      = %i" ulpsLog)
  $ counterexample (printf "ulps[lin]      = %i" ulpsLin)
  $ or [ p == 0     && logP == (-1/0)
       , p < 1e-308 && logP < 609
         -- To avoid problems with roundtripping error in case
         -- when density is computed as exponent of logDensity we
         -- accept either inequality
       ,  (ulpsLog <= n) || (ulpsLin <= n)
       ]
  where
    p    = probability d x
    logP = logProbability d x
    n    = prec_logDensity d
    ulpsLog = ulpDistance (log p) logP
    ulpsLin = ulpDistance p       (exp logP)


-- | Parameters for distribution testing. Some distribution require
--   relaxing parameters a bit
class Param a where
  -- | Whether quantileIsInvCDF is enabled
  quantileIsInvCDF_enabled :: T a -> Bool
  quantileIsInvCDF_enabled _ = True
  -- | Whether cdfLimitAtNegInfinity is enabled
  cdfLimitAtNegInfinity_enabled :: T a -> Bool
  cdfLimitAtNegInfinity_enabled _ = True
  -- | Precision for 'quantileIsInvCDF' test
  prec_quantile_CDF :: a -> (Double,Double)
  prec_quantile_CDF _ = (16,16)
  -- |
  prec_discreteCDF :: a -> Double
  prec_discreteCDF _ = 32 * m_epsilon
  -- | Precision of CDF's complement
  prec_complementCDF :: a -> Double
  prec_complementCDF _ = 1e-14
  -- | Precision for logDensity check
  prec_logDensity :: a -> Word64
  prec_logDensity _ = 32

instance Param StudentT where
  -- FIXME: disabled unless incompleteBeta troubles are sorted out
  quantileIsInvCDF_enabled _ = False

instance Param BetaDistribution where
  -- FIXME: See https://github.com/haskell/statistics/issues/161 for details
  quantileIsInvCDF_enabled _ = False

instance Param FDistribution where
  -- FIXME: disabled unless incompleteBeta troubles are sorted out
  quantileIsInvCDF_enabled _ = False
  -- We compute CDF and complement using same method so precision
  -- should be very good here.
  prec_complementCDF _ = 64 * m_epsilon

instance Param ChiSquared where
  prec_quantile_CDF _ = (32,32)

instance Param BinomialDistribution where
  prec_discreteCDF _ = 1e-12
  prec_logDensity  _ = 48
instance Param CauchyDistribution where
  -- Distribution is long-tailed enough that we may never get to zero
  cdfLimitAtNegInfinity_enabled _ = False

instance Param DiscreteUniform
instance Param ExponentialDistribution
instance Param GammaDistribution where
  -- We lose precision near `incompleteGamma 10` because of error
  -- introduced by exp . logGamma.  This could only be fixed in
  -- math-function by implementing gamma
  prec_quantile_CDF _ = (24,24)
  prec_logDensity   _ = 512
instance Param GeometricDistribution
instance Param GeometricDistribution0
instance Param HypergeometricDistribution
instance Param LaplaceDistribution
instance Param LognormalDistribution where
  prec_quantile_CDF _ = (64,64)
instance Param NegativeBinomialDistribution where
  prec_discreteCDF  _ = 1e-12
  prec_logDensity   _ = 48
instance Param NormalDistribution
instance Param PoissonDistribution
instance Param UniformDistribution
instance Param WeibullDistribution
instance Param a => Param (LinearTransform a)

----------------------------------------------------------------
-- Unit tests
----------------------------------------------------------------

unitTests :: TestTree
unitTests = testGroup "Unit tests"
  [ testAssertion "density (gammaDistr 150 1/150) 1 == 4.883311" $
      4.883311418525483 =~ density (gammaDistr 150 (1/150)) 1
    -- Student-T
  , testStudentPDF 0.3  1.34  0.0648215  -- PDF
  , testStudentPDF 1    0.42  0.27058
  , testStudentPDF 4.4  0.33  0.352994
  , testStudentCDF 0.3  3.34  0.757146   -- CDF
  , testStudentCDF 1    0.42  0.626569
  , testStudentCDF 4.4  0.33  0.621739
    -- Student-T General
  , testStudentUnstandardizedPDF 0.3    1.2  4      0.45 0.0533456  -- PDF
  , testStudentUnstandardizedPDF 4.3  (-2.4) 3.22 (-0.6) 0.0971141
  , testStudentUnstandardizedPDF 3.8    0.22 7.62   0.14 0.0490523
  , testStudentUnstandardizedCDF 0.3    1.2  4      0.45 0.458035   -- CDF
  , testStudentUnstandardizedCDF 4.3  (-2.4) 3.22 (-0.6) 0.698001
  , testStudentUnstandardizedCDF 3.8    0.22 7.62   0.14 0.496076
    -- F-distribution
  , testFdistrPDF  1  3   3     (1/(6 * pi)) -- PDF
  , testFdistrPDF  2  2   1.2   0.206612
  , testFdistrPDF  10 12  8     0.000385613179281892790166
  , testFdistrCDF  1  3   3     0.81830988618379067153 -- CDF
  , testFdistrCDF  2  2   1.2   0.545455
  , testFdistrCDF  10 12  8     0.99935509863451408041
  ]
  where
    -- Student-T
    testStudentPDF ndf x exact
      = testAssertion (printf "density (studentT %f) %f ~ %f" ndf x exact)
      $ eq 1e-5  exact  (density (studentT ndf) x)
    testStudentCDF ndf x exact
      = testAssertion (printf "cumulative (studentT %f) %f ~ %f" ndf x exact)
      $ eq 1e-5  exact  (cumulative (studentT ndf) x)
    -- Student-T General
    testStudentUnstandardizedPDF ndf mu sigma x exact
      = testAssertion (printf "density (studentTUnstandardized %f %f %f) %f ~ %f" ndf mu sigma x exact)
      $ eq 1e-5  exact  (density (studentTUnstandardized ndf mu sigma) x)
    testStudentUnstandardizedCDF ndf mu sigma x exact
      = testAssertion (printf "cumulative (studentTUnstandardized %f %f %f) %f ~ %f" ndf mu sigma x exact)
      $ eq 1e-5  exact  (cumulative (studentTUnstandardized ndf mu sigma) x)
    -- F-distribution
    testFdistrPDF n m x exact
      = testAssertion (printf "density (fDistribution %i %i) %f ~ %f [got %f]" n m x exact d)
      $ eq 1e-5  exact d
      where d = density (fDistribution n m) x
    testFdistrCDF n m x exact
      = testAssertion (printf "cumulative (fDistribution %i %i) %f ~ %f [got %f]" n m x exact d)
      $ eq 1e-5  exact d
      where d = cumulative (fDistribution n m) x

----------------------------------------------------------------
-- Tests for linear transformation of distribution
----------------------------------------------------------------

-- Unit tests for LinearTransform. Exponential distribution is used as
-- a base since all its moments are known analytically, it's not
-- symmetric and it provides accurate complCumulative.
--
-- > LinearTransform loc sc d  ~  loc + sc·X   where X ~ d
linearTransformTests :: TestTree
linearTransformTests = testGroup "Linear transform"
  [ -- Mean is affected by scale: E[loc + sc·X] = loc + sc·E[X]
    testAssertion "mean (LinearTransform 0 3 e) == 3"
  $ 3 == mean (scaled 3)
  , testAssertion "maybeMean (LinearTransform 0 3 e) == Just 3"
  $ Just 3 == maybeMean (scaled 3)
    -- Differential entropy is not scale invariant: h(loc + sc·X) = h(X) + log|sc|
  , testAssertion "entropy (LinearTransform 0 3 e) ~ 1 + log 3"
  $ eq 1e-12 (1 + log 3) (entropy (scaled 3))
  , testAssertion "maybeEntropy (LinearTransform 0 3 e) ~ Just (1 + log 3)"
  $ maybe False (eq 1e-12 (1 + log 3)) (maybeEntropy (scaled 3))
    -- Negative scale reverses real line and swaps tails:
    -- P(loc + sc·X ≤ x) = P(X ≥ (x-loc)/sc)
  , testAssertion "cumulative (LinearTransform 0 (-1) e) (-1) ~ exp(-1)"
  $ eq 1e-12 (exp (-1)) (cumulative (scaled (-1)) (-1))
  , testAssertion "complCumulative (LinearTransform 0 (-1) e) (-1) ~ 1 - exp(-1)"
  $ eq 1e-12 (1 - exp (-1)) (complCumulative (scaled (-1)) (-1))
  , testAssertion "quantile (LinearTransform 0 (-1) e) 0.25 ~ -log 4"
  $ eq 1e-12 (negate (log 4)) (quantile (scaled (-1)) 0.25)
  , testAssertion "complQuantile (LinearTransform 0 (-1) e) 0.25 ~ -log (4/3)"
  $ eq 1e-12 (negate (log (4/3))) (complQuantile (scaled (-1)) 0.25)
  , testAssertion "density (LinearTransform 0 (-1) e) (-1) ~ exp(-1)"
  $ eq 1e-12 (exp (-1)) (density (scaled (-1)) (-1))
  , testAssertion "logDensity (LinearTransform 0 (-1) e) (-1) ~ -1"
  $ eq 1e-12 (-1) (logDensity (scaled (-1)) (-1))
  , testAssertion "mean (LinearTransform 0 (-1) e) == -1"
  $ -1 == mean (scaled (-1))
  , testAssertion "stdDev (LinearTransform 0 (-1) e) == 1"
  $ 1 == stdDev (scaled (-1))
  , testAssertion "maybeStdDev (LinearTransform 0 (-1) e) == Just 1"
  $ Just 1 == maybeStdDev (scaled (-1))
    -- Complement of CDF of underlying distribution should be used in
    -- order to not lose precision in far tail
  , testAssertion "complCumulative (LinearTransform 2 1 e) 42 == complCumulative e 40"
  $ complCumulative e 40 == complCumulative (shifted 2) 42
    -- Transformations which are already correct and should stay so
  , testAssertion "cumulative (LinearTransform 2 1 e) 3 ~ 1 - exp(-1)"
  $ eq 1e-12 (1 - exp (-1)) (cumulative (shifted 2) 3)
  , testAssertion "density (LinearTransform 0 3 e) 3 == density e 1 / 3"
  $ density e 1 / 3 == density (scaled 3) 3
  , testAssertion "variance (LinearTransform 0 3 e) == 9"
  $ 9 == variance (scaled 3)
  , testAssertion "stdDev (LinearTransform 0 3 e) == 3"
  $ 3 == stdDev (scaled 3)
  , testAssertion "mean (LinearTransform 2 1 e) == 3"
  $ 3 == mean (shifted 2)
  ]
  where
    e           = exponential 1
    scaled  sc  = LinearTransform 0   sc e
    shifted loc = LinearTransform loc 1  e
