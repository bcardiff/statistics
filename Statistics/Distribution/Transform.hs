{-# LANGUAGE DeriveDataTypeable, DeriveGeneric, FlexibleContexts,
    FlexibleInstances, UndecidableInstances #-}
-- |
-- Module    : Statistics.Distribution.Transform
-- Copyright : (c) 2013 John McDonnell;
-- License   : BSD3
--
-- Maintainer  : bos@serpentine.com
-- Stability   : experimental
-- Portability : portable
--
-- Transformations over distributions
module Statistics.Distribution.Transform
    (
      LinearTransform (..)
    , linTransFixedPoint
    , scaleAround
    ) where

import Data.Aeson (FromJSON, ToJSON)
import Data.Binary (Binary)
import Data.Binary (put, get)
import Data.Data (Data, Typeable)
import GHC.Generics (Generic)
import qualified Statistics.Distribution as D

-- | Linear transformation applied to distribution.
--
-- > LinearTransform μ σ _
-- > x' = μ + σ·x
data LinearTransform d = LinearTransform
  { linTransLocation :: {-# UNPACK #-} !Double
    -- ^ Location parameter.
  , linTransScale    :: {-# UNPACK #-} !Double
    -- ^ Scale parameter.
  , linTransDistr    :: d
    -- ^ Distribution being transformed.
  } deriving (Eq, Show, Read, Typeable, Data, Generic)

instance (FromJSON d) => FromJSON (LinearTransform d)
instance (ToJSON d) => ToJSON (LinearTransform d)

instance (Binary d) => Binary (LinearTransform d) where
    get = LinearTransform <$> get <*> get <*> get
    put (LinearTransform x y z) = put x >> put y >> put z

-- | Apply linear transformation to distribution.
scaleAround :: Double           -- ^ Fixed point
            -> Double           -- ^ Scale parameter
            -> d                -- ^ Distribution
            -> LinearTransform d
scaleAround x0 sc = LinearTransform (x0 * (1 - sc)) sc

-- | Get fixed point of linear transformation
linTransFixedPoint :: LinearTransform d -> Double
linTransFixedPoint (LinearTransform loc sc _) = loc / (1 - sc)

instance Functor LinearTransform where
  fmap f (LinearTransform loc sc dist) = LinearTransform loc sc (f dist)

-- Note that negative scale reverses order of real line, so mass below
-- a point becomes mass above its preimage:
--
-- > P(loc + sc·X ≤ x) = P(X ≥ (x-loc)/sc)   when sc < 0
--
-- Thus complement of CDF of underlying distribution has to be used in
-- that case and same for quantiles.
instance D.Distribution d => D.Distribution (LinearTransform d) where
  cumulative (LinearTransform loc sc dist) x
    | sc < 0    = D.complCumulative dist y
    | otherwise = D.cumulative      dist y
    where y = (x - loc) / sc
  complCumulative (LinearTransform loc sc dist) x
    | sc < 0    = D.cumulative      dist y
    | otherwise = D.complCumulative dist y
    where y = (x - loc) / sc

instance D.ContDistr d => D.ContDistr (LinearTransform d) where
  density    (LinearTransform loc sc dist) x = D.density    dist ((x-loc) / sc) / abs sc
  logDensity (LinearTransform loc sc dist) x = D.logDensity dist ((x-loc) / sc) - log (abs sc)
  quantile      (LinearTransform loc sc dist) p
    | sc < 0    = loc + sc * D.complQuantile dist p
    | otherwise = loc + sc * D.quantile      dist p
  complQuantile (LinearTransform loc sc dist) p
    | sc < 0    = loc + sc * D.quantile      dist p
    | otherwise = loc + sc * D.complQuantile dist p

instance D.MaybeMean d => D.MaybeMean (LinearTransform d) where
  maybeMean (LinearTransform loc sc dist) = (\m -> loc + sc * m) <$> D.maybeMean dist

instance (D.Mean d) => D.Mean (LinearTransform d) where
  mean (LinearTransform loc sc dist) = loc + sc * D.mean dist

instance D.MaybeVariance  d => D.MaybeVariance (LinearTransform d) where
  maybeVariance (LinearTransform _ sc dist) = (*(sc*sc)) <$> D.maybeVariance dist
  maybeStdDev   (LinearTransform _ sc dist) = (*abs sc)  <$> D.maybeStdDev dist

instance (D.Variance d) => D.Variance (LinearTransform d) where
  variance (LinearTransform _ sc dist) = sc * sc * D.variance dist
  stdDev   (LinearTransform _ sc dist) = abs sc * D.stdDev dist

instance (D.MaybeEntropy d) => D.MaybeEntropy (LinearTransform d) where
  maybeEntropy (LinearTransform _ sc dist) = (+ log (abs sc)) <$> D.maybeEntropy dist

instance (D.Entropy d) => D.Entropy (LinearTransform d) where
  entropy (LinearTransform _ sc dist) = D.entropy dist + log (abs sc)

instance D.ContGen d => D.ContGen (LinearTransform d) where
  genContVar (LinearTransform loc sc d) g = do
    x <- D.genContVar d g
    return $! loc + sc * x
