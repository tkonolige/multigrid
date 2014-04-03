{-# LANGUAGE BangPatterns, FlexibleContexts, ScopedTypeVariables #-}

-- a couple of tests using the poisson problem with a 2d stencil

import Prelude hiding (replicate)
import Math.Solver.Multigrid
import Data.Array.Repa
import Data.Array.Repa.Stencil
import Data.Array.Repa.Stencil.Dim2
import Data.Array.Repa.Arbitrary
import Control.Monad.ST

import Test.Tasty
import Test.Tasty.QuickCheck

-- simple poisson matrix in two dimensions
matMul :: DIM2 -> Maybe Float
matMul !ix = case ix of
                Z :. -1 :. 0  -> Just $ -1
                Z :. 0  :. -1 -> Just $ -1
                Z :. 1  :. 0  -> Just $ -1
                Z :. 0  :. 1  -> Just $ -1
                Z :. 0  :. 0  -> Just $  4
                _             -> Nothing
{-# INLINE matMul #-}

problem_size :: Int
problem_size = 127*127

tests :: TestTree
tests = testGroup "Tests" [properties]

properties = testGroup "QuickCheck"
  [ testProperty "vcycle should improve the residual" $
      forAllUShaped (ix2 problem_size problem_size) (\vec ->
        let zeros = replicate (ix2 problem_size problem_size) 0
            res = vcycle (ix2 3 3) matMul zeros vec 32 1 2
            r_before = residual (ix2 3 3) matMul zeros vec
            r_end = residual (ix2 3 3) matMul res vec
        in  r_before >= r_end
      )
  ]

main = defaultMain tests
