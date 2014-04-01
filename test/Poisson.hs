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
matMul :: (Source a Float) => Array a DIM1 Float -> Float -> Array D DIM1 Float
matMul !x !h = to1 $ mapStencil2 (BoundConst 0) sten $ to2 x
    where
      sten = makeStencil2 3 3 func
      func = \ !ix -> case ix of
                          Z :. -1 :. 0  -> Just $ -1/(h*h)
                          Z :. 0  :. -1 -> Just $ -1/(h*h)
                          Z :. 1  :. 0  -> Just $ -1/(h*h)
                          Z :. 0  :. 1  -> Just $ -1/(h*h)
                          Z :. 0  :. 0  -> Just $ 4/(h*h)
                          _             -> Nothing
      {-# INLINE func #-}
{-# INLINE matMul #-}

problem_size :: Int
problem_size = 127*127

tests :: TestTree
tests = testGroup "Tests" [properties]

properties = testGroup "QuickCheck"
  [ testProperty "vcycle should improve the residual" $
      forAllUShaped (ix1 problem_size) (\vec ->
        let zeros = replicate (ix1 problem_size) 0
            res = vcycle matMul zeros vec 32 1 2
            r_before = residual (\a -> matMul a 1) zeros vec
            r_end = residual (\a -> matMul a 1) res vec
        in  r_before >= r_end
      )
  ]

main = defaultMain tests
