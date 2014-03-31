{-# LANGUAGE BangPatterns, FlexibleContexts #-}

module Main where

import Prelude hiding (replicate)
import Math.Solver.Multigrid
import Data.Array.Repa
import Data.Array.Repa.Stencil
import Data.Array.Repa.Stencil.Dim2
import Data.Array.Repa.Repr.Unboxed
import Control.Monad.ST

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
problem_size = 127

test_zero = runST $ equalsP (replicate (Z :. problem_size) 0) (vcycle matMul (replicate (Z :. problem_size) 0) (replicate (Z :. problem_size) 0) 32 1 2)

test_poisson = vcycle matMul (replicate (Z :. problem_size) 0) (replicate (Z :. problem_size) 1) 32 1 2

main = do
    res <- computeUnboxedP test_poisson
    print res
