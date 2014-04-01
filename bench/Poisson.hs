{-# LANGUAGE BangPatterns, FlexibleContexts #-}

import Prelude hiding (replicate)
import Math.Solver.Multigrid
import Data.Array.Repa
import Data.Array.Repa.Stencil
import Data.Array.Repa.Stencil.Dim2
import Data.Array.Repa.Arbitrary
import Control.Monad.ST

import Criterion
import Criterion.Main

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

bench_vcycle :: Int -> Int -> Int -> Array U DIM1 Float
bench_vcycle size min_size num_smooth = runST $ computeUnboxedP $ vcycle matMul (replicate (ix1 (size * size)) 0) (replicate (ix1 (size * size)) 1) min_size 1 num_smooth

bench_vcycle_size size = bench_vcycle size 32 2

bench_vcycle_smoothing num_smooth = bench_vcycle 127 32 num_smooth

main = defaultMain
  [ bgroup "Problem Size"
    [ bench "127x127" $ whnf bench_vcycle_size 127
    , bench "255x255" $ whnf bench_vcycle_size 255
    , bench "511x511" $ whnf bench_vcycle_size 511
    ]
  , bgroup "Number of relaxation operations"
    [ bench "1" $ whnf bench_vcycle_smoothing 1
    , bench "2" $ whnf bench_vcycle_smoothing 2
    , bench "3" $ whnf bench_vcycle_smoothing 3
    , bench "4" $ whnf bench_vcycle_smoothing 4
    ]
  ]
