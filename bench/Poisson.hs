{-# LANGUAGE BangPatterns, FlexibleContexts, ForeignFunctionInterface, ScopedTypeVariables #-}

import Prelude hiding (replicate)
import Math.Solver.Multigrid
import Data.Array.Repa as R
import Data.Array.Repa.Stencil
import Data.Array.Repa.Stencil.Dim2
import Data.Array.Repa.Arbitrary
import Data.Array.Repa.Repr.ForeignPtr
import Control.Monad.ST

import Criterion
import Criterion.Main

import Foreign hiding (unsafePerformIO)
import Foreign.C.Types
import Foreign.Marshal.Alloc

-- simple poisson matrix in two dimensions
matMul !ix = case ix of
                Z :. -1 :. 0  -> Just $ -1
                Z :. 0  :. -1 -> Just $ -1
                Z :. 1  :. 0  -> Just $ -1
                Z :. 0  :. 1  -> Just $ -1
                Z :. 0  :. 0  -> Just $  4
                _             -> Nothing
{-# INLINE matMul #-}

offdiag' !ix = case ix of
                Z :. -1 :. 0  -> Just $ -1
                Z :. 0  :. -1 -> Just $ -1
                Z :. 1  :. 0  -> Just $ -1
                Z :. 0  :. 1  -> Just $ -1
                _             -> Nothing
{-# INLINE offdiag' #-}

offdiag :: Stencil DIM2 Float
offdiag = makeStencil (ix2 3 3) offdiag'
{-# INLINE offdiag #-}

foreign import ccall unsafe "jacobi.c jacobi"
    c_jacobi :: CInt -> Ptr CFloat -> Ptr CFloat -> Ptr CFloat -> IO ()

foreign import ccall unsafe "jacobi.c diagMul"
    c_offdiagmul :: CInt -> Ptr CFloat -> Ptr CFloat -> IO ()

bench_c_offdiagmul size vf rf = do
    withForeignPtr (toForeignPtr vf) (\vfp ->
      withForeignPtr (toForeignPtr rf) (\rfp ->
          c_offdiagmul (fromIntegral $ size * size) vfp rfp
      ))
    return rf

bench_offdiagmul size v = return $ computeUnboxedS $ mapStencil2 (BoundConst 0) offdiag v

bench_c_jacobi size vf ff rf = do
    withForeignPtr (toForeignPtr vf) (\vfp ->
      withForeignPtr (toForeignPtr ff) (\ffp ->
        withForeignPtr (toForeignPtr rf) (\rfp ->
          c_jacobi (fromIntegral $ size * size) vfp ffp rfp
        )))
    return rf

jacobi' !diag !offdiag !v !f = (smap (* (w / diag)) (f -^- d)) +^+ (R.map (* (1-w)) v)
  where
    d = mapStencil2 (BoundConst 0) offdiag v
    w = 2/3
    a -^- b = szipWith (-) a b
    a +^+ b = szipWith (+) b a

bench_jacobi size v f = return $ computeUnboxedS $ jacobi 4 offdiag v f

bench_vcycle size min_size num_smooth = do
    v <- computeUnboxedP $ replicate (ix2 size size) 0
    f <- computeUnboxedP $ replicate (ix2 size size) 1
    computeUnboxedP $ vcycle (ix2 3 3) matMul v f min_size 1 num_smooth

bench_vcycle_size size = bench_vcycle size 32 2

bench_vcycle_smoothing num_smooth = bench_vcycle 127 32 num_smooth

main = do
  let size = 1024
  v :: Array U DIM2 Float <- computeP (replicate (ix2 size size) 0)
  f :: Array U DIM2 Float <- computeP (replicate (ix2 size size) 1)
  vf :: Array F DIM2 CFloat  <- computeP (replicate (ix2 size size) 0)
  ff :: Array F DIM2 CFloat <- computeP (replicate (ix2 size size) 1)
  rf :: Array F DIM2 CFloat <- computeP (replicate (ix2 size size) 0)
  defaultMain
    [ bgroup "Jacobi size"
      [ bench "1024" $ whnfIO $ bench_jacobi size v f
      , bench "1024c" $ whnfIO $ bench_c_jacobi size vf ff rf
      ]
    , bgroup "Offdiag Mult"
      [ bench "1024" $ whnfIO $ bench_offdiagmul size v
      , bench "1024c" $ whnfIO $ bench_c_offdiagmul size vf rf
      ]
    , bgroup "Problem Size"
      [ bench "127x127" $ whnfIO $ bench_vcycle_size 127
      , bench "255x255" $ whnfIO $ bench_vcycle_size 255
      , bench "511x511" $ whnfIO $ bench_vcycle_size 511
      ]
    , bgroup "Number of relaxation operations"
      [ bench "1" $ whnfIO $ bench_vcycle_smoothing 1
      , bench "2" $ whnfIO $ bench_vcycle_smoothing 2
      , bench "3" $ whnfIO $ bench_vcycle_smoothing 3
      , bench "4" $ whnfIO $ bench_vcycle_smoothing 4
      ]
    ]
