{-# LANGUAGE BangPatterns, FlexibleContexts, ForeignFunctionInterface, ScopedTypeVariables #-}

import Prelude hiding (replicate)
import Math.Solver.Multigrid hiding (jacobi)
import Data.Array.Repa as R
import Data.Array.Repa.Stencil
import Data.Array.Repa.Stencil.Dim2
import Data.Array.Repa.Arbitrary
import Control.Monad.ST

offdiag' !ix = case ix of
                Z :. -1 :. 0  -> Just $ -1
                Z :. 0  :. -1 -> Just $ -1
                Z :. 1  :. 0  -> Just $ -1
                Z :. 0  :. 1  -> Just $ -1
                _             -> Nothing
{-# INLINE offdiag' #-}

offdiag = makeStencil (ix2 3 3) offdiag'
{-# INLINE offdiag #-}

jacobi !diag !offdiag !v !f = {-# SCC "jacobi" #-} (R.map (* (w / diag)) (f -^ mapStencil2 (BoundConst 0) offdiag v)) +^ (R.map (* (1-w)) v)
  where
    w = 2/3

{-# NOINLINE jacobi #-}

main =  do
  v <- computeUnboxedP (replicate (ix2 127 127) 0)
  f <- computeUnboxedP (replicate (ix2 127 127) 1)
  res :: Array U DIM2 Float <- computeUnboxedP $ jacobi 4 offdiag v f
  return ()
