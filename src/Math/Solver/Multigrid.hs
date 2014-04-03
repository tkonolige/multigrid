{-# LANGUAGE FlexibleContexts,
             BangPatterns,
             RankNTypes
  #-}

-----------------------------------------------------------------------------
-- |
-- Copyright   :  (C) 2014 Tristan Konolige
-- License     :  MIT (see the file LICENSE)
-- Maintainer  :  Tristan Konolige (tristan.konolige@gmail.com)
-- Stability   :  experimental
--
-- A multigrid solver for systems of partial differential equations.
--
----------------------------------------------------------------------------

module Math.Solver.Multigrid
  (
    vcycle
  , to2
  , to1
  , replicate
  , residual
  ) where

import Prelude hiding (replicate)
import Data.Array.Repa as R
import Data.Array.Repa.Eval
import Data.Array.Repa.Stencil
import Data.Array.Repa.Stencil.Dim2
import Data.Array.Repa.Repr.ForeignPtr
import Data.Array.Repa.Repr.Unboxed
import Data.Array.Repa.Operators.Traversal
import Data.Array.Repa.Operators.IndexSpace
import Control.Monad.ST
import Control.Applicative

-- | residual of Ax=b
residual :: (Source a1 b, Source a2 b, Floating b, Elt b, Unbox b)
         => DIM2
         -> (DIM2 -> Maybe b)
         -> Array a1 DIM2 b
         -> Array a2 DIM2 b
         -> b
residual !sh !matMul !x !b = runST $ foldAllP (+) 0 $ R.map (**2) (b -^ mapStencil2 (BoundConst 0) (makeStencil sh matMul) x)

-- | convert a one dimensional array to a square grid
to2 :: Source a b => Array a DIM1 b -> Array D DIM2 b
to2 !v = reshape (Z :. width :. width) v
  where
    (Z :. l) = extent v
    width = floor $ sqrt $ fromIntegral l
{-# INLINE to2 #-}

-- | convert a square grid to a one dimensional array
to1 :: Source a b => Array a DIM2 b -> Array D DIM1 b
to1 !v = reshape (Z :. l) v
  where
    (Z :. width :. height) = extent v
    l = height * width
{-# INLINE to1 #-}

-- | an array full of a single value
replicate :: (Shape sh) => sh -> a -> Array D sh a
replicate sh a = fromFunction sh (\_ -> a)

-- solve Ax=b directly TODO: use different ways to solve
solve :: (Source a Float, Source b Float, Shape sh) => Array a sh Float -> Array b sh Float -> Float -> Array D sh Float
solve !v !f !h = delay v
-- solve !v !f !h = sol'
--   where
--     (Z :. s) = extent v
--     mat = buildMat s h
--     b = buildMatrix s 1 (\ !(r,c) -> f `unsafeIndex` (ix1 r))
--     sol = single $ linearSolveR (double mat) (double b)
--     sol' = fromFunction (Z :. s) (\ !(Z :. i) -> sol @@> (i, 0))
-- {-# INLINE solve #-}
--
-- buildMat !ss !h = buildMatrix ss ss (\ !coord -> case coord of
--                           !(x,y) | x == y -> let cc = to2D x
--                                             in (valNeuman (top cc) + valNeuman (bot cc) +
--                                                 valNeuman (left cc) + valNeuman (right cc)) / (h*h)
--                           !(x,y) | (x `div` s == y `div` s && (abs((x `mod` s) - (y `mod` s)) == s ||
--                                    abs((x `mod` s) - (y `mod` s)) == 1)) ||
--                                    (x `mod` s == y `mod` s && (abs((x `div` s) - (y `div` s)) == s ||
--                                    abs((x `div` s) - (y `div` s)) == 1)) -> let cc = to2D x
--                                                                           in -1 * valNeuman cc / (h*h)
--                           otherwise -> 0
--                           )
--   where
--     s = floor $ sqrt $ fromIntegral ss
--     isNeuman !(x, y) = x < 0 || x >= s || y < 0 || y >= s
--     valNeuman !c = if isNeuman c then 1 else 1
--     to2D !x = (x `mod` s, x `div` s)
--     top !(x, y) = (x, y+1)
--     bot !(x, y) = (x, y-1)
--     left !(x, y) = (x+1, y)
--     right !(x, y) = (x-1, y)

-- | Perform one multigrid iteration using the vcycle scheme
vcycle :: (Source b Float, Source c Float)
       => DIM2                  -- * size of stencil
       -> (DIM2 -> Maybe Float) -- * stencil function
       -> Array b DIM2 Float    -- * initial guess at solution
       -> Array c DIM2 Float    -- * right hand side of Ax=b
       -> Int                   -- * size to perform direct solve on
       -> Float                 -- * initial size of each grid cell
       -> Int                   -- * number of time to relax at each level
       -> Array D DIM2 Float    -- * solution to linear system of equations
vcycle !stenSize !sten !v !f !size !h !nRelax = 
    let (Z :. s :. _) = extent v
    in case s <= size of
      False -> runST $ do
        v_ <- computeUnboxedP $ relax diag offdiag v f nRelax
        f_ <- computeUnboxedP $ f -^ mapStencil2 (BoundConst 0) matMul v_ -- residual
        f1 <- computeUnboxedP $ restrict f_ -- next levels residual
        v1 <- computeUnboxedP $ restrict v_ -- next levels solution
        vN <- computeUnboxedP $ vcycle stenSize sten v1 f1 size (h*2) nRelax -- recurse on next level
        v' <- computeUnboxedP $ v_ +^ interpolate vN -- take next level and add it to current level
        return $ relax diag offdiag v' f_ nRelax -- relax again
      True -> solve v f h -- direct solve
  where
    matMul = makeStencil stenSize (\ !i -> (/(h*h)) <$> sten i)
    Just diag = (/(h*h)) <$> sten (Z :. 0 :. 0)
    offdiag = makeStencil stenSize (\ !i@(Z :. x :. y) -> case x == y of
                                                         True -> Nothing
                                                         False -> (/(h*h)) <$> sten i
                                   )

relax :: (Fractional c, Unbox c, Source a c, Source b c)
      => c
      -> Stencil DIM2 c
      -> Array a DIM2 c
      -> Array b DIM2 c
      -> Int
      -> Array D DIM2 c
relax  _     _       !v  _  0 = delay v
relax !diag !offdiag !v !f !n = runST $ do
  relaxed <- computeUnboxedP $ jacobi diag offdiag v f
  return $ relax diag offdiag relaxed f (n-1)

-- weighted jacobi relaxation
-- TODO: consider structed map
-- TODO: make sure first array is fully evaluated
jacobi :: (Fractional c, Unbox c, Source a c, Source b c)
       => c
       -> Stencil DIM2 c
       -> Array a DIM2 c
       -> Array b DIM2 c
       -> Array D DIM2 c
jacobi !diag !offdiag !v !f = (R.map (* (w / diag)) (f -^ mapStencil2 (BoundConst 0) offdiag v)) +^ (R.map (* (1-w)) v)
  where
    w = 2/3
    -- (Z :. s)   = extent v
    -- width      = floor $ sqrt $ fromIntegral s
    -- val        = -1/(h*h)
    -- func !ix = case ix of
    --               Z :. -1 :. 0  -> Just val
    --               Z :. 0  :. -1 -> Just val
    --               Z :. 1  :. 0  -> Just val
    --               Z :. 0  :. 1  -> Just val
    --               _             -> Nothing
    -- {-# INLINE func #-}
    -- sten = makeStencil2 3 3 func

-- sor smoothing
-- V'(i,j)
-- = (1-w)V(i,j)+w/4*(V(i+1,j)+V'(i-1,j)+V(i,j+1)+V'(i,j-1)+h^2*p(i,j))
-- sor :: (Fractional c, Unbox c, SOurce a c, Source b c) => Array a DIM1 c -> Array b DIM1 c -> Int -> c -> Array D DIM1 c
-- sor !v !f !n !h =
--   where
--     relaxedOdd =

restrict :: (Fractional b, Source a b) => Array a DIM2 b -> Array D DIM2 b
restrict !x = shrinkVec x

interpolate :: (Fractional b, Source a b) => Array a DIM2 b -> Array D DIM2 b
interpolate !x = growVec x

shrinkVec :: (Fractional a, Source b a) => Array b DIM2 a -> Array D DIM2 a
shrinkVec !v = unsafeTraverse v newDim (\ !f !(Z :. x :. y) ->
                    (f (Z :. 2 * x :. 2 * y) + f (Z :. 2 * x :. 2 * y + 2) +
                     f (Z :. 2 * x + 2 :. 2 * y) + f (Z :. 2 * x + 2 :. 2 * y + 2) +
                     2 * (f (Z :. 2 * x + 1 :. 2 * y) + f (Z :. 2 * x + 1:. 2 * y + 2) +
                          f (Z :. 2 * x :. 2 * y + 1) + f (Z :. 2 * x + 2 :. 2 * y + 1))
                     + 4 * f (Z :. 2 * x + 1 :. 2 * y + 1)) / 16)
  where
    newDim !(Z :. x :. y) = (Z :. x `div` 2 :. y `div` 2)

growVec :: (Fractional b, Source r b) => Array r DIM2 b -> Array D DIM2 b
growVec !v = unsafeTraverse v newDim (\ !f !(Z :. x :. y) ->
              let f' (Z :. x :. y) | x < 0 || x >= width || y < 0 || y >= height = 0
                  f' c = f c
                  {-# INLINE f' #-}
              in 1/4 * case (x `mod` 2, y `mod` 2) of
                (0, 0) -> f' (Z :. x `div` 2 - 1 :. y `div` 2) +
                          f' (Z :. x `div` 2 :. y `div` 2 - 1) +
                          f' (Z :. x `div` 2 - 1:. y `div` 2 - 1) +
                          f' (Z :. x `div` 2 :. y `div` 2)
                (1, 0) -> 2 * (f' (Z :. x `div` 2 :. y `div` 2) +
                               f' (Z :. x `div` 2 :. y `div` 2 - 1))
                (0, 1) -> 2 * (f' (Z :. x `div` 2 :. y `div` 2) +
                               f' (Z :. x `div` 2 - 1 :. y `div` 2))
                (1, 1) -> 4 * f' (Z :. x `div` 2 :. y `div` 2)
                )
  where
    (Z :. width :. height) = extent v
    newDim !(Z :. x :. y) = (Z :. x * 2 + 1 :. y * 2 + 1)
