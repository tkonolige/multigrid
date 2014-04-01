{-# LANGUAGE QuasiQuotes,
             NoMonomorphismRestriction,
             FlexibleContexts,
             BangPatterns,
             ForeignFunctionInterface,
             ScopedTypeVariables,
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

-- TODO: try structed zip with

-- | residual of Ax=b
residual :: (Source a2 b, Source a3 b, Source a4 b, Floating b, Elt b, Unbox b)
         => (forall a. Source a b => Array a DIM1 b -> Array a2 DIM1 b)
         -> Array a3 DIM1 b
         -> Array a4 DIM1 b
         -> b
residual !matMul !x !b = runST $ foldAllP (+) 0 $ R.map (**2) (b -^ matMul x)

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
solve :: (Source a Float, Source b Float) => Array a DIM1 Float -> Array b DIM1 Float -> Float -> Array D DIM1 Float
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
       => (forall a. (Source a Float)
          => Array a DIM1 Float
          -> Float
          -> Array D DIM1 Float) -- * matrix multiplication using array and grid size
       -> Array b DIM1 Float     -- * initial guess at solution
       -> Array c DIM1 Float     -- * right hand side of Ax=b
       -> Int                    -- * size to perform direct solve on
       -> Float                  -- * initial size of each grid cell
       -> Int                    -- * number of time to relax at each level
       -> Array D DIM1 Float     -- * solution to linear system of equations
vcycle !matMul !v !f !size !h !nRelax = let (Z :. s) = extent v
                                        in case s <= size of
  False -> runST $ do
    v_ <- computeUnboxedP $ relax v f nRelax h
    f1 <- computeUnboxedP $ restrict $ f -^ matMul v_ h -- next levels residual
    vN <- computeUnboxedP $ restrict v_
    v1 <- computeUnboxedP $ vcycle matMul vN f1 size (h*2) nRelax -- recurse on next level
    v' <- computeUnboxedP $ v_ +^ interpolate v1 -- take next level and add it to current level
    return $ relax v' f nRelax h -- relax again
  True -> solve v f h -- direct solve

relax :: (Fractional c, Unbox c, Source a c, Source b c)
      => Array a DIM1 c
      -> Array b DIM1 c
      -> Int
      -> c
      -> Array D DIM1 c
relax !v  _  0  _ = delay v
relax !v !f !n !h = runST $ do
  relaxed <- computeUnboxedP $ jacobi v f h
  return $ relax relaxed f (n-1) h
{-# INLINE relax #-}

-- weighted jacobi relaxation
-- TODO: consider structed map
-- TODO: make sure first array is fully evaluated
jacobi :: (Fractional c, Unbox c, Source a c, Source b c)
       => Array a DIM1 c
       -> Array b DIM1 c
       -> c
       -> Array D DIM1 c
jacobi !v !f !h = (R.map (* (w / (4/(h*h)))) (f -^ (to1 $ mapStencil2 (BoundConst 0) sten (to2 v)))) +^ (R.map (* (1-w)) v)
  where
    w = 2/3
    (Z :. s)   = extent v
    width      = floor $ sqrt $ fromIntegral s
    val        = -1/(h*h)
    func !ix = case ix of
                  Z :. -1 :. 0  -> Just val
                  Z :. 0  :. -1 -> Just val
                  Z :. 1  :. 0  -> Just val
                  Z :. 0  :. 1  -> Just val
                  _             -> Nothing
    {-# INLINE func #-}
    sten = makeStencil2 3 3 func
{-# INLINE jacobi #-}

-- sor smoothing
-- V'(i,j)
-- = (1-w)V(i,j)+w/4*(V(i+1,j)+V'(i-1,j)+V(i,j+1)+V'(i,j-1)+h^2*p(i,j))
-- sor :: (Fractional c, Unbox c, SOurce a c, Source b c) => Array a DIM1 c -> Array b DIM1 c -> Int -> c -> Array D DIM1 c
-- sor !v !f !n !h =
--   where
--     relaxedOdd =

restrict :: (Fractional b, Source a b) => Array a DIM1 b -> Array D DIM1 b
restrict !x = to1 $ shrinkVec $ to2 x
{-# INLINE restrict #-}

interpolate :: (Fractional b, Source a b) => Array a DIM1 b -> Array D DIM1 b
interpolate !x = to1 $ growVec $ to2 x
{-# INLINE interpolate #-}

shrinkVec :: (Fractional a, Source b a) => Array b DIM2 a -> Array D DIM2 a
shrinkVec !v = unsafeTraverse v newDim (\ !f !(Z :. x :. y) ->
                    (f (Z :. 2 * x :. 2 * y) + f (Z :. 2 * x :. 2 * y + 2) +
                     f (Z :. 2 * x + 2 :. 2 * y) + f (Z :. 2 * x + 2 :. 2 * y + 2) +
                     2 * (f (Z :. 2 * x + 1 :. 2 * y) + f (Z :. 2 * x + 1:. 2 * y + 2) +
                          f (Z :. 2 * x :. 2 * y + 1) + f (Z :. 2 * x + 2 :. 2 * y + 1))
                     + 4 * f (Z :. 2 * x + 1 :. 2 * y + 1)) / 16)
  where
    newDim !(Z :. x :. y) = (Z :. x `div` 2 :. y `div` 2)
{-# INLINE shrinkVec #-}

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
{-# INLINE growVec #-}
