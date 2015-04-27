{-# LANGUAGE OverloadedStrings #-}

module Utils where

import Plover.Types
import Plover.Macros
import Plover.Compile

constants :: CExpr
constants = seqList [
  Extern "GPS_L1_LAMBDA_NO_VAC" NumType
 ]

-- | Measure the integer ambiguity just from the code and carrier measurements.
-- The expectation value of @carrier + code / lambda@ is
-- @integer ambiguity + bias@.
--
-- Pseudorange bias can many wavelengths, so averaging @carrier + code@ isn't
-- sufficient for determining the ambiguity. Regardless of bias, this is an
-- important measurement. It is especially useful as a simple initialization
-- for filters.
simpleAmbMeas :: CExpr -- ^ Carrier phase measurement in wavelengths
              -> CExpr -- ^ Pseudorange (code) measurement in meters
              -> CExpr -- ^ Estimate of the integer ambiguity
simpleAmbMeas carrier code = carrier + code / "GPS_L1_LAMBDA_NO_VAC"

