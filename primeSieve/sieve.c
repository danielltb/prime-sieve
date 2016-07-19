//
//  sieve.c
//  primeSieve
//
//  Created by Thomas Daniell on 17/07/2016.
//  Copyright © 2016 Thomas Daniell. All rights reserved.
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "global.h"
#include "sieve.h"

#define WHEEL_PRIMES 46
#define NUM_WHEELS 48

#define FIRST_SQUARE_GAP 60
#define SECOND_WHEEL 11

#define HALF_MOD MOD/2
#define MOD_SQUARE HALF_MOD*HALF_MOD*2
#define ADD_GRAD MOD_SQUARE*2
#define ADD_INT MOD_SQUARE*3

// Lookup table for wheels used in SoE algorithm
static const byte wheelGaps[NUM_WHEELS] = {
   10, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4, 2, 6, 4, 6, 8, 4, 2, 4,
   2, 4, 8, 6, 4, 6, 2, 4, 6, 2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2, 10, 2
};

// Lookup table for (wheel - 1)/2
static const byte halfGaps[NUM_WHEELS] = {
   5, 1, 2, 1, 2, 3, 1, 3, 2, 1, 2, 3, 3, 1, 3, 2, 1, 3, 2, 3, 4, 2, 1, 2,
   1, 2, 4, 3, 2, 3, 1, 2, 3, 1, 3, 3, 2, 1, 2, 3, 1, 3, 2, 1, 2, 1, 5, 1
};

// Lookup table for (wheel^2 - 1)/2
static const smallInt squareGaps[NUM_WHEELS] = {
   60, 24, 60, 36, 84, 156, 60, 204, 156, 84, 180, 300, 336, 120, 384, 276, 144,
   456, 324, 516, 744, 396, 204, 420, 216, 444, 936, 744, 516, 804, 276, 564, 876,
   300, 924, 960, 660, 336, 684, 1056, 360, 1104, 756, 384, 780, 396, 2040, 0
};

static inline void clearMults(byte* sieve, bigInt prime, bigInt bytePos,
                              bigInt endByte, byte maxBit, int len);

static void prepSieve(byte* sieve, bigInt endByte, byte maxBit, int len);

/* Interface functions */
int getAllocSize(bigInt range) {
   return range/(2*BYTE_SIZE) + MOD;
}

void runSieve(byte* sieve, bigInt range, int len) {
   // 1) Range calculations
   bigInt halfRange = (range - 1)/2;
   byte maxBit = (byte)(1 << halfRange/len);
   bigInt sqrtHalfRange = (bigInt)((sqrt(range) - 1)/2);
   int endByte = halfRange % len;
   
   // 2) Eliminate all wheel multiples
   prepSieve(sieve, endByte, maxBit, len);
   
   // 3) Apply SoE algorithm on sieve to flag all composite numbers as 0
   bigInt incrSquare = MOD_SQUARE;
   bigInt prime = 1 + MOD;
   bigInt minPrimePos = 0;
   bigInt squareBase;
   bigInt count = 0;
   int bytePos = HALF_MOD;
   
   // Loop through all values of n (wheel fixed) in 'wheel + MOD*n'
   for (int incr = HALF_MOD; bytePos <= sqrtHalfRange; incr += HALF_MOD) {
      squareBase = 0;
      
      // Loop through all wheels (n fixed) in 'wheel + MOD*n'
      for (int indx = 0; indx < NUM_WHEELS; ++indx) {
         /* Since the largest prime candidate is sqrt(range) all their values
          lie in the first bit segment so we only need to check their primality
          by performing a bitwise and operation with value 1. */
         
         if (!(sieve[bytePos] & 1)) {
            // If number is prime, flag all its multiples as not prime (0)
            minPrimePos = incrSquare + squareBase + incr*2*(prime - incr*2); // This probs can be incremented
            clearMults(sieve, prime, minPrimePos, endByte, maxBit, len);
         }
         
         bytePos += halfGaps[indx];
         prime += wheelGaps[indx];
         squareBase += squareGaps[indx];
      }
      
      incrSquare += ADD_GRAD*count + ADD_INT;
      ++count;
   }
}

// Optimise this function!!
bigInt countPrimes(const byte* sieve, bigInt range, int len) {
   bigInt primeCount = WHEEL_PRIMES;
   int doLoop = 1;
   
   // Range calculations
   bigInt halfRange = (range - 1)/2;
   byte maxBit = (byte)(1 << halfRange/len);
   int endByte = halfRange % len;
   int maxByte = (maxBit == 1) ? endByte : len - 1;
   
   // Sieve position tracking variables
   bigInt wrapper = 0;
   int bytePos = 0;
   byte bit = 1;
   
   // Iterate through wheel factorised sieve and count all primes
   for (bigInt incr = HALF_MOD; doLoop; incr += HALF_MOD) {
      bytePos = incr - wrapper;
      
      for (int wheelIndx = 0; wheelIndx < NUM_WHEELS; ++wheelIndx) {
         if (bytePos > maxByte) {
            bytePos -= len;
            wrapper += len;
            
            if (bit < maxBit/2) {
               bit <<= 1;
            } else if (bit < maxBit) {
               maxByte = endByte;
               bit <<= 1;
            } else {
               doLoop = 0;
               break;
            }
         }
         
         primeCount += !(sieve[bytePos] & bit);
         bytePos += halfGaps[wheelIndx]; // Do I update this too?
      }
   }
   
   return primeCount;
}

/* Helper functions */
static inline void clearMults(byte* sieve, bigInt prime, bigInt bytePos,
                              bigInt endByte, byte maxBit, int len) {
   // Bit segments to be fully sieved
   for (byte bit = 1; bit < maxBit; bit <<= 1) {
      for (; bytePos < len; bytePos += prime) {
         sieve[bytePos] |= bit;
      }
      
      bytePos -= len; // Wrap around to start
   }
   
   // Final segment that is partially sieved
   for (; bytePos <= endByte; bytePos += prime) {
      sieve[bytePos] |= maxBit;
   }
}

static void prepSieve(byte* sieve, bigInt endByte, byte maxBit, int len) {
   // Eliminate all direct wheel multiples (skipping 1)
   int bytePos = 0, wheel = SECOND_WHEEL;
   bigInt squareBase = FIRST_SQUARE_GAP;
   
   for (int wheelIndx = 1; wheelIndx < NUM_WHEELS; ++wheelIndx) {
      bytePos = squareBase; // Min value to be sieved
      
      // Bit segments to be fully sieved
      for (byte bit = 1; bit < maxBit; bit <<= 1) {
         for (; bytePos < len; bytePos += wheel) {
            sieve[bytePos] |= bit;
         }
         
         bytePos -= len; // Wrap around to start
      }
      
      // Final segment that is partially sieved
      for (byte bit = maxBit; bytePos <= endByte; bytePos += wheel) {
         sieve[bytePos] |= bit;
      }
      
      wheel += wheelGaps[wheelIndx];
      squareBase += squareGaps[wheelIndx];
   }
}
