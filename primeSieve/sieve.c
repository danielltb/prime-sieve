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
#define SECOND_WHEEL 11
#define HALF_MOD MOD/2

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
static const smallInt squareTable[NUM_WHEELS] = {
   0, 60, 84, 144, 180, 264, 420, 480, 684, 840, 924, 1104, 1404,
   1740, 1860, 2244, 2520, 2664, 3120, 3444, 3960, 4704, 5100, 5304,
   5724, 5940, 6384, 7320, 8064, 8580, 9384, 9660, 10224, 11100,
   11400, 12324, 13284, 13944, 14280, 14964, 16020, 16380, 17484,
   18240, 18624, 19404, 19800, 21840
};

static const byte halfTable[NUM_WHEELS] = {
   105, 110, 111, 113, 114, 116, 119, 120, 123, 125, 126, 128, 131, 134,
   135, 138, 140, 141, 144, 146, 149, 153, 155, 156, 158, 159, 161, 165,
   168, 170, 173, 174, 176, 179, 180, 183, 186, 188, 189, 191, 194, 195,
   198, 200, 201, 203, 204, 209
};

static void wheelFactorise(byte* sieve, bigInt endByte, byte maxBit, int len);

/* Interface functions */
int getAllocSize(bigInt range) {
   return range/(2*BYTE_SIZE) + MOD;
}

void runSieve(byte* sieve, bigInt range, int len) {
   // 1) Range calculations
   bigInt sqrtHalfRange = (sqrt(range) - 1)/2;
   bigInt halfRange = (range - 1)/2;
   bigInt endByte = halfRange % len;
   byte maxBit = 1 << halfRange/len;
   
   // 2) Eliminate all wheel multiples from sieve
   wheelFactorise(sieve, endByte, maxBit, len);
   
   // Idea: Storing each offset in a different bit
   // 3) Apply SoE algorithm on sieve to flag all composite numbers as 0
   for (int indx = 0; indx < NUM_WHEELS; ++indx) {
      for (bigInt bytePos = halfTable[indx]; bytePos <= sqrtHalfRange; bytePos += HALF_MOD) {
         if (!(sieve[bytePos] & 1)) {
            bigInt prime = 2*bytePos + 1;
            bigInt multBytePos = 2*bytePos*(bytePos + 1);
            
            for (byte bit = 1; bit < maxBit; bit <<= 1) {
               for (; multBytePos < len; multBytePos += prime) {
                  sieve[multBytePos] |= bit;
               }
               
               multBytePos -= len; // Wrap around to start
            }
            
            for (; multBytePos <= endByte; multBytePos += prime) {
               sieve[multBytePos] |= maxBit;
            }
         }
      }
   }
}

// Optimise this function!!
bigInt countPrimes(const byte* sieve, bigInt range, int len) {
   bigInt primeCount = WHEEL_PRIMES;
   bigInt bytePos = HALF_MOD;
   byte bit = 1;
   
   // Range calculations
   bigInt halfRange = (range - 1)/2;
   byte maxBit = 1 << halfRange/len;
   int endByte = halfRange % len;
   int maxByte = (maxBit == 1) ? endByte : len - 1;
   
   // Iterate through wheel factorised sieve and count all primes
   while (TRUE) { // Swapping inner and outer loops?
      for (int indx = 0; indx < NUM_WHEELS; ++indx) {
         if (bytePos > maxByte) {
            bytePos -= len;
            
            if (bit < maxBit/2) {
               bit <<= 1;
            } else if (bit < maxBit) {
               maxByte = endByte;
               bit <<= 1;
            } else {
               return primeCount;
            }
         }
         
         primeCount += !(sieve[bytePos] & bit);
         bytePos += halfGaps[indx];
      }
   }
   
   return primeCount;
}

// Idea: Assume all numbers are not prime and only set ones lying in the wheel as candidates
static void wheelFactorise(byte* sieve, bigInt endByte, byte maxBit, int len) {
   // Eliminate all direct wheel multiples (skipping 1)
   int wheel = SECOND_WHEEL;
   int bytePos = 0;
   
   for (int indx = 1; indx < NUM_WHEELS; ++indx) {
      bytePos = squareTable[indx]; // Min value to be sieved
      
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
      
      wheel += wheelGaps[indx];
   }
}
