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

// Lookup table for 2*halfGaps^2 (not 100% whether this is less/more efficient)
static const byte halfGapSquare[NUM_WHEELS] = {
   50, 2, 8, 2, 8, 18, 2, 18, 8, 2, 8, 18, 18, 2, 18, 8, 2, 18, 8, 18, 32, 8, 2, 8,
   2, 8, 32, 18, 8, 18, 2, 8, 18, 2, 18, 18, 8, 2, 8, 18, 2, 18, 8, 2, 8, 2, 50, 2
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
   
   // 3) Eliminate composite numbers using SoE
   bigInt minPrimePos = MOD*(HALF_MOD + 1);
   bigInt bytePos = HALF_MOD;
   bigInt prime = 1 + MOD;
   
   while (bytePos <= sqrtHalfRange) {
      for (int indx = 0; indx < NUM_WHEELS; ++indx) {
         
         // If candidate (all lie in 1st bit segment) is prime, flag its multiples as !prime
         if (!(sieve[bytePos] & 1)) {
            bigInt multBytePos = minPrimePos;
            
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
         
         // Only primes on wheel spokes considered
         minPrimePos += prime*wheelGaps[indx] + halfGapSquare[indx];
         bytePos += halfGaps[indx];
         prime += wheelGaps[indx];
      }
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

static void wheelFactorise(byte* sieve, bigInt endByte, byte maxBit, int len) {
   // Eliminate all direct wheel multiples (skipping 1)
   int wheel = SECOND_WHEEL;
   int bytePos = 0;
   
   for (int wheelIndx = 1; wheelIndx < NUM_WHEELS; ++wheelIndx) {
      bytePos = squareTable[wheelIndx]; // Min value to be sieved
      
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
   }
}
