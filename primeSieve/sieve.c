//
//  sieve.c
//  primeSieve
//
//  Created by Thomas Daniell on 17/07/2016.
//  Copyright © 2016 Thomas Daniell. All rights reserved.
//  Do testing, allocator and small sieve next!!

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "global.h"
#include "sieve.h"

#define WHEEL_PRIMES 46
#define NUM_WHEELS 48

#define HALF_MOD MOD/2
#define MOD_SQUARE HALF_MOD*HALF_MOD*2
#define ADD_GRAD MOD_SQUARE*2
#define ADD_INT MOD_SQUARE*3

// Lookup table for wheels used in SoE algorithm
static const byte wheelTable[NUM_WHEELS] = {
   1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
   71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 121, 127, 131,
   137, 139, 143, 149, 151, 157, 163, 167, 169, 173, 179, 181, 187,
   191, 193, 197, 199, 209
};

// Lookup table for (wheel - 1)/2
static const byte halfTable[NUM_WHEELS] = {
   0, 5, 6, 8, 9, 11, 14, 15, 18, 20, 21, 23, 26, 29, 30, 33,
   35, 36, 39, 41, 44, 48, 50, 51, 53, 54, 56, 60, 63, 65, 68,
   69, 71, 74, 75, 78, 81, 83, 84, 86, 89, 90, 93, 95, 96, 98,
   99, 104
};

// Lookup table for (wheel^2 - 1)/2
static const smallInt squareTable[NUM_WHEELS] = {
   0, 60, 84, 144, 180, 264, 420, 480, 684, 840, 924, 1104, 1404,
   1740, 1860, 2244, 2520, 2664, 3120, 3444, 3960, 4704, 5100, 5304,
   5724, 5940, 6384, 7320, 8064, 8580, 9384, 9660, 10224, 11100,
   11400, 12324, 13284, 13944, 14280, 14964, 16020, 16380, 17484,
   18240, 18624, 19404, 19800, 21840
};

static inline void clearMults(byte* sieve, bigInt prime, bigInt bytePos,
                              bigInt endByte, byte maxBit);

static void prepSieve(byte* sieve, bigInt endByte, byte maxBit);

/* Interface functions */

void runSieve(byte* sieve, bigInt range) {
   // 1) Range calculations
   bigInt halfRange = (range - 1)/2;
   byte maxBit = (byte)(1 << halfRange/BYTES_ALLOC);
   bigInt endByte = halfRange % BYTES_ALLOC;
   bigInt sqrtHalfRange = (bigInt)((sqrt(range) - 1)/2);
   
   // 2) Eliminate all wheel multiples
   prepSieve(sieve, endByte, maxBit);
   
   // 3) Apply SoE algorithm on sieve to flag all composite numbers as 0
   bigInt incrSquare = MOD_SQUARE;
   bigInt bytePos = 0, count = 0;
   bigInt minPrimePos = 0;
   
   // Loop through all values of n (wheel fixed) in 'wheel + MOD*n'
   for (bigInt incr = HALF_MOD; bytePos <= sqrtHalfRange; incr += HALF_MOD) {
      
      // Loop through all wheels (n fixed) in 'wheel + MOD*n'
      for (int indx = 0; indx < NUM_WHEELS; ++indx) {
         bytePos = halfTable[indx] + incr;
         
         /* Since the largest prime candidate is sqrt(range) all their values
          lie in the first bit segment so we only need to check their primality
          by performing a bitwise and operation with value 1. */
         
         if (!(sieve[bytePos] & 1)) {
            // If number is prime, flag all its multiples as not prime (0)
            minPrimePos = incrSquare + squareTable[indx] + incr*wheelTable[indx]*2;
            clearMults(sieve, wheelTable[indx] + incr*2, minPrimePos, endByte, maxBit);
         }
      }
      
      incrSquare += ADD_GRAD*count + ADD_INT;
      ++count;
   }
}

bigInt countPrimes(const byte* sieve, bigInt range) {
   // Range calculations
   bigInt halfRange = (range - 1)/2;
   byte maxBit = (byte)(1 << halfRange/BYTES_ALLOC);
   bigInt endByte = halfRange % BYTES_ALLOC;
   bigInt maxByte = (maxBit == 1) ? endByte: BYTES_ALLOC - 1;
   
   bigInt primeCount = WHEEL_PRIMES;
   bool continueSieve = true;
   bigInt bytePos = 0;
   bigInt wrapper = 0;
   byte bit = 1;
   
   for (bigInt incr = HALF_MOD; continueSieve; incr += HALF_MOD) {
      for (int wheelIndx = 0; wheelIndx < NUM_WHEELS; ++wheelIndx) {
         bytePos = halfTable[wheelIndx] + incr - wrapper;
         
         if (bytePos > maxByte) {
            bytePos -= BYTES_ALLOC;
            wrapper += BYTES_ALLOC;
            
            if (bit < maxBit/2) {
               bit <<= 1;
            } else if (bit < maxBit) {
               maxByte = endByte;
               bit <<= 1;
            } else {
               continueSieve = false;
               break;
            }
         }
         
         primeCount += !(sieve[bytePos] & bit);
      }
   }
   
   return primeCount;
}

/* Helper functions */
static inline void clearMults(byte* sieve, bigInt prime, bigInt bytePos,
                              bigInt endByte, byte maxBit) {
   // Bit segments to be fully sieved
   for (byte bit = 1; bit < maxBit; bit <<= 1) {
      for (; bytePos < BYTES_ALLOC; bytePos += prime) {
         sieve[bytePos] |= bit;
      }
      
      bytePos -= BYTES_ALLOC; // Wrap around to start
   }
   
   // Final segment that is partially sieved
   for (; bytePos <= endByte; bytePos += prime) {
      sieve[bytePos] |= maxBit;
   }
}

static void prepSieve(byte* sieve, bigInt endByte, byte maxBit) {
   // Eliminate all direct wheel multiples (skipping 1)
   bigInt bytePos = 0, wheel = 0;
   
   for (int wheelIndx = 1; wheelIndx < NUM_WHEELS; ++wheelIndx) {
      wheel = wheelTable[wheelIndx];
      bytePos = squareTable[wheelIndx]; // Min value to be sieved
      
      // Bit segments to be fully sieved
      for (byte bit = 1; bit < maxBit; bit <<= 1) {
         for (; bytePos < BYTES_ALLOC; bytePos += wheel) {
            sieve[bytePos] |= bit;
         }
         
         bytePos -= BYTES_ALLOC; // Wrap around to start
      }
      
      // Final segment that is partially sieved
      for (byte bit = maxBit; bytePos <= endByte; bytePos += wheel) {
         sieve[bytePos] |= bit;
      }
   }
}
