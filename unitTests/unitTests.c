//
//  unitTests.c
//  primeSieve
//
//  Created by Thomas Daniell on 18/07/2016.
//  Copyright © 2016 Thomas Daniell. All rights reserved.
//

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#include "global.h"
#include "unitTests.h"
#include "smallSieve.h"
#include "sieve.h"

static void smallSieveTests(void);
static void primeSieveTests(void);

static void smallSieveTest(bigInt range, bigInt correctCount);
static void primeSieveTest(bigInt range, bigInt correctCount);

void runUnitTests(void) {
   printf("1) Conducting small sieve tests...\n");
   smallSieveTests();
   
   printf("2) Conducting main prime sieve tests...\n");
   primeSieveTests();
   
   printf("All tests passed. You are super awesome!\n\n");
}

static void smallSieveTests(void) {
   smallSieveTest(2, 1);
   smallSieveTest(3, 2);
   smallSieveTest(4, 2);
   smallSieveTest(5, 3);
   smallSieveTest(6, 3);
   smallSieveTest(7, 4);
   smallSieveTest(8, 4);
   smallSieveTest(9, 4);
   smallSieveTest(10, 4);
   smallSieveTest(11, 5);
   smallSieveTest(12, 5);
   smallSieveTest(20, 8);
   smallSieveTest(30, 10);
   smallSieveTest(50, 15);
   smallSieveTest(100, 25);
   smallSieveTest(150, 35);
   smallSieveTest(200, 46);
   smallSieveTest(205, 46);
   smallSieveTest(206, 46);
   smallSieveTest(207, 46);
   smallSieveTest(208, 46);
   smallSieveTest(209, 46);
   smallSieveTest(210, 46);
   smallSieveTest(211, 47);
   smallSieveTest(420, 81);
   smallSieveTest(500, 95);
   smallSieveTest(1000, 168);
   smallSieveTest(2310, 343);
   smallSieveTest(30030, 3248);
   smallSieveTest(100000, 9592);
   smallSieveTest(1000000, 78498);
   
   printf("- All small sieve tests passed.\n\n");
}

static void primeSieveTests(void) {
   for (bigInt i = MOD*MOD + 1; i <= MOD*MOD + 10; ++i) {
      primeSieveTest(i, 4590);
   }
   
   primeSieveTest(MOD*MOD + MOD + 1, 4612);
   primeSieveTest(MOD*MOD + 2*MOD + 1, 4627);
   
   primeSieveTest(211, 47);
   primeSieveTest(420, 81);
   primeSieveTest(500, 95);
   primeSieveTest(1000, 168);
   primeSieveTest(2310, 343);
   primeSieveTest(30030, 3248);
   primeSieveTest(100000, 9592);
   primeSieveTest(1000000, 78498);
   primeSieveTest(3437869, 246094);
   primeSieveTest(10000000, 664579);
   primeSieveTest(100000000, 5761455);
   primeSieveTest(347534534, 18678164);
   primeSieveTest(1000000000, 50847534);
   primeSieveTest(1124344445, 56831090);
   //primeSieveTest(10000000000, 455052511);
   
   printf("- All prime sieve tests passed.\n\n");
}

static void smallSieveTest(bigInt range, bigInt correctCount) {
   bigInt len = (range <= MOD) ? (MOD + 1) : (range + 1);
   smallInt* sieve = calloc(len, sizeof(smallInt));
   assert(sieve);
   
   runSmallSieve(sieve, range);
   assert(countSmallPrimes(sieve, range) == correctCount);
   free(sieve);
}

static void primeSieveTest(bigInt range, bigInt correctCount) {
   int len = getAllocSize(range);
   byte* sieve = calloc(len, sizeof(byte));
   assert(sieve);
   
   runSieve(sieve, range, len);
   int count = countPrimes(sieve, range, len);
   
   printf("Range %llu: count = %d & should be %llu\n", range, count, correctCount);
   assert(count == correctCount);
   free(sieve);
}



