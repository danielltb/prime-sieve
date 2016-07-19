//
//  smallSieve.c
//  primeSieve
//
//  Created by Thomas Daniell on 18/07/2016.
//  Copyright © 2016 Thomas Daniell. All rights reserved.
//

#include <stdio.h>
#include <math.h>

#include "global.h"
#include "smallSieve.h"

void runSmallSieve(smallInt* sieve, int range) {
   int sqrtRange = sqrt(range);
   int rangeCast = range;
   
   for (int prime = 2; prime <= sqrtRange; ++prime) {
      if (!sieve[prime]) {
         for (int comp = prime*prime; comp <= rangeCast; comp += prime) {
            sieve[comp] = 1;
         }
      }
   }
}

int countSmallPrimes(smallInt* sieve, int range) {
   int rangeCast = (int)range;
   int primeCount = 0;
   
   for (int prime = 2; prime <= rangeCast; ++prime) {
      primeCount += !sieve[prime];
   }
   
   return primeCount;
}
