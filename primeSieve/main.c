//
//  main.c
//  primeSieve
//
//  Created by Thomas Daniell on 17/07/2016.
//  Copyright © 2016 Thomas Daniell. All rights reserved.
//

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "global.h"
#include "sieve.h"
#include "smallSieve.h"
#include "unitTests.h"

#define LIMIT 2*(BYTES_ALLOC*8 - MOD) - 1

int main(void) {
   runUnitTests();
   
   clock_t start, end;
   bigInt primeCount;
   bigInt range;
   
   printf("Enter a bound to sieve up to: ");
   scanf("%u", &range);
      
   if ((range > 1) && (range <= LIMIT)) {
      printf("\nGenerating prime numbers up to the bound %u...\n", range);
      
      if (range > MOD) {
         byte* sieve = calloc(BYTES_ALLOC, sizeof(byte));
      
         if (!sieve) {
            fputs("Sieve memory allocation failed.\n", stderr);
            exit(EXIT_FAILURE);
         }
         
         start = clock();
         runSieve(sieve, range);
         end = clock();
         
         primeCount = countPrimes(sieve, range);
         free(sieve);
         
      } else {
         smallInt sieve[MOD + 1] = {0};
         
         start = clock();
         runSmallSieve(sieve, range);
         end = clock();
         
         primeCount = countSmallPrimes(sieve, range);
      }
      
      double diff = (double)(end - start)/CLOCKS_PER_SEC;
      printf("\nSieve generated in %lf seconds, counting primes...\n", diff);
      printf("%d primes counted between 0 and %u\n", primeCount, range);
   } else if (range == 1) {
      printf("Bound must be at least two to count any prime numbers.\n");
   } else {
      printf("Bound specified is not in the valid range.\n");
   }
   
   return EXIT_SUCCESS;
}
