//
//  sieve.h
//  primeSieve
//
//  Created by Thomas Daniell on 17/07/2016.
//  Copyright © 2016 Thomas Daniell. All rights reserved.
//

#ifndef SIEVE_H
#define SIEVE_H

bigInt countPrimes(const byte* sieve, bigInt range);
void runSieve(byte* sieve, bigInt range);

#endif
