//
//  sieve.h
//  primeSieve
//
//  Created by Thomas Daniell on 17/07/2016.
//  Copyright © 2016 Thomas Daniell. All rights reserved.
//

#ifndef SIEVE_H
#define SIEVE_H

bigInt countPrimes(const byte* sieve, bigInt range, bigInt len);
void runSieve(byte* sieve, bigInt range, bigInt len);
bigInt getAllocSize(bigInt range);

#endif
