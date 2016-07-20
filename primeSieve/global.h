//
//  global.h
//  primeSieve
//
//  Created by Thomas Daniell on 28/06/2016.
//  Copyright © 2016 Thomas Daniell. All rights reserved.
//

#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdint.h>

#define BYTE_SIZE 8
#define MOD 210

typedef uint64_t bigInt;
typedef uint16_t smallInt;
typedef uint8_t byte;

#endif

/* Sieve segmentation plan: 
 - Allocate a fixed segment size based on the given input range
 -
 
 */

// Consider endianness??

/* Working out for prime^2:
  2*bytePos*bytePos - 1;
 
*/

/* Original sieve specifications:
 For the time being, we always want to allocate 1 000 000 bytes of memory. The sieve is divided into 8 segments,
 each of which is stored in a unique bit across all 1 000 000 bytes. Every number in the sieve is represented by
 a unique bit and byte 'coordinate'. The numbers are divided as follows):
 - Bit #0    (int = 1): Numbers 0       to 999999  (Byte 0: 0       Byte 999999: 999999)  < 1000000
 - Bit #1    (int = 2): Numbers 1000000 to 1999999 (Byte 0: 1000000 Byte 999999: 1999999) < 2000000
 - Bit #2    (int = 4): Numbers 2000000 to 2999999 (Byte 0: 2000000 Byte 999999: 2999999) < 3000000
 - Bit #3    (int = 8): Numbers 3000000 to 3999999 (Byte 0: 3000000 Byte 999999: 3999999) < 4000000
 - Bit #4   (int = 16): Numbers 4000000 to 4999999 (Byte 0: 4000000 Byte 999999: 4999999) < 5000000
 - Bit #5   (int = 32): Numbers 5000000 to 5999999 (Byte 0: 5000000 Byte 999999: 5999999) < 6000000
 - Bit #6   (int = 64): Numbers 6000000 to 6999999 (Byte 0: 6000000 Byte 999999: 6999999) < 7000000
 - Bit #7  (int = 128): Numbers 7000000 to 7999999 (Byte 0: 7000000 Byte 999999: 7999999) < 8000000
 
 - Wheel modulus: 210
 - Minimum sieve range: 0
 - Maximum sieve range: 7999789 = 7999999 - 210 (makes overflow impossible)
 - Number of segments: 8
 - Maximum bit: 7 (128 as an int)
 - Maximum segment byte: 999999 (< 1000000)
 
 - Range value calculations: (let bytesAlloc = 1000000, modulus = 210)
 - Assert input range > 0 && < (bytesAlloc - modulus)
 - Base value of the last bit segment accessed: maxBit = 1 << range/bytesAlloc
 - Base value of the last bit segment fully sieved: maxFullBit = maxBit/2
 - Last byte accessed in sieve: endByte = range % bytesAlloc
 - Calculate sqrt(range) */


