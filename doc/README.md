# prime-sieve
Optimised sieve of eratosthenes project:

I am writing a program that counts all prime numbers up to a given range, using a sieve of eratosthenes algorithm,
for fun :D. Below is a list of the optimisations I have currently made to the sieve:
- Only iterating up to the square root of the input range in the outer loop
- When eliminating multiples of a prime in the inner loop, I start from prime^2
- Storing the sieve in the form of a list of boolean values (0 if prime, 1 if not)
- Using bitwise operations to store 1 boolean value per bit
- Storing each eighth of the sieve in a separate bit (less efficient to store bits contiguously)
- Using Modulo 210 wheel factorisation and only storing odd numbers in the sieve, enabling ~77%
  all numbers up to the bound to be skipped without missing any primes
- Wheel lookup tables to cut down on costly arithmetic operations

I am in the process of fixing a bug that causes the prime count to be off by one. When finished I
will optimise the algorithm further by using sieve segmentation (making the program more cache friendly).
