/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#ifndef MPHEAAN_NUMB_H_
#define MPHEAAN_NUMB_H_

#pragma once

#include <iostream>
#include <vector>
#include <math.h>
#include <stdint.h>

#include "Common.h"

using namespace std;

void mulMod(uint64_t& r, uint64_t a, uint64_t b, uint64_t p);

void mulModBarrett(uint64_t& r, uint64_t a, uint64_t b, uint64_t p, uint64_t pr, long twok);

uint64_t powMod(uint64_t x, uint64_t y, uint64_t p);

uint64_t inv(uint64_t x);

uint32_t bitReverse(uint32_t x);

void findPrimeFactors(vector<uint64_t> &s, uint64_t number);

uint64_t findPrimitiveRoot(uint64_t m);

uint64_t findMthRootOfUnity(uint64_t M, uint64_t p);

#endif
