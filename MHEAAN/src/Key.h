/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef MPHEAAN_KEY_H_
#define MPHEAAN_KEY_H_

#include "Params.h"
#include <NTL/ZZ.h>

using namespace NTL;

class Key {
public:

	uint64_t rax[Nnprimes];
	uint64_t rbx[Nnprimes];

	Key();

	virtual ~Key();
};

#endif
