/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#ifndef MHEAAN_PLAINTEXT_H_
#define MHEAAN_PLAINTEXT_H_

#include <NTL/ZZ.h>
#include "Params.h"

using namespace std;
using namespace NTL;

class Plaintext {
public:

	ZZ* mx = new ZZ[N];

	long logp;
	long n0;
	long n1;

	Plaintext(long logp = 0, long n0 = 0, long n1 = 0);

	virtual ~Plaintext();
};

#endif
