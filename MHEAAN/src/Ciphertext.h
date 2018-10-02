/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef MPHEAAN_CIPHERTEXT_H_
#define MPHEAAN_CIPHERTEXT_H_

#include <NTL/ZZ.h>
#include "Params.h"

using namespace std;
using namespace NTL;

class Ciphertext {
public:

	ZZ ax[N];
	ZZ bx[N];

	long logp;
	long logq;

	long n0;
	long n1;

	Ciphertext(long logp, long logq, long n0, long n1);

	Ciphertext(const Ciphertext* o);

	virtual ~Ciphertext();
};

#endif
