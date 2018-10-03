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

	Ciphertext(long logp = 0, long logq = 0, long n0 = 0, long n1 = 0);

	Ciphertext(const Ciphertext* o);

	void copyParams(Ciphertext* o);

	void copy(Ciphertext* o);

	void free();

	virtual ~Ciphertext();
};

#endif
