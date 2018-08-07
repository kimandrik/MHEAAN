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

using namespace std;
using namespace NTL;

/**
 * Ciphertext if MRLWE(RLWE) instanc
 */
class Ciphertext {
public:

	ZZ* axy;
	ZZ* bxy;

	long logp;
	long logq;

	long Nx;
	long Ny;

	long nx;
	long ny;

	Ciphertext(ZZ* axy = NULL, ZZ* bxy = NULL, long logp = 0, long logq = 0, long Nx = 0, long Ny = 0, long nx = 0, long ny = 0);

	Ciphertext(const Ciphertext& o);

	Ciphertext& operator=(const Ciphertext &o);

	~Ciphertext();
};

#endif
