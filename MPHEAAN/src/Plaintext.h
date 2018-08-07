/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef MPHEAAN_PLAINTEXT_H_
#define MPHEAAN_PLAINTEXT_H_

#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;

class Plaintext {
public:

	ZZ* mxy;

	long logp;
	long logq;

	long Nx;
	long Ny;

	long nx;
	long ny;


	Plaintext(ZZ* mxy = NULL, long logp = 0, long logq = 0, long Nx = 0, long Ny = 0, long nx = 0, long ny = 0);

	Plaintext(const Plaintext& o);

	Plaintext& operator=(const Plaintext &o);

	~Plaintext();

};

#endif
