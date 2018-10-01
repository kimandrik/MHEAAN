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

	ZZ* mx;

	long logp;
	long logq;

	long N0;
	long N1;

	long n0;
	long n1;


	Plaintext(ZZ* mx = NULL, long logp = 0, long logq = 0, long N0 = 0, long N1 = 0, long n0 = 0, long n1 = 0);

	Plaintext(const Plaintext* o);

	~Plaintext();

};

#endif
