/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef MPHEAAN_PLAINTEXT_H_
#define MPHEAAN_PLAINTEXT_H_

#include "Params.h"
#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;

class Plaintext {
public:

	ZZ mx[N];

	long logp;
	long logq;

	long n0;
	long n1;

	Plaintext(long logp, long logq, long n0, long n1);

	~Plaintext();

};

#endif
