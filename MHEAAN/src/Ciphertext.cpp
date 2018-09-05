/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "Ciphertext.h"

Ciphertext::Ciphertext(ZZ* ax, ZZ* bx, long logp, long logq, long N0, long N1, long n0, long n1) :
		ax(ax), bx(bx), logp(logp), logq(logq), N0(N0), N1(N1), n0(n0), n1(n1) {
}

Ciphertext::Ciphertext(const Ciphertext& o) : logp(o.logp), logq(o.logq), N0(o.N0), N1(o.N1), n0(o.n0), n1(o.n1) {
	long N = N0 * N1;
	ax = new ZZ[N];
	bx = new ZZ[N];
	for (long i = 0; i < N; ++i) {
		ax[i] = o.ax[i];
		bx[i] = o.bx[i];
	}
}

Ciphertext& Ciphertext::operator=(const Ciphertext& o) {
	if(this == &o) return *this;
	delete[] ax;
	delete[] bx;
	N0 = o.N0;
	N1 = o.N1;
	n0 = o.n0;
	n1 = o.n1;
	logp = o.logp;
	logq = o.logq;
	long N = N0 * N1;
	ax = new ZZ[N];
	bx = new ZZ[N];
	for (long i = 0; i < N; ++i) {
		ax[i] = o.ax[i];
		bx[i] = o.bx[i];
	}
	return *this;
}

Ciphertext::~Ciphertext() {
	delete[] ax;
	delete[] bx;
}
