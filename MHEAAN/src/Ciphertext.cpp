/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "Ciphertext.h"

Ciphertext::Ciphertext(long logp, long logq, long n0, long n1) : logp(logp), logq(logq), n0(n0), n1(n1) {
}

Ciphertext::Ciphertext(const Ciphertext* o) : logp(o->logp), logq(o->logq), n0(o->n0), n1(o->n1) {
	for (long i = 0; i < N; ++i) {
		ax[i] = o->ax[i];
		bx[i] = o->bx[i];
	}
}

Ciphertext::~Ciphertext() {
}
