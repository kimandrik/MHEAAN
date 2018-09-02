/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "Ciphertext.h"

Ciphertext::Ciphertext(ZZ* axy, ZZ* bxy, long logp, long logq, long Nx, long Ny, long nx, long ny) :
		ax(axy), bx(bxy), logp(logp), logq(logq), Nx(Nx), Ny(Ny), nx(nx), ny(ny) {
}

Ciphertext::Ciphertext(const Ciphertext& o) : logp(o.logp), logq(o.logq), Nx(o.Nx), Ny(o.Ny), nx(o.nx), ny(o.ny) {
	long N = Nx * Ny;
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
	Nx = o.Nx;
	Ny = o.Ny;
	nx = o.nx;
	ny = o.ny;
	logp = o.logp;
	logq = o.logq;
	long N = Nx * Ny;
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
