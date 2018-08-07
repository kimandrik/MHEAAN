/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "Plaintext.h"

Plaintext::Plaintext(ZZ* mxy, long logp, long logq, long Nx, long Ny, long nx, long ny)
		: mxy(mxy), logp(logp), logq(logq), Nx(Nx), Ny(Ny), nx(nx), ny(ny) {
}

Plaintext::Plaintext(const Plaintext& o) : mxy(o.mxy), logp(o.logp), logq(o.logq), Nx(o.Nx), Ny(o.Ny), nx(o.nx), ny(o.ny) {
	long N = Nx * Ny;
	mxy = new ZZ[N];
	for (long i = 0; i < N; ++i) {
		mxy[i] = o.mxy[i];
	}
}

Plaintext& Plaintext::operator=(const Plaintext& o) {
	if(this == &o) return *this;
	delete[] mxy;
	Nx = o.Nx;
	Ny = o.Ny;
	nx = o.nx;
	ny = o.ny;
	logp = o.logp;
	logq = o.logq;
	long N = Nx * Ny;
	mxy = new ZZ[N];
	for (long i = 0; i < N; ++i) {
		mxy[i] = o.mxy[i];
	}
	return *this;
}

Plaintext::~Plaintext() {
	delete[] mxy;
}



