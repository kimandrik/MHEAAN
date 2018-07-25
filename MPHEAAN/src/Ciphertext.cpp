/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#include "Ciphertext.h"

Ciphertext::Ciphertext(ZZ* axy, ZZ* bxy, long logp, long logq, long Nx, long Ny, long nx, long ny, bool isComplex) :
		axy(axy), bxy(bxy), logp(logp), logq(logq), Nx(Nx), Ny(Ny), nx(nx), ny(ny), isComplex(isComplex) {
}

Ciphertext::Ciphertext(const Ciphertext& o) : logp(o.logp), logq(o.logq), Nx(o.Nx), Ny(o.Ny), nx(o.nx), ny(o.ny), isComplex(o.isComplex) {
	long N = Nx * Ny;
	axy = new ZZ[N];
	bxy = new ZZ[N];
	for (long i = 0; i < N; ++i) {
		axy[i] = o.axy[i];
		bxy[i] = o.bxy[i];
	}
}

Ciphertext& Ciphertext::operator=(const Ciphertext& o) {
	if(this == &o) return *this;
	delete[] axy;
	delete[] bxy;
	Nx = o.Nx;
	Ny = o.Ny;
	nx = o.nx;
	ny = o.ny;
	isComplex = o.isComplex;
	logp = o.logp;
	logq = o.logq;
	long N = Nx * Ny;
	axy = new ZZ[N];
	bxy = new ZZ[N];
	for (long i = 0; i < N; ++i) {
		axy[i] = o.axy[i];
		bxy[i] = o.bxy[i];
	}
	return *this;
}

Ciphertext::~Ciphertext() {
	delete[] axy;
	delete[] bxy;
}
