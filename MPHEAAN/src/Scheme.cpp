/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "Scheme.h"

#include "NTL/BasicThreadPool.h"

#include "StringUtils.h"
//-----------------------------------------

Scheme::Scheme(SecretKey& secretKey, Ring& ring) : ring(ring) {
	addEncKey(secretKey);
	addMultKey(secretKey);
};


//----------------------------------------------------------------------------------
//   KEYS GENERATION
//----------------------------------------------------------------------------------


void Scheme::addEncKey(SecretKey& secretKey) {
	ZZ* ex = new ZZ[ring.N];
	ZZ* ax = new ZZ[ring.N];
	ZZ* bx = new ZZ[ring.N];

	ring.sampleUniform(ax, ring.logQQ);
	ring.sampleGauss(ex);

	long np = ceil((1 + ring.logQQ + ring.logN + 3)/59.0);
	ring.mult(bx, secretKey.sx, ax, np, ring.QQ);
	ring.sub(bx, ex, bx, ring.QQ);

	np = ceil((2 * ring.logN + 4 * ring.logQ + 2)/59.0);
	uint64_t* rax = ring.toNTT(ax, np);
	uint64_t* rbx = ring.toNTT(bx, np);

	keyMap.insert(pair<long, Key>(ENCRYPTION, Key(rax, rbx)));

	delete[] ex;
	delete[] ax;
	delete[] bx;
}

void Scheme::addMultKey(SecretKey& secretKey) {
	ZZ* ex = new ZZ[ring.N];
	ZZ* ax = new ZZ[ring.N];
	ZZ* bx = new ZZ[ring.N];
	ZZ* sxsx = new ZZ[ring.N];

	long np = ceil((1 + 1 + ring.logN + 3)/59.0);
	ring.mult(sxsx, secretKey.sx, secretKey.sx, np, ring.Q);
	ring.leftShiftAndEqual(sxsx, ring.logQ, ring.QQ);
	ring.sampleUniform(ax, ring.logQQ);
	ring.sampleGauss(ex);
	ring.addAndEqual(ex, sxsx, ring.QQ);

	np = ceil((1 + ring.logQQ + ring.logN + 3)/59.0);
	ring.mult(bx, secretKey.sx, ax, np, ring.QQ);
	ring.sub(bx, ex, bx, ring.QQ);

	np = ceil((2 * ring.logN + 4 * ring.logQ + 2)/59.0);
	uint64_t* rax = ring.toNTT(ax, np);
	uint64_t* rbx = ring.toNTT(bx, np);

	keyMap.insert(pair<long, Key>(MULTIPLICATION, Key(rax, rbx)));

	delete[] ex;
	delete[] ax;
	delete[] bx;
	delete[] sxsx;

}

void Scheme::addLeftRotKey(SecretKey& secretKey, long rx, long ry) {
	ZZ* exy = new ZZ[ring.N];
	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ZZ* rysxy = ring.leftRotate(secretKey.sx, rx, ry);

	ring.leftShiftAndEqual(rysxy, ring.logQ, ring.QQ);
	ring.sampleUniform(axy, ring.logQQ);
	ring.sampleGauss(exy);
	ring.addAndEqual(exy, rysxy, ring.QQ);

	long np = ceil((1 + ring.logQQ + ring.logN + 3)/59.0);
	ring.mult(bxy, secretKey.sx, axy, np, ring.QQ);
	ring.sub(bxy, exy, bxy, ring.QQ);

	np = ceil((2 * ring.logQQ + ring.logN + 3)/59.0);
	uint64_t* rax = ring.toNTT(axy, np);
	uint64_t* rbx = ring.toNTT(bxy, np);

	leftRotKeyMap.insert(pair<pair<long, long>, Key>({rx, ry}, Key(rax, rbx)));

	delete[] exy;
	delete[] axy;
	delete[] bxy;
	delete[] rysxy;
}

void Scheme::addLeftXRotKeys(SecretKey& secretKey) {
	for (long ix = 1; ix < ring.N0h; ix <<= 1) {
		if(leftRotKeyMap.find({ix, 0}) == leftRotKeyMap.end()) {
			addLeftRotKey(secretKey, ix, 0);
		}
	}
}

void Scheme::addLeftYRotKeys(SecretKey& secretKey) {
	for (long iy = 1; iy < ring.N1; iy <<=1) {
		if(leftRotKeyMap.find({0, iy}) == leftRotKeyMap.end()) {
			addLeftRotKey(secretKey, 0, iy);
		}
	}
}

void Scheme::addRightXRotKeys(SecretKey& secretKey) {
	for (long ix = 1; ix < ring.N0h; ix <<=1) {
		long idx = ring.N0h - ix;
		if(leftRotKeyMap.find({idx, 0}) == leftRotKeyMap.end()) {
			addLeftRotKey(secretKey, idx, 0);
		}
	}
}

void Scheme::addRightYRotKeys(SecretKey& secretKey) {
	for (long iy = 1; iy < ring.N1; iy<<=1) {
		long idx = ring.N1 - iy;
		if(leftRotKeyMap.find({0, idx}) == leftRotKeyMap.end()) {
			addLeftRotKey(secretKey, 0, idx);
		}
	}
}

void Scheme::addConjKey(SecretKey& secretKey) {
	ZZ* ex = new ZZ[ring.N];
	ZZ* ax = new ZZ[ring.N];
	ZZ* bx = new ZZ[ring.N];

	ZZ* conjsx = ring.conjugate(secretKey.sx);

	ring.leftShiftAndEqual(conjsx, ring.logQ, ring.QQ);
	ring.sampleUniform(ax, ring.logQQ);
	ring.sampleGauss(ex);
	ring.addAndEqual(ex, conjsx, ring.QQ);

	long np = ceil((1 + ring.logQQ + ring.logN + 3)/59.0);
	ring.mult(bx, secretKey.sx, ax, np, ring.QQ);
	ring.sub(bx, ex, bx, ring.QQ);

	np = ceil((2 * ring.logQQ + ring.logN + 3)/59.0);
	uint64_t* rax = ring.toNTT(ax, np);
	uint64_t* rbx = ring.toNTT(bx, np);

	keyMap.insert(pair<long, Key>(CONJUGATION, Key(rax, rbx)));

	delete[] ex;
	delete[] ax;
	delete[] bx;
	delete[] conjsx;
}

void Scheme::addBootKey(SecretKey& secretKey, long lognx, long logny, long logp) {
	ring.addBootContext(lognx, logny, logp);

	addConjKey(secretKey);
	addLeftXRotKeys(secretKey);
	addLeftYRotKeys(secretKey);

	long lognxh = lognx / 2;
	long kx = 1 << lognxh;
	long mx = 1 << (lognx - lognxh);

	for (long i = 1; i < kx; ++i) {
		if(leftRotKeyMap.find({i,0}) == leftRotKeyMap.end()) {
			addLeftRotKey(secretKey, i, 0);
		}
	}

	for (long i = 1; i < mx; ++i) {
		long idx = i * kx;
		if(leftRotKeyMap.find({idx, 0}) == leftRotKeyMap.end()) {
			addLeftRotKey(secretKey, idx, 0);
		}
	}

	long lognyh = logny / 2;
	long ky = 1 << lognyh;
	long my = 1 << (logny - lognyh);

	for (long i = 1; i < ky; ++i) {
		if(leftRotKeyMap.find({0,i}) == leftRotKeyMap.end()) {
			addLeftRotKey(secretKey, 0, i);
		}
	}

	for (long i = 1; i < my; ++i) {
		long idx = i * ky;
		if(leftRotKeyMap.find({0, idx}) == leftRotKeyMap.end()) {
			addLeftRotKey(secretKey, 0, idx);
		}
	}
}

void Scheme::addSquareMatrixKeys(SecretKey& secretKey, long lognx) {
	ring.addMatrixContext(lognx);
	long nx = 1 << lognx;
	for (long ix = 1; ix < nx; ++ix) {
		long idx = ring.N0h - ix;
		if(leftRotKeyMap.find({idx, 0}) == leftRotKeyMap.end()) {
			addLeftRotKey(secretKey, idx, 0);
		}
	}
	addLeftYRotKeys(secretKey);
}


//----------------------------------------------------------------------------------
//   ENCODING & DECODING
//----------------------------------------------------------------------------------


Plaintext Scheme::encode(complex<double>* vals, long nx, long ny, long logp, long logq) {
	ZZ* mx = new ZZ[ring.N];
	ring.encode(mx, vals, nx, ny, logp + ring.logQ);
	return Plaintext(mx, logp, logq, ring.N0, ring.N1, nx, ny);
}

Plaintext Scheme::encode(double* vals, long nx, long ny, long logp, long logq) {
	ZZ* mxy = new ZZ[ring.N];
	ring.encode(mxy, vals, nx, ny, logp + ring.logQ);
	return Plaintext(mxy, logp, logq, ring.N0, ring.N1, nx, ring.N1);
}

Plaintext Scheme::encodeSingle(complex<double> val, long logp, long logq) {
	ZZ* mxy = new ZZ[ring.N];

	mxy[0] = EvaluatorUtils::scaleUpToZZ(val.real(), logp + ring.logQ);
	mxy[ring.Nh] = EvaluatorUtils::scaleUpToZZ(val.imag(), logp + ring.logQ);

	return Plaintext(mxy, logp, logq, ring.N0, ring.N1, 1, 1);
}

Plaintext Scheme::encodeSingle(double val, long logp, long logq) {
	ZZ* mx = new ZZ[ring.N];

	mx[0] = EvaluatorUtils::scaleUpToZZ(val, logp + ring.logQ);

	return Plaintext(mx, logp, logq, ring.N0, ring.N1, 1, 1);
}

complex<double>* Scheme::decode(Plaintext& msg) {
	complex<double>* vals = new complex<double>[msg.n0 * msg.n1];
	ring.decode(msg.mx, vals, msg.n0, msg.n1, msg.logp, msg.logq);
	return vals;
}

complex<double> Scheme::decodeSingle(Plaintext& msg) {
	complex<double> res;
	ZZ tmp;
	ZZ q = ring.qvec[msg.logq];
	ZZ qh = ring.qvec[msg.logq - 1];
	tmp = msg.mx[0];
	while(tmp < 0) tmp += q;
	while(tmp > qh) tmp -= q;
	res.real(EvaluatorUtils::scaleDownToReal(tmp, msg.logp));

	tmp = msg.mx[ring.Nh];
	while(tmp < 0) tmp += q;
	while(tmp > qh) tmp -= q;
	res.imag(EvaluatorUtils::scaleDownToReal(tmp, msg.logp));
	return res;
}


//----------------------------------------------------------------------------------
//   ENCRYPTION & DECRYPTION
//----------------------------------------------------------------------------------


Ciphertext Scheme::encryptMsg(Plaintext& msg) {
	ZZ qQ = ring.qvec[msg.logq + ring.logQ];

	ZZ* exy = new ZZ[ring.N];
	ZZ* vxy = new ZZ[ring.N];
	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.sampleZO(vxy);

	Key key = keyMap.at(ENCRYPTION);

	long np = ceil((1 + ring.logQQ + ring.logN + 3)/59.0);
	ring.multNTT(axy, vxy, key.rax, np, qQ);
	ring.sampleGauss(exy);
	ring.addAndEqual(axy, exy, qQ);

	ring.multNTT(bxy, vxy, key.rbx, np, qQ);
	ring.sampleGauss(exy);
	ring.addAndEqual(bxy, exy, qQ);

	ring.addAndEqual(bxy, msg.mx, qQ);
	ring.rightShiftAndEqual(axy, ring.logQ);
	ring.rightShiftAndEqual(bxy, ring.logQ);

	delete[] exy;
	delete[] vxy;

	return Ciphertext(axy, bxy, msg.logp, msg.logq, msg.N0, msg.N1, msg.n0, msg.n1);
}

Ciphertext Scheme::encrypt(complex<double>* vals, long n0, long n1, long logp, long logq) {
	Plaintext msg = encode(vals, n0, n1, logp, logq);
	return encryptMsg(msg);
}

Ciphertext Scheme::encrypt(double* vals, long n0, long n1, long logp, long logq) {
	Plaintext msg = encode(vals, n0, n1, logp, logq);
	return encryptMsg(msg);
}

Ciphertext Scheme::encryptSingle(complex<double> val, long logp, long logq) {
	Plaintext msg = encodeSingle(val, logp, logq);
	return encryptMsg(msg);
}

Ciphertext Scheme::encryptSingle(double val, long logp, long logq) {
	Plaintext msg = encodeSingle(val, logp, logq);
	return encryptMsg(msg);
}

Ciphertext Scheme::encryptZeros(long nx, long ny, long logp, long logq) {
	Ciphertext czeros = encryptSingle(0.0, logp, logq);
	czeros.n0 = nx;
	czeros.n1 = ny;
	return czeros;
}

Plaintext Scheme::decryptMsg(SecretKey& secretKey, Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ* mx = new ZZ[ring.N];

	long np = ceil((1 + cipher.logq + ring.logN + 3)/59.0);
	ring.mult(mx, cipher.ax, secretKey.sx, np, q);
	ring.addAndEqual(mx, cipher.bx, q);
	return Plaintext(mx, cipher.logp, cipher.logq, cipher.N0, cipher.N1, cipher.n0, cipher.n1);
}

complex<double>* Scheme::decrypt(SecretKey& secretKey, Ciphertext& cipher) {
	Plaintext msg = decryptMsg(secretKey, cipher);
	return decode(msg);
}

complex<double> Scheme::decryptSingle(SecretKey& secretKey, Ciphertext& cipher) {
	Plaintext msg = decryptMsg(secretKey, cipher);
	return decodeSingle(msg);
}


//----------------------------------------------------------------------------------
//   HOMOMORPHIC OPERATIONS
//----------------------------------------------------------------------------------


Ciphertext Scheme::negate(Ciphertext& cipher) {
	ZZ* ax = new ZZ[ring.N];
	ZZ* bx = new ZZ[ring.N];

	ring.negate(ax, cipher.ax);
	ring.negate(bx, cipher.bx);

	return Ciphertext(ax, bx, cipher.logp, cipher.logq, cipher.N0, cipher.N1, cipher.n0, cipher.n1);
}

void Scheme::negateAndEqual(Ciphertext& cipher) {
	ring.negateAndEqual(cipher.ax);
	ring.negateAndEqual(cipher.bx);
}

Ciphertext Scheme::add(Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZ q = ring.qvec[cipher1.logq];

	ZZ* ax = new ZZ[ring.N];
	ZZ* bx = new ZZ[ring.N];

	ring.add(ax, cipher1.ax, cipher2.ax, q);
	ring.add(bx, cipher1.bx, cipher2.bx, q);

	return Ciphertext(ax, bx, cipher1.logp, cipher1.logq, cipher1.N0, cipher1.N1, cipher1.n0, cipher1.n1);
}

void Scheme::addAndEqual(Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZ q = ring.qvec[cipher1.logq];

	ring.addAndEqual(cipher1.ax, cipher2.ax, q);
	ring.addAndEqual(cipher1.bx, cipher2.bx, q);
}

Ciphertext Scheme::addConst(Ciphertext& cipher, double cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	Ciphertext res = cipher;
	ZZ cnstZZ = logp < 0 ? EvaluatorUtils::scaleUpToZZ(cnst, cipher.logp) : EvaluatorUtils::scaleUpToZZ(cnst, logp);
	AddMod(res.bx[0], res.bx[0], cnstZZ, q);
	return res;
}

Ciphertext Scheme::addConst(Ciphertext& cipher, RR& cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	Ciphertext res = cipher;
	ZZ cnstZZ = logp < 0 ? EvaluatorUtils::scaleUpToZZ(cnst, cipher.logp) : EvaluatorUtils::scaleUpToZZ(cnst, logp);
	AddMod(res.bx[0], res.bx[0], cnstZZ, q);
	return res;
}

void Scheme::addConstAndEqual(Ciphertext& cipher, double cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ cnstZZ = logp < 0 ? EvaluatorUtils::scaleUpToZZ(cnst, cipher.logp) : EvaluatorUtils::scaleUpToZZ(cnst, logp);

	AddMod(cipher.bx[0], cipher.bx[0], cnstZZ, q);
}

void Scheme::addConstAndEqual(Ciphertext& cipher, RR& cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ cnstZZ = logp < 0 ? EvaluatorUtils::scaleUpToZZ(cnst, cipher.logp) : EvaluatorUtils::scaleUpToZZ(cnst, logp);

	AddMod(cipher.bx[0], cipher.bx[0], cnstZZ, q);
}

Ciphertext Scheme::addPoly(Ciphertext& cipher, ZZ* poly, long logp) {
	ZZ q = ring.qvec[cipher.logq];
	Ciphertext res = cipher;
	ring.addAndEqual(res.bx, poly, q);
	return res;
}

void Scheme::addPolyAndEqual(Ciphertext& cipher, ZZ* poly, long logp) {
	ZZ q = ring.qvec[cipher.logq];
	ring.addAndEqual(cipher.bx, poly, q);
}

Ciphertext Scheme::sub(Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZ q = ring.qvec[cipher1.logq];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.sub(axy, cipher1.ax, cipher2.ax, q);
	ring.sub(bxy, cipher1.bx, cipher2.bx, q);

	return Ciphertext(axy, bxy, cipher1.logp, cipher1.logq, cipher1.N0, cipher1.N1, cipher1.n0, cipher1.n1);
}

void Scheme::subAndEqual(Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZ q = ring.qvec[cipher1.logq];

	ring.subAndEqual(cipher1.ax, cipher2.ax, q);
	ring.subAndEqual(cipher1.bx, cipher2.bx, q);
}

void Scheme::subAndEqual2(Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZ q = ring.qvec[cipher1.logq];

	ring.subAndEqual2(cipher1.ax, cipher2.ax, q);
	ring.subAndEqual2(cipher1.bx, cipher2.bx, q);
}

Ciphertext Scheme::imult(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ* axy = ring.multByMonomial(cipher.ax, ring.N0h, 0, q);
	ZZ* bxy = ring.multByMonomial(cipher.bx, ring.N0h, 0, q);

	return Ciphertext(axy, bxy, cipher.logp, cipher.logq, cipher.N0, cipher.N1, cipher.n0, cipher.n1);
}

Ciphertext Scheme::idiv(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ* axy = ring.multByMonomial(cipher.ax, 3 * ring.N0h, 0, q);
	ZZ* bxy = ring.multByMonomial(cipher.bx, 3 * ring.N0h, 0, q);

	return Ciphertext(axy, bxy, cipher.logp, cipher.logq, cipher.N0, cipher.N1, cipher.n0, cipher.n1);
}

void Scheme::imultAndEqual(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];

	ring.multByMonomialAndEqual(cipher.ax, ring.N0h, 0, q);
	ring.multByMonomialAndEqual(cipher.bx, ring.N0h, 0, q);
}

void Scheme::idivAndEqual(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	ring.multByMonomialAndEqual(cipher.ax, 3 * ring.N0h, 0, q);
	ring.multByMonomialAndEqual(cipher.bx, 3 * ring.N0h, 0, q);
}

Ciphertext Scheme::mult(Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZ q = ring.qvec[cipher1.logq];
	ZZ qQ = ring.qvec[cipher1.logq + ring.logQ];

	long np = ceil((2 + cipher1.logq + cipher2.logq + ring.logN + 3)/59.0);
	uint64_t* ra1 = ring.toNTT(cipher1.ax, np);
	uint64_t* rb1 = ring.toNTT(cipher1.bx, np);

	uint64_t* ra2 = ring.toNTT(cipher2.ax, np);
	uint64_t* rb2 = ring.toNTT(cipher2.bx, np);

	uint64_t* rab1 = ring.addNTT(ra1, rb1, np);
	uint64_t* rab2 = ring.addNTT(ra2, rb2, np);

	ZZ* abx = new ZZ[ring.N];
	ring.multDNTT(abx, rab1, rab2, np, q);

	ZZ* aax = new ZZ[ring.N];
	ZZ* bbx = new ZZ[ring.N];
	np = ceil((cipher1.logq + cipher2.logq + ring.logN + 3)/59.0);
	ring.multDNTT(aax, ra1, ra2, np, q);
	ring.multDNTT(bbx, rb1, rb2, np, q);

	Key key = keyMap.at(MULTIPLICATION);

	ZZ* ax = new ZZ[ring.N];
	ZZ* bx = new ZZ[ring.N];

	np = ceil((cipher1.logq + ring.logQQ + ring.logN + 3)/59.0);
	uint64_t* raa = ring.toNTT(aax, np);
	ring.multDNTT(ax, raa, key.rax, np, qQ);
	ring.multDNTT(bx, raa, key.rbx, np, qQ);

	ring.rightShiftAndEqual(ax, ring.logQ);
	ring.rightShiftAndEqual(bx, ring.logQ);

	ring.addAndEqual(ax, abx, q);
	ring.subAndEqual(ax, bbx, q);
	ring.subAndEqual(ax, aax, q);
	ring.addAndEqual(bx, bbx, q);

	delete[] abx;
	delete[] aax;
	delete[] bbx;

	delete[] ra1;
	delete[] ra2;
	delete[] rb1;
	delete[] rb2;
	delete[] rab1;
	delete[] rab2;
	delete[] raa;

	return Ciphertext(ax, bx, cipher1.logp + cipher2.logp, cipher1.logq, cipher1.N0, cipher1.N1, cipher1.n0, cipher1.n1);
}

void Scheme::multAndEqual(Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZ q = ring.qvec[cipher1.logq];
	ZZ qQ = ring.qvec[cipher1.logq + ring.logQ];

	long np = ceil((2 + cipher1.logq + cipher2.logq + ring.logN + 3)/59.0);

	uint64_t* ra1 = ring.toNTT(cipher1.ax, np);
	uint64_t* rb1 = ring.toNTT(cipher1.bx, np);

	uint64_t* ra2 = ring.toNTT(cipher2.ax, np);
	uint64_t* rb2 = ring.toNTT(cipher2.bx, np);

	uint64_t* rab1 = ring.addNTT(ra1, rb1, np);
	uint64_t* rab2 = ring.addNTT(ra2, rb2, np);


	ZZ* abx = new ZZ[ring.N];
	ring.multDNTT(abx, rab1, rab2, np, q);

	ZZ* aax = new ZZ[ring.N];
	ZZ* bbx = new ZZ[ring.N];
	np = ceil((cipher1.logq + cipher2.logq + ring.logN + 3)/59.0);
	ring.multDNTT(aax, ra1, ra2, np, q);
	ring.multDNTT(bbx, rb1, rb2, np, q);

	Key key = keyMap.at(MULTIPLICATION);

	np = ceil((cipher1.logq + ring.logQQ + ring.logN + 3)/59.0);
	uint64_t* raa = ring.toNTT(aax, np);
	ring.multDNTT(cipher1.ax, raa, key.rax, np, qQ);
	ring.multDNTT(cipher1.bx, raa, key.rbx, np, qQ);

	ring.rightShiftAndEqual(cipher1.ax, ring.logQ);
	ring.rightShiftAndEqual(cipher1.bx, ring.logQ);

	ring.addAndEqual(cipher1.ax, abx, q);
	ring.subAndEqual(cipher1.ax, bbx, q);
	ring.subAndEqual(cipher1.ax, aax, q);
	ring.addAndEqual(cipher1.bx, bbx, q);

	cipher1.logp += cipher2.logp;
	delete[] abx;
	delete[] aax;
	delete[] bbx;

	delete[] ra1;
	delete[] ra2;
	delete[] rb1;
	delete[] rb2;
	delete[] rab1;
	delete[] rab2;
	delete[] raa;
}

Ciphertext Scheme::square(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ qQ = ring.qvec[cipher.logq + ring.logQ];

	long np = ceil((2 * cipher.logq + ring.logN + 3)/59.0);
	uint64_t* ra = ring.toNTT(cipher.ax, np);
	uint64_t* rb = ring.toNTT(cipher.bx, np);

	ZZ* aax = new ZZ[ring.N];
	ZZ* bbx = new ZZ[ring.N];
	ZZ* abx = new ZZ[ring.N];

	ring.squareNTT(bbx, rb, np, q);
	ring.multDNTT(abx, ra, rb, np, q);
	ring.addAndEqual(abx, abx, q);
	ring.squareNTT(aax, ra, np, q);

	Key key = keyMap.at(MULTIPLICATION);

	ZZ* ax = new ZZ[ring.N];
	ZZ* bx = new ZZ[ring.N];

	np = ceil((cipher.logq + ring.logQQ + ring.logN + 3)/59.0);
	uint64_t* raa = ring.toNTT(aax, np);
	ring.multDNTT(ax, raa, key.rax, np, qQ);
	ring.multDNTT(bx, raa, key.rbx, np, qQ);

	ring.rightShiftAndEqual(ax, ring.logQ);
	ring.rightShiftAndEqual(bx, ring.logQ);

	ring.addAndEqual(ax, abx, q);
	ring.addAndEqual(bx, bbx, q);

	delete[] abx;
	delete[] aax;
	delete[] bbx;

	delete[] ra;
	delete[] rb;
	delete[] raa;

	return Ciphertext(ax, bx, 2 * cipher.logp, cipher.logq, cipher.N0, cipher.N1, cipher.n0, cipher.n1);
}

void Scheme::squareAndEqual(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ qQ = ring.qvec[cipher.logq + ring.logQ];

	long np = ceil((2 * cipher.logq + ring.logN + 3)/59.0);
	uint64_t* ra = ring.toNTT(cipher.ax, np);
	uint64_t* rb = ring.toNTT(cipher.bx, np);

	ZZ* aax = new ZZ[ring.N];
	ZZ* bbx = new ZZ[ring.N];
	ZZ* abx = new ZZ[ring.N];

	ring.squareNTT(bbx, rb, np, q);
	ring.multDNTT(abx, ra, rb, np, q);
	ring.addAndEqual(abx, abx, q);
	ring.squareNTT(aax, ra, np, q);

	Key key = keyMap.at(MULTIPLICATION);

	np = ceil((cipher.logq + ring.logQQ + ring.logN + 3)/59.0);
	uint64_t* raa = ring.toNTT(aax, np);
	ring.multDNTT(cipher.ax, raa, key.rax, np, qQ);
	ring.multDNTT(cipher.bx, raa, key.rbx, np, qQ);

	ring.rightShiftAndEqual(cipher.ax, ring.logQ);
	ring.rightShiftAndEqual(cipher.bx, ring.logQ);

	ring.addAndEqual(cipher.ax, abx, q);
	ring.addAndEqual(cipher.bx, bbx, q);

	cipher.logp *= 2;

	delete[] abx;
	delete[] aax;
	delete[] bbx;

	delete[] ra;
	delete[] rb;
	delete[] raa;
}

Ciphertext Scheme::multConst(Ciphertext& cipher, RR& cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ZZ cnstZZ = EvaluatorUtils::scaleUpToZZ(cnst, logp);

	ring.multByConst(axy, cipher.ax, cnstZZ, q);
	ring.multByConst(bxy, cipher.bx, cnstZZ, q);

	return Ciphertext(axy, bxy, cipher.logp + logp, cipher.logq, cipher.N0, cipher.N1, cipher.n0, cipher.n1);
}


Ciphertext Scheme::multConst(Ciphertext& cipher, double cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ZZ cnstZZ = EvaluatorUtils::scaleUpToZZ(cnst, logp);

	ring.multByConst(axy, cipher.ax, cnstZZ, q);
	ring.multByConst(bxy, cipher.bx, cnstZZ, q);

	return Ciphertext(axy, bxy, cipher.logp + logp, cipher.logq, cipher.N0, cipher.N1, cipher.n0, cipher.n1);
}

Ciphertext Scheme::multConst(Ciphertext& cipher, complex<double> cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ZZ* axyi = ring.multByMonomial(cipher.ax, ring.N0h, 0, q);
	ZZ* bxyi = ring.multByMonomial(cipher.bx, ring.N0h, 0, q);

	ZZ cnstrZZ = EvaluatorUtils::scaleUpToZZ(cnst.real(), logp);
	ZZ cnstiZZ = EvaluatorUtils::scaleUpToZZ(cnst.imag(), logp);

	ring.multByConst(axy, cipher.ax, cnstrZZ, q);
	ring.multByConst(bxy, cipher.bx, cnstrZZ, q);

	ring.multByConstAndEqual(axyi, cnstiZZ, q);
	ring.multByConstAndEqual(bxyi, cnstiZZ, q);

	ring.addAndEqual(axy, axyi, q);
	ring.addAndEqual(bxy, bxyi, q);

	return Ciphertext(axy, bxy, cipher.logp + logp, cipher.logq, cipher.N0, cipher.N1, cipher.n0, cipher.n1);
}

void Scheme::multConstAndEqual(Ciphertext& cipher, RR& cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ cnstZZ = EvaluatorUtils::scaleUpToZZ(cnst, logp);

	ring.multByConstAndEqual(cipher.ax, cnstZZ, q);
	ring.multByConstAndEqual(cipher.bx, cnstZZ, q);

	cipher.logp += logp;
}

void Scheme::multConstAndEqual(Ciphertext& cipher, double cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ cnstZZ = EvaluatorUtils::scaleUpToZZ(cnst, logp);

	ring.multByConstAndEqual(cipher.ax, cnstZZ, q);
	ring.multByConstAndEqual(cipher.bx, cnstZZ, q);

	cipher.logp += logp;
}

void Scheme::multConstAndEqual(Ciphertext& cipher, complex<double> cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];


	ZZ cnstrZZ = EvaluatorUtils::scaleUpToZZ(cnst.real(), logp);
	ZZ cnstiZZ = EvaluatorUtils::scaleUpToZZ(cnst.imag(), logp);

	ZZ* axyi = ring.multByMonomial(cipher.ax, ring.N0h, 0, q);
	ZZ* bxyi = ring.multByMonomial(cipher.bx, ring.N0h, 0, q);

	ring.multByConstAndEqual(axyi, cnstiZZ, q);
	ring.multByConstAndEqual(bxyi, cnstiZZ, q);

	ring.multByConstAndEqual(cipher.ax, cnstrZZ, q);
	ring.multByConstAndEqual(cipher.bx, cnstrZZ, q);

	ring.addAndEqual(cipher.ax, axyi, q);
	ring.addAndEqual(cipher.bx, bxyi, q);

	cipher.logp += logp;
}


Ciphertext Scheme::multPolyX0(Ciphertext& cipher, ZZ* poly, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ* ax = new ZZ[ring.N];
	ZZ* bx = new ZZ[ring.N];

	long np = ceil((cipher.logq + logp + ring.logN0 + 3)/59.0);
	uint64_t* rpoly = ring.toNTTX0(poly, np);
	ring.multNTTX0(ax, cipher.ax, rpoly, np, q);
	ring.multNTTX0(bx, cipher.bx, rpoly, np, q);
	delete[] rpoly;

	return Ciphertext(ax, bx, cipher.logp + logp, cipher.logq, cipher.N0, cipher.N1, cipher.n0, cipher.n1);
}

void Scheme::multPolyX0AndEqual(Ciphertext& cipher, ZZ* poly, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	long np = ceil((cipher.logq + logp + ring.logN0 + 3)/59.0);
	uint64_t* rpoly = ring.toNTTX0(poly, np);
	ring.multNTTX0AndEqual(cipher.ax, rpoly, np, q);
	ring.multNTTX0AndEqual(cipher.bx, rpoly, np, q);
	delete[] rpoly;

	cipher.logp += logp;
}

Ciphertext Scheme::multPolyNTTX0(Ciphertext& cipher, uint64_t* rpoly, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ* ax = new ZZ[ring.N];
	ZZ* bx = new ZZ[ring.N];

	long np = ceil((cipher.logq + logp + ring.logN0 + 3)/59.0);
	ring.multNTTX0(ax, cipher.ax, rpoly, np, q);
	ring.multNTTX0(bx, cipher.bx, rpoly, np, q);
	return Ciphertext(ax, bx, cipher.logp + logp, cipher.logq, cipher.N0, cipher.N1, cipher.n0, cipher.n1);

}

void Scheme::multPolyNTTX0AndEqual(Ciphertext& cipher, uint64_t* rpoly, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	long np = ceil((cipher.logq + logp + ring.logN0 + 3)/59.0);
	ring.multNTTX0AndEqual(cipher.ax, rpoly, np, q);
	ring.multNTTX0AndEqual(cipher.bx, rpoly, np, q);
	cipher.logp += logp;
}

Ciphertext Scheme::multPolyX1(Ciphertext& cipher, ZZ* rpoly, ZZ* ipoly, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ* axyr = new ZZ[ring.N];
	ZZ* bxyr = new ZZ[ring.N];

	ZZ* axyi = ring.multByMonomial(cipher.ax, ring.N0h, 0, q);
	ZZ* bxyi = ring.multByMonomial(cipher.bx, ring.N0h, 0, q);

	long np = ceil((cipher.logq + logp + ring.logN1 + 3)/59.0);
	uint64_t* ripoly = ring.toNTT(ipoly, np);
	ring.multNTTX1AndEqual(axyi, ripoly, np, q);
	ring.multNTTX1AndEqual(bxyi, ripoly, np, q);

	np = ceil((cipher.logq + logp + ring.logN1 + 3)/59.0);
	uint64_t* rrpoly = ring.toNTT(rpoly, np);
	ring.multNTTX1(axyr, cipher.ax, rrpoly, np, q);
	ring.multNTTX1(bxyr, cipher.bx, rrpoly, np, q);

	ring.addAndEqual(axyr, axyi, q);
	ring.addAndEqual(bxyr, bxyi, q);

	delete[] axyi;
	delete[] bxyi;
	delete[] ripoly;
	delete[] rrpoly;

	return Ciphertext(axyr, bxyr, cipher.logp + logp, cipher.logq, cipher.N0, cipher.N1, cipher.n0, cipher.n1);
}

void Scheme::multPolyX1AndEqual(Ciphertext& cipher, ZZ* rpoly, ZZ* ipoly, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ* axyi = ring.multByMonomial(cipher.ax, ring.N0h, 0, q);
	ZZ* bxyi = ring.multByMonomial(cipher.bx, ring.N0h, 0, q);

	long np = ceil((cipher.logq + logp + ring.logN1 + 3)/59.0);
	uint64_t* rrpoly = ring.toNTT(rpoly, np);
	ring.multNTTX1AndEqual(cipher.ax, rrpoly, np, q);
	ring.multNTTX1AndEqual(cipher.bx, rrpoly, np, q);

	np = ceil((cipher.logq + logp + ring.logN1 + 3)/59.0);
	uint64_t* ripoly = ring.toNTT(ipoly, np);
	ring.multNTTX1AndEqual(axyi, ripoly, np, q);
	ring.multNTTX1AndEqual(bxyi, ripoly, np, q);

	ring.addAndEqual(cipher.ax, axyi, q);
	ring.addAndEqual(cipher.bx, bxyi, q);

	delete[] axyi;
	delete[] bxyi;
	delete[] ripoly;
	delete[] rrpoly;

	cipher.logp += logp;
}

Ciphertext Scheme::multPoly(Ciphertext& cipher, ZZ* poly, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	long np = ceil((cipher.logq + logp + ring.logN + 3)/59.0);
	uint64_t* rpoly = ring.toNTT(poly, np);
	ring.multNTT(axy, cipher.ax, rpoly, np, q);
	ring.multNTT(bxy, cipher.bx, rpoly, np, q);
	delete[] rpoly;
	return Ciphertext(axy, bxy, cipher.logp + logp, cipher.logq, cipher.N0, cipher.N1, cipher.n0, cipher.n1);
}

void Scheme::multPolyAndEqual(Ciphertext& cipher, ZZ* poly, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	long np = ceil((cipher.logq + logp + ring.logN + 3)/59.0);
	uint64_t* rpoly = ring.toNTT(poly, np);
	ring.multNTTAndEqual(cipher.ax, rpoly, np, q);
	ring.multNTTAndEqual(cipher.bx, rpoly, np, q);
	delete[] rpoly;

	cipher.logp += logp;
}


Ciphertext Scheme::multPolyNTT(Ciphertext& cipher, uint64_t* rpoly, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	long np = ceil((cipher.logq + logp + ring.logN + 3)/59.0);
	ring.multNTT(axy, cipher.ax, rpoly, np, q);
	ring.multNTT(bxy, cipher.bx, rpoly, np, q);

	return Ciphertext(axy, bxy, cipher.logp + logp, cipher.logq, cipher.N0, cipher.N1, cipher.n0, cipher.n1);
}

void Scheme::multPolyNTTAndEqual(Ciphertext& cipher, uint64_t* rpoly, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	long np = ceil((cipher.logq + logp + ring.logN + 3)/59.0);
	ring.multNTTAndEqual(cipher.ax, rpoly, np, q);
	ring.multNTTAndEqual(cipher.bx, rpoly, np, q);

	cipher.logp += logp;
}

Ciphertext Scheme::multByMonomial(Ciphertext& cipher, const long dx, const long dy) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ* axy = ring.multByMonomial(cipher.ax, dx, dy, q);
	ZZ* bxy = ring.multByMonomial(cipher.bx, dx, dy, q);

	return Ciphertext(axy, bxy, cipher.logp, cipher.logq, cipher.N0, cipher.N1, cipher.n0, cipher.n1);
}

void Scheme::multByMonomialAndEqual(Ciphertext& cipher, const long dx, const long dy) {
	ZZ q = ring.qvec[cipher.logq];

	ring.multByMonomialAndEqual(cipher.ax, dx, dy, q);
	ring.multByMonomialAndEqual(cipher.bx, dx, dy, q);
}

Ciphertext Scheme::multPo2(Ciphertext& cipher, long degree) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.leftShift(axy, cipher.ax, degree, q);
	ring.leftShift(bxy, cipher.bx, degree, q);

	return Ciphertext(axy, bxy, cipher.logp, cipher.logq, cipher.N0, cipher.N1, cipher.n0, cipher.n1);
}

void Scheme::multPo2AndEqual(Ciphertext& cipher, long degree) {
	ZZ q = ring.qvec[cipher.logq];

	ring.leftShiftAndEqual(cipher.ax, degree, q);
	ring.leftShiftAndEqual(cipher.bx, degree, q);
}

void Scheme::doubleAndEqual(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];

	ring.doubleAndEqual(cipher.ax, q);
	ring.doubleAndEqual(cipher.bx, q);
}

Ciphertext Scheme::divPo2(Ciphertext& cipher, long logd) {
	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.rightShift(axy, cipher.ax, logd);
	ring.rightShift(bxy, cipher.bx, logd);

	return Ciphertext(axy, bxy, cipher.logp, cipher.logq - logd, cipher.N0, cipher.N1, cipher.n0, cipher.n1);
}

void Scheme::divPo2AndEqual(Ciphertext& cipher, long logd) {
	ring.rightShiftAndEqual(cipher.ax, logd);
	ring.rightShiftAndEqual(cipher.bx, logd);
	cipher.logq -= logd;
}


//----------------------------------------------------------------------------------
//   RESCALING & MODULUS DOWN
//----------------------------------------------------------------------------------


Ciphertext Scheme::reScaleBy(Ciphertext& cipher, long dlogq) {
	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.rightShift(axy, cipher.ax, dlogq);
	ring.rightShift(bxy, cipher.bx, dlogq);

	return Ciphertext(axy, bxy, cipher.logp - dlogq, cipher.logq - dlogq, cipher.N0, cipher.N1, cipher.n0, cipher.n1);
}

Ciphertext Scheme::reScaleTo(Ciphertext& cipher, long logq) {
	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	long dlogq = cipher.logq - logq;
	ring.rightShift(axy, cipher.ax, dlogq);
	ring.rightShift(bxy, cipher.bx, dlogq);

	return Ciphertext(axy, bxy, cipher.logp - dlogq, logq, cipher.N0, cipher.N1, cipher.n0, cipher.n1);
}

void Scheme::reScaleByAndEqual(Ciphertext& cipher, long dlogq) {
	ring.rightShiftAndEqual(cipher.ax, dlogq);
	ring.rightShiftAndEqual(cipher.bx, dlogq);
	cipher.logq -= dlogq;
	cipher.logp -= dlogq;
}

void Scheme::reScaleToAndEqual(Ciphertext& cipher, long logq) {
	long dlogq = cipher.logq - logq;
	ring.rightShiftAndEqual(cipher.ax, dlogq);
	ring.rightShiftAndEqual(cipher.bx, dlogq);
	cipher.logq = logq;
	cipher.logp -= dlogq;
}

Ciphertext Scheme::modDownBy(Ciphertext& cipher, long dlogq) {
	long logq = cipher.logq - dlogq;
	ZZ q = ring.qvec[logq];
	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.mod(axy, cipher.ax, q);
	ring.mod(bxy, cipher.bx, q);
	return Ciphertext(axy, bxy, cipher.logp, logq, cipher.N0, cipher.N1, cipher.n0, cipher.n1);
}

void Scheme::modDownByAndEqual(Ciphertext& cipher, long dlogq) {
	ZZ q = ring.qvec[cipher.logq - dlogq];

	ring.modAndEqual(cipher.ax, q);
	ring.modAndEqual(cipher.bx, q);

	cipher.logq -= dlogq;
}

Ciphertext Scheme::modDownTo(Ciphertext& cipher, long logq) {
	ZZ q = ring.qvec[logq];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.mod(axy, cipher.ax, q);
	ring.mod(bxy, cipher.bx, q);

	return Ciphertext(axy, bxy, cipher.logp, logq, cipher.N0, cipher.N1, cipher.n0, cipher.n1);
}

void Scheme::modDownToAndEqual(Ciphertext& cipher, long newlogq) {
	ZZ q = ring.qvec[newlogq];

	ring.modAndEqual(cipher.ax, q);
	ring.modAndEqual(cipher.bx, q);

	cipher.logq = newlogq;
}


//----------------------------------------------------------------------------------
//   ROTATIONS & CONJUGATIONS
//----------------------------------------------------------------------------------


Ciphertext Scheme::leftRotateFast(Ciphertext& cipher, long rx, long ry) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ qQ = ring.qvec[cipher.logq + ring.logQ];

	ZZ* rxbx = ring.leftRotate(cipher.bx, rx, ry);
	ZZ* rxax = ring.leftRotate(cipher.ax, rx, ry);

	Key key = leftRotKeyMap.at({rx, ry});

	ZZ* ax = new ZZ[ring.N];
	ZZ* bx = new ZZ[ring.N];
	long np = ceil((cipher.logq + ring.logQQ + ring.logN + 3)/59.0);
	uint64_t* rrxa = ring.toNTT(rxax, np);
	ring.multDNTT(ax, rrxa, key.rax, np, qQ);
	ring.multDNTT(bx, rrxa, key.rbx, np, qQ);

	ring.rightShiftAndEqual(ax, ring.logQ);
	ring.rightShiftAndEqual(bx, ring.logQ);

	ring.addAndEqual(bx, rxbx, q);

	delete[] rxbx;
	delete[] rxax;
	delete[] rrxa;

	return Ciphertext(ax, bx, cipher.logp, cipher.logq, cipher.N0, cipher.N1, cipher.n0, cipher.n1);
}

Ciphertext Scheme::rightRotateFast(Ciphertext& cipher, long rx, long ry) {
	Ciphertext res = cipher;
	rightRotateFastAndEqual(res, rx, ry);
	return res;
}

void Scheme::leftRotateFastAndEqual(Ciphertext& cipher, long rx, long ry) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ qQ = ring.qvec[cipher.logq + ring.logQ];

	ZZ* rxax = ring.leftRotate(cipher.ax, rx, ry);
	ZZ* rxbx = ring.leftRotate(cipher.bx, rx, ry);

	Key key = leftRotKeyMap.at({rx, ry});

	long np = ceil((cipher.logq + ring.logQQ + ring.logN + 3)/59.0);
	uint64_t* rrxa = ring.toNTT(rxax, np);
	ring.multDNTT(cipher.bx, rrxa, key.rbx, np, qQ);
	ring.multDNTT(cipher.ax, rrxa, key.rax, np, qQ);

	ring.rightShiftAndEqual(cipher.ax, ring.logQ);
	ring.rightShiftAndEqual(cipher.bx, ring.logQ);

	ring.addAndEqual(cipher.bx, rxbx, q);

	delete[] rxax;
	delete[] rxbx;
	delete[] rrxa;
}

void Scheme::rightRotateFastAndEqual(Ciphertext& cipher, long rx, long ry) {
	long rrx = rx == 0 ? 0 : ring.N0h - rx;
	long rry = ry == 0 ? 0 : ring.N1 - ry;
	leftRotateFastAndEqual(cipher, rrx, rry);
}

Ciphertext Scheme::leftRotate(Ciphertext& cipher, long rx, long ry) {
	Ciphertext res = cipher;
	leftRotateAndEqual(res, rx, ry);
	return res;
}

Ciphertext Scheme::rightRotate(Ciphertext& cipher, long rx, long ry) {
	Ciphertext res = cipher;
	rightRotateAndEqual(res, rx, ry);
	return res;
}


void Scheme::leftRotateAndEqual(Ciphertext& cipher, long r0, long r1) {
	r0 %= cipher.n0;
	long logrxRem = log2((double)r0) + 1;
	long ipow;
	for (long i = 0; i < logrxRem; ++i) {
		if(bit(r0, i)) {
			ipow = 1 << i;
			leftRotateFastAndEqual(cipher, ipow, 0);
		}
	}
	long ryRem = r1 % cipher.n1;
	long logryRem = log2((double)ryRem) + 1;
	for (long i = 0; i < logryRem; ++i) {
		if(bit(ryRem, i)) {
			ipow = 1 << i;
			leftRotateFastAndEqual(cipher, 0, ipow);
		}
	}
}

void Scheme::rightRotateAndEqual(Ciphertext& cipher, long rx, long ry) {
	long rxRem = rx % cipher.n0;
	long logrxRem = log2((double)rxRem) + 1;
	long ipow;
	for (long i = 0; i < logrxRem; ++i) {
		if(bit(rxRem, i)) {
			ipow = 1 << i;
			rightRotateFastAndEqual(cipher, ipow, 0);
		}
	}
	long ryRem = ry % cipher.n1;
	long logryRem = log2((double)ryRem) + 1;
	for (long i = 0; i < logryRem; ++i) {
		if(bit(ryRem, i)) {
			ipow = 1 << i;
			rightRotateFastAndEqual(cipher, 0, ipow);
		}
	}
}

Ciphertext Scheme::conjugate(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ qQ = ring.qvec[cipher.logq + ring.logQ];

	ZZ* ax = new ZZ[ring.N];
	ZZ* bx = new ZZ[ring.N];

	ZZ* conjbx = ring.conjugate(cipher.bx);
	ZZ* conjax = ring.conjugate(cipher.ax);

	Key key = keyMap.at(CONJUGATION);

	long np = ceil((cipher.logq + ring.logQQ + ring.logN + 3)/59.0);
	uint64_t* rconja = ring.toNTT(conjax, np);
	ring.multDNTT(ax, rconja, key.rax, np, qQ);
	ring.multDNTT(bx, rconja, key.rbx, np, qQ);

	ring.rightShiftAndEqual(ax, ring.logQ);
	ring.rightShiftAndEqual(bx, ring.logQ);

	ring.addAndEqual(bx, conjbx, q);

	delete[] conjax;
	delete[] conjbx;
	delete[] rconja;

	return Ciphertext(ax, bx, cipher.logp, cipher.logq, cipher.N0, cipher.N1, cipher.n0, cipher.n1);
}

void Scheme::conjugateAndEqual(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ qQ = ring.qvec[cipher.logq + ring.logQ];

	ZZ* conjax = ring.conjugate(cipher.ax);
	ZZ* conjbx = ring.conjugate(cipher.bx);

	Key key = keyMap.at(CONJUGATION);

	long np = ceil((cipher.logq + ring.logQQ + ring.logN + 3)/59.0);
	uint64_t* rconja = ring.toNTT(conjax, np);
	ring.multDNTT(cipher.ax, rconja, key.rax, np, qQ);
	ring.multDNTT(cipher.bx, rconja, key.rbx, np, qQ);

	ring.rightShiftAndEqual(cipher.ax, ring.logQ);
	ring.rightShiftAndEqual(cipher.bx, ring.logQ);

	ring.addAndEqual(cipher.bx, conjbx, q);

	delete[] conjax;
	delete[] conjbx;
	delete[] rconja;
}


//----------------------------------------------------------------------------------
//   BOOTSTRAPPING
//----------------------------------------------------------------------------------


void Scheme::normalizeAndEqual(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ qh = ring.qvec[cipher.logq - 1];
	for (long i = 0; i < ring.N; ++i) {
		if (cipher.ax[i] > qh) {
			cipher.ax[i] -= q;
		} else if (cipher.ax[i] < -qh) {
			cipher.ax[i] += q;
		}
		if (cipher.bx[i] > qh) {
			cipher.bx[i] -= q;
		} else if (cipher.bx[i] < -qh) {
			cipher.bx[i] += q;
		}
	}
}

void Scheme::coeffToSlotX0AndEqual(Ciphertext& cipher) {
	long nx = cipher.n0;
	long ny = cipher.n1;

	long lognx = log2(nx);
	long logny = log2(ny);

	long logkx = lognx / 2;
	long kx = 1 << logkx;

	Ciphertext* rotvec = new Ciphertext[kx];

	NTL_EXEC_RANGE(kx, first, last);
	for (long j = first; j < last; ++j) {
		rotvec[j] = (j == 0) ? cipher : leftRotateFast(cipher, j, 0);
	}
	NTL_EXEC_RANGE_END;

	BootContext bootContext = ring.bootContextMap.at({lognx, logny});

	Ciphertext* tmpvec = new Ciphertext[kx];

	NTL_EXEC_RANGE(kx, first, last);
	for (long j = first; j < last; ++j) {
		tmpvec[j] = multPolyNTTX0(rotvec[j], bootContext.rpxVec[j], bootContext.logp);
	}
	NTL_EXEC_RANGE_END;

	for (long j = 1; j < kx; ++j) {
		addAndEqual(tmpvec[0], tmpvec[j]);
	}
	cipher = tmpvec[0];

	for (long ki = kx; ki < nx; ki += kx) {
		NTL_EXEC_RANGE(kx, first, last);
		for (long j = first; j < last; ++j) {
			tmpvec[j] = multPolyNTTX0(rotvec[j], bootContext.rpxVec[j + ki], bootContext.logp);
		}
		NTL_EXEC_RANGE_END;
		for (long j = 1; j < kx; ++j) {
			addAndEqual(tmpvec[0], tmpvec[j]);
		}
		leftRotateFastAndEqual(tmpvec[0], ki, 0);
		addAndEqual(cipher, tmpvec[0]);
	}
	reScaleByAndEqual(cipher, bootContext.logp);
	delete[] rotvec;
	delete[] tmpvec;
}

void Scheme::coeffToSlotX1AndEqual(Ciphertext& cipher) {
	long n0 = cipher.n0;
	long n1 = cipher.n1;
	long logn0 = log2(n0);
	long logn1 = log2(n1);

	long logk1 = logn1 / 2;
	long k1 = 1 << logk1;

	BootContext bootContext = ring.bootContextMap.at({logn0, logn1});

	multByMonomialAndEqual(cipher, 0, 1);
	Ciphertext rot = cipher;
	for (long i = 1; i < n1; i <<= 1) {
		Ciphertext tmp = leftRotateFast(rot, 0, i);
		addAndEqual(rot, tmp);
	}

	Ciphertext* rotvec = new Ciphertext[k1];

	NTL_EXEC_RANGE(k1, first, last);
	for (long j = first; j < last; ++j) {
		rotvec[j] = (j == 0) ? cipher : leftRotateFast(cipher, 0, j);
	}
	NTL_EXEC_RANGE_END;

	Ciphertext* tmpvec = new Ciphertext[k1];

	NTL_EXEC_RANGE(k1, first, last);
	for (long j = first; j < last; ++j) {
		complex<double> cnst = ring.ksiM1Pows[ring.M1 - ring.gM1Pows[j]] * 256./257.;
		tmpvec[j] = multConst(rotvec[j], cnst, bootContext.logp);
	}
	NTL_EXEC_RANGE_END;

	for (long j = 1; j < k1; ++j) {
		addAndEqual(tmpvec[0], tmpvec[j]);
	}
	cipher = tmpvec[0];

	for (long ki = k1; ki < n1; ki += k1) {
		NTL_EXEC_RANGE(k1, first, last);
		for (long j = first; j < last; ++j) {
			complex<double> cnst = ring.ksiM1Pows[ring.M1 - ring.gM1Pows[j + ki]] * 256./257.;
			tmpvec[j] = multConst(rotvec[j], cnst, bootContext.logp);
		}
		NTL_EXEC_RANGE_END;
		for (long j = 1; j < k1; ++j) {
			addAndEqual(tmpvec[0], tmpvec[j]);
		}
		leftRotateFastAndEqual(tmpvec[0], 0, ki);
		addAndEqual(cipher, tmpvec[0]);
	}

	multConstAndEqual(rot, 256./257., bootContext.logp);
	subAndEqual(cipher, rot);
	reScaleByAndEqual(cipher, bootContext.logp);
	delete[] rotvec;
	delete[] tmpvec;
}

void Scheme::coeffToSlotAndEqual(Ciphertext& cipher) {
	coeffToSlotX1AndEqual(cipher);
	coeffToSlotX0AndEqual(cipher);
}

void Scheme::slotToCoeffX0AndEqual(Ciphertext& cipher) {
	long nx = cipher.n0;
	long ny = cipher.n1;
	long lognx = log2(nx);
	long logny = log2(ny);

	long logk = lognx / 2;
	long k = 1 << logk;

	Ciphertext* rotvec = new Ciphertext[k];

	NTL_EXEC_RANGE(k, first, last);
	for (long j = first; j < last; ++j) {
		rotvec[j] = (j == 0) ? cipher : leftRotateFast(cipher, j, 0);
	}
	NTL_EXEC_RANGE_END;

	BootContext bootContext = ring.bootContextMap.at({lognx, logny});

	Ciphertext* tmpvec = new Ciphertext[k];

	NTL_EXEC_RANGE(k, first, last);
	for (long j = first; j < last; ++j) {
		tmpvec[j] = multPolyNTTX0(rotvec[j], bootContext.rpxInvVec[j], bootContext.logp);
	}
	NTL_EXEC_RANGE_END;

	for (long j = 1; j < k; ++j) {
		addAndEqual(tmpvec[0], tmpvec[j]);
	}
	cipher = tmpvec[0];

	for (long ki = k; ki < nx; ki+=k) {
		NTL_EXEC_RANGE(k, first, last);
		for (long j = first; j < last; ++j) {
			tmpvec[j] = multPolyNTTX0(rotvec[j], bootContext.rpxInvVec[j + ki], bootContext.logp);
		}
		NTL_EXEC_RANGE_END;

		for (long j = 1; j < k; ++j) {
			addAndEqual(tmpvec[0], tmpvec[j]);
		}

		leftRotateFastAndEqual(tmpvec[0], ki, 0);
		addAndEqual(cipher, tmpvec[0]);
	}
	reScaleByAndEqual(cipher, bootContext.logp);
	delete[] rotvec;
	delete[] tmpvec;
}

void Scheme::slotToCoeffX1AndEqual(Ciphertext& cipher) {
	long nx = cipher.n0;
	long ny = cipher.n1;
	long lognx = log2(nx);
	long logny = log2(ny);

	long logk = logny / 2;
	long k = 1 << logk;

	Ciphertext* rotvec = new Ciphertext[k];

	NTL_EXEC_RANGE(k, first, last);
	for (long j = first; j < last; ++j) {
		rotvec[j] = (j == 0) ? cipher : leftRotateFast(cipher, 0, j);
	}
	NTL_EXEC_RANGE_END;

	BootContext bootContext = ring.bootContextMap.at({lognx, logny});

	Ciphertext* tmpvec = new Ciphertext[k];

	NTL_EXEC_RANGE(k, first, last);
	for (long j = first; j < last; ++j) {
		complex<double> cnst = ring.ksiM1Pows[ring.gM1Pows[(ring.N1 - j)%ring.N1]];
		tmpvec[j] = multConst(rotvec[j], cnst, bootContext.logp);
	}
	NTL_EXEC_RANGE_END;

	for (long j = 1; j < k; ++j) {
		addAndEqual(tmpvec[0], tmpvec[j]);
	}
	cipher = tmpvec[0];

	for (long ki = k; ki < ny; ki+=k) {
		NTL_EXEC_RANGE(k, first, last);
		for (long j = first; j < last; ++j) {
			complex<double> cnst = ring.ksiM1Pows[ring.gM1Pows[(ring.N1 - j - ki)%ring.N1]];
			tmpvec[j] = multConst(rotvec[j], cnst, bootContext.logp);
		}
		NTL_EXEC_RANGE_END;

		for (long j = 1; j < k; ++j) {
			addAndEqual(tmpvec[0], tmpvec[j]);
		}

		leftRotateFastAndEqual(tmpvec[0], 0, ki);
		addAndEqual(cipher, tmpvec[0]);
	}
	multByMonomialAndEqual(cipher, 0, ring.M1 - 1);

	reScaleByAndEqual(cipher, bootContext.logp);
	delete[] rotvec;
	delete[] tmpvec;
}

void Scheme::slotToCoeffAndEqual(Ciphertext& cipher) {
	slotToCoeffX0AndEqual(cipher);
	slotToCoeffX1AndEqual(cipher);
}

void Scheme::exp2piAndEqual(Ciphertext& cipher, long logp) {
	Ciphertext cipher2 = square(cipher);
	reScaleByAndEqual(cipher2, logp);

	Ciphertext cipher4 = square(cipher2);
	reScaleByAndEqual(cipher4, logp);

	RR c = 1/(2*Pi);
	Ciphertext cipher01 = addConst(cipher, c, logp);

	c = 2*Pi;
	multConstAndEqual(cipher01, c, logp);
	reScaleByAndEqual(cipher01, logp);

	c = 3/(2*Pi);
	Ciphertext cipher23 = addConst(cipher, c, logp);

	c = 4*Pi*Pi*Pi/3;
	multConstAndEqual(cipher23, c, logp);
	reScaleByAndEqual(cipher23, logp);

	multAndEqual(cipher23, cipher2);
	reScaleByAndEqual(cipher23, logp);

	addAndEqual(cipher23, cipher01);

	c = 5/(2*Pi);
	Ciphertext cipher45 = addConst(cipher, c, logp);

	c = 4*Pi*Pi*Pi*Pi*Pi/15;
	multConstAndEqual(cipher45, c, logp);
	reScaleByAndEqual(cipher45, logp);

	c = 7/(2*Pi);
	addConstAndEqual(cipher, c, logp);

	c = 8*Pi*Pi*Pi*Pi*Pi*Pi*Pi/315;
	multConstAndEqual(cipher, c, logp);
	reScaleByAndEqual(cipher, logp);

	multAndEqual(cipher, cipher2);
	reScaleByAndEqual(cipher, logp);

	modDownByAndEqual(cipher45, logp);
	addAndEqual(cipher, cipher45);

	multAndEqual(cipher, cipher4);
	reScaleByAndEqual(cipher, logp);

	modDownByAndEqual(cipher23, logp);
	addAndEqual(cipher, cipher23);
}

void Scheme::evalExpAndEqual(Ciphertext& cipher, long logT, long logI) {
	//TODO change method for different cases
	long logn0 = log2(cipher.n0);
	long logn1 = log2(cipher.n1);

	BootContext bootContext = ring.bootContextMap.at({logn0, logn1});

	Ciphertext tmp = conjugate(cipher);
	Ciphertext cimag = sub(cipher, tmp);
	addAndEqual(cipher, tmp);
	imultAndEqual(cipher);

	divPo2AndEqual(cipher, logT + ring.logN);
	divPo2AndEqual(cimag, logT + ring.logN);

	exp2piAndEqual(cipher, bootContext.logp);
	exp2piAndEqual(cimag, bootContext.logp);

	for (long i = 0; i < logI + logT; ++i) {
		squareAndEqual(cipher);
		squareAndEqual(cimag);
		reScaleByAndEqual(cipher, bootContext.logp);
		reScaleByAndEqual(cimag, bootContext.logp);
	}

	tmp = conjugate(cimag);
	subAndEqual(cimag, tmp);
	tmp = conjugate(cipher);
	subAndEqual(cipher, tmp);
	imultAndEqual(cipher);
	subAndEqual2(cimag, cipher);

	RR c = 0.25/Pi;
	multConstAndEqual(cipher, c, bootContext.logp);
	reScaleByAndEqual(cipher, bootContext.logp + logI);
}

void Scheme::bootstrapXAndEqual(Ciphertext& cipher, long logq, long logQ, long logT, long logI) {
	long lognx = log2(cipher.n0);

	long logp = cipher.logp;

	modDownToAndEqual(cipher, logq);
	normalizeAndEqual(cipher);

	cipher.logq = logQ;
	cipher.logp = logq + logI;

	for (long i = lognx; i < ring.logN0h; ++i) {
		Ciphertext rot = leftRotateFast(cipher, (1 << i), 0);
		addAndEqual(cipher, rot);
	}
	coeffToSlotX0AndEqual(cipher);
	divPo2AndEqual(cipher, ring.logN0h);
	evalExpAndEqual(cipher, logT, logI);
	slotToCoeffX0AndEqual(cipher);
	cipher.logp = logp;
}

void Scheme::bootstrapYAndEqual(Ciphertext& cipher, long logq, long logQ, long logT, long logI) {
	long logny = log2(cipher.n1);

	long logp = cipher.logp;

	modDownToAndEqual(cipher, logq);
	normalizeAndEqual(cipher);

	cipher.logq = logQ;
	cipher.logp = logq + logI;

	for (long i = logny; i < ring.logN1; ++i) {
		Ciphertext rot = leftRotateFast(cipher, 0, (1 << i));
		addAndEqual(cipher, rot);
	}
	coeffToSlotX1AndEqual(cipher);
	divPo2AndEqual(cipher, ring.logN1);
	evalExpAndEqual(cipher, logT, logI);
	slotToCoeffX1AndEqual(cipher);
	cipher.logp = logp;
}

void Scheme::bootstrapAndEqual(Ciphertext& cipher, long logq, long logQ, long logT, long logI) {
	long lognx = log2(cipher.n0);
	long logny = log2(cipher.n1);

	long logp = cipher.logp;

	modDownToAndEqual(cipher, logq);
	normalizeAndEqual(cipher);

	cipher.logq = logQ;
	cipher.logp = logq + logI;

	for (long i = lognx; i < ring.logN0h; ++i) {
		Ciphertext rot = leftRotateFast(cipher, (1 << i), 0);
		addAndEqual(cipher, rot);
	}
	for (long i = logny; i < ring.logN1; ++i) {
		Ciphertext rot = leftRotateFast(cipher, 0, (1 << i));
		addAndEqual(cipher, rot);
	}
	coeffToSlotAndEqual(cipher);
	divPo2AndEqual(cipher, ring.logN0h + ring.logN1);
	evalExpAndEqual(cipher, logT, logI);
	slotToCoeffAndEqual(cipher);
	cipher.logp = logp;
}
