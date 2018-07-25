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

Scheme::Scheme(SecretKey& secretKey, Ring2XY& ring) : ring(ring) {
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

	ring.mult(bx, secretKey.sxy, ax, ring.QQ);

	ring.sub(bx, ex, bx, ring.QQ);

	delete[] ex;

	keyMap.insert(pair<long, Key>(ENCRYPTION, Key(ax, bx)));
}

void Scheme::addMultKey(SecretKey& secretKey) {
	ZZ* ex = new ZZ[ring.N];
	ZZ* ax = new ZZ[ring.N];
	ZZ* bx = new ZZ[ring.N];
	ZZ* sxsx = new ZZ[ring.N];

	ring.mult(sxsx, secretKey.sxy, secretKey.sxy, ring.Q);
	ring.leftShiftAndEqual(sxsx, ring.logQ, ring.QQ);
	ring.sampleUniform(ax, ring.logQQ);
	ring.sampleGauss(ex);
	ring.addAndEqual(ex, sxsx, ring.QQ);
	ring.mult(bx, secretKey.sxy, ax, ring.QQ);
	ring.sub(bx, ex, bx, ring.QQ);

	delete[] ex;
	delete[] sxsx;

	keyMap.insert(pair<long, Key>(MULTIPLICATION, Key(ax, bx)));
}

void Scheme::addLeftRotKey(SecretKey& secretKey, long rx, long ry) {
	ZZ* exy = new ZZ[ring.N];
	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];
	ZZ* rysxy = new ZZ[ring.N];
	ring.leftRotate(rysxy, secretKey.sxy, rx, ry);
	ring.leftShiftAndEqual(rysxy, ring.logQ, ring.QQ);
	ring.sampleUniform(axy, ring.logQQ);
	ring.sampleGauss(exy);
	ring.addAndEqual(exy, rysxy, ring.QQ);
	ring.mult(bxy, secretKey.sxy, axy, ring.QQ);
	ring.sub(bxy, exy, bxy, ring.QQ);

	delete[] exy;
	delete[] rysxy;

	leftRotKeyMap.insert(pair<pair<long, long>, Key>({rx, ry}, Key(axy, bxy)));
}

void Scheme::addLeftXRotKeys(SecretKey& secretKey) {
	for (long ix = 1; ix < ring.Nxh; ix <<= 1) {
		if(leftRotKeyMap.find({ix, 0}) == leftRotKeyMap.end()) {
			addLeftRotKey(secretKey, ix, 0);
		}
	}
}

void Scheme::addLeftYRotKeys(SecretKey& secretKey) {
	for (long iy = 1; iy < ring.Ny; iy <<=1) {
		if(leftRotKeyMap.find({0, iy}) == leftRotKeyMap.end()) {
			addLeftRotKey(secretKey, 0, iy);
		}
	}
}

void Scheme::addRightXRotKeys(SecretKey& secretKey) {
	for (long ix = 1; ix < ring.Nxh; ix <<=1) {
		long idx = ring.Nxh - ix;
		if(leftRotKeyMap.find({idx, 0}) == leftRotKeyMap.end()) {
			addLeftRotKey(secretKey, idx, 0);
		}
	}
}

void Scheme::addRightYRotKeys(SecretKey& secretKey) {
	for (long iy = 1; iy < ring.Ny; iy<<=1) {
		long idx = ring.Ny - iy;
		if(leftRotKeyMap.find({0, idx}) == leftRotKeyMap.end()) {
			addLeftRotKey(secretKey, 0, idx);
		}
	}
}

void Scheme::addConjKey(SecretKey& secretKey) {
	ZZ* exy = new ZZ[ring.N];
	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];
	ZZ* conjsxy = new ZZ[ring.N];

	ring.conjugate(conjsxy, secretKey.sxy);
	ring.leftShiftAndEqual(conjsxy, ring.logQ, ring.QQ);
	ring.sampleUniform(axy, ring.logQQ);
	ring.sampleGauss(exy);
	ring.addAndEqual(exy, conjsxy, ring.QQ);
	ring.mult(bxy, secretKey.sxy, axy, ring.QQ);
	ring.sub(bxy, exy, bxy, ring.QQ);

	delete[] exy;
	delete[] conjsxy;

	keyMap.insert(pair<long, Key>(CONJUGATION, Key(axy, bxy)));
}

void Scheme::addTranspKey(SecretKey& secretKey) {
	ZZ* exy = new ZZ[ring.N];
	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];
	ZZ* trsxy = new ZZ[ring.N];

	ring.transpose(trsxy, secretKey.sxy);
	ring.leftShiftAndEqual(trsxy, ring.logQ, ring.QQ);
	ring.sampleUniform(axy, ring.logQQ);
	ring.sampleGauss(exy);
	ring.addAndEqual(exy, trsxy, ring.QQ);
	ring.mult(bxy, secretKey.sxy, axy, ring.QQ);
	ring.sub(bxy, exy, bxy, ring.QQ);

	delete[] exy;
	delete[] trsxy;

	keyMap.insert(pair<long, Key>(TRANSPOSITION, Key(axy, bxy)));
}

void Scheme::addSquareMatrixKeys(SecretKey& secretKey, long lognx) {
	ring.addMatrixContext(lognx);
	long nx = 1 << lognx;
	for (long ix = 1; ix < nx; ++ix) {
		long idx = ring.Nxh - ix;
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
//	ring.encode(mx, vals, nx, logp + ring.logQ);
	ring.encode(mx, vals, nx, logp);
	return Plaintext(mx, logp, logq, ring.Nx, ring.Ny, nx, ny, true);
}

Plaintext Scheme::encode(double* vals, long nx, long ny, long logp, long logq) {
	ZZ* mxy = new ZZ[ring.N];
	ring.encode(mxy, vals, nx, logp + ring.logQ);
	return Plaintext(mxy, logp, logq, ring.Nx, ring.Ny, nx, ny, true);
}

Plaintext Scheme::encodeSingle(complex<double> val, long logp, long logq) {
	ZZ* mxy = new ZZ[ring.N];

	mxy[0] = EvaluatorUtils::scaleUpToZZ(val.real(), logp + ring.logQ);
	mxy[ring.Nh] = EvaluatorUtils::scaleUpToZZ(val.imag(), logp + ring.logQ);

	return Plaintext(mxy, logp, logq, ring.Nx, ring.Ny, 1, 1, true);
}

Plaintext Scheme::encodeSingle(double val, long logp, long logq) {
	ZZ* mx = new ZZ[ring.N];

	mx[0] = EvaluatorUtils::scaleUpToZZ(val, logp + ring.logQ);

	return Plaintext(mx, logp, logq, ring.Nx, ring.Ny, 1, 1, false);
}

complex<double>* Scheme::decode(Plaintext& msg) {
	complex<double>* vals = new complex<double>[msg.nx * ring.Ny];
	ring.decode(msg.mxy, vals, msg.nx, msg.logp, msg.logq);
	return vals;
}

complex<double> Scheme::decodeSingle(Plaintext& msg) {
	complex<double> res;
	ZZ tmp;
	ZZ q = ring.qvec[msg.logq];
	ZZ qh = ring.qvec[msg.logq - 1];
	AddMod(tmp, msg.mxy[0], -msg.mxy[ring.Nxh + ring.Nh], q);
	while(tmp < 0) tmp += q;
	while(tmp > qh) tmp -= q;
	res.real(EvaluatorUtils::scaleDownToReal(tmp, msg.logp));

	if(msg.isComplex) {
		AddMod(tmp, msg.mxy[ring.Nh], msg.mxy[ring.Nxh], q);
		while(tmp < 0) tmp += q;
		while(tmp > qh) tmp -= q;
		res.imag(EvaluatorUtils::scaleDownToReal(tmp, msg.logp));
	}
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

	ring.mult(axy, vxy, key.axy, qQ);
	ring.sampleGauss(exy);
	ring.addAndEqual(axy, exy, qQ);

	ring.mult(bxy, vxy, key.bxy, qQ);
	ring.sampleGauss(exy);
	ring.addAndEqual(bxy, exy, qQ);

	ring.addAndEqual(bxy, msg.mxy, qQ);
	ring.rightShiftAndEqual(axy, ring.logQ);
	ring.rightShiftAndEqual(bxy, ring.logQ);

	delete[] exy;
	delete[] vxy;

	return Ciphertext(axy, bxy, msg.logp, msg.logq, msg.Nx, msg.Ny, msg.nx, msg.ny, msg.isComplex);
}

Ciphertext Scheme::encrypt(complex<double>* vals, long nx, long ny, long logp, long logq) {
	Plaintext msg = encode(vals, nx, ny, logp, logq);
	return encryptMsg(msg);
}

Ciphertext Scheme::encrypt(double* vals, long nx, long ny, long logp, long logq) {
	Plaintext msg = encode(vals, nx, ny, logp, logq);
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
	czeros.isComplex = true;
	czeros.nx = nx;
	czeros.ny = ny;
	return czeros;
}

Plaintext Scheme::decryptMsg(SecretKey& secretKey, Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ* mxy = new ZZ[ring.N];
	ring.mult(mxy, cipher.axy, secretKey.sxy, q);
	ring.addAndEqual(mxy, cipher.bxy, q);
	return Plaintext(mxy, cipher.logp, cipher.logq, cipher.Nx, cipher.Ny, cipher.nx, cipher.ny, cipher.isComplex);
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
	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.negate(axy, cipher.axy);
	ring.negate(bxy, cipher.bxy);

	return Ciphertext(axy, bxy, cipher.logp, cipher.logq, cipher.Nx, cipher.Ny, cipher.nx, cipher.ny, cipher.isComplex);
}

void Scheme::negateAndEqual(Ciphertext& cipher) {
	ring.negateAndEqual(cipher.axy);
	ring.negateAndEqual(cipher.bxy);
}

Ciphertext Scheme::add(Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZ q = ring.qvec[cipher1.logq];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.add(axy, cipher1.axy, cipher2.axy, q);
	ring.add(bxy, cipher1.bxy, cipher2.bxy, q);

	return Ciphertext(axy, bxy, cipher1.logp, cipher1.logq, cipher1.Nx, cipher1.Ny, cipher1.nx, cipher1.ny, cipher1.isComplex);
}

void Scheme::addAndEqual(Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZ q = ring.qvec[cipher1.logq];

	ring.addAndEqual(cipher1.axy, cipher2.axy, q);
	ring.addAndEqual(cipher1.bxy, cipher2.bxy, q);
}

Ciphertext Scheme::addConst(Ciphertext& cipher, double cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	Ciphertext res = cipher;
	ZZ cnstZZ = logp < 0 ? EvaluatorUtils::scaleUpToZZ(cnst, cipher.logp) : EvaluatorUtils::scaleUpToZZ(cnst, logp);
	AddMod(res.bxy[0], res.bxy[0], cnstZZ, q);
	return res;
}

Ciphertext Scheme::addConst(Ciphertext& cipher, RR& cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	Ciphertext res = cipher;
	ZZ cnstZZ = logp < 0 ? EvaluatorUtils::scaleUpToZZ(cnst, cipher.logp) : EvaluatorUtils::scaleUpToZZ(cnst, logp);
	AddMod(res.bxy[0], res.bxy[0], cnstZZ, q);
	return res;
}

void Scheme::addConstAndEqual(Ciphertext& cipher, double cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ cnstZZ = logp < 0 ? EvaluatorUtils::scaleUpToZZ(cnst, cipher.logp) : EvaluatorUtils::scaleUpToZZ(cnst, logp);

	AddMod(cipher.bxy[0], cipher.bxy[0], cnstZZ, q);
}

void Scheme::addConstAndEqual(Ciphertext& cipher, RR& cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ cnstZZ = logp < 0 ? EvaluatorUtils::scaleUpToZZ(cnst, cipher.logp) : EvaluatorUtils::scaleUpToZZ(cnst, logp);

	AddMod(cipher.bxy[0], cipher.bxy[0], cnstZZ, q);
}

Ciphertext Scheme::addPoly(Ciphertext& cipher, ZZ* poly, long logp) {
	ZZ q = ring.qvec[cipher.logq];
	Ciphertext res = cipher;
	ring.addAndEqual(res.bxy, poly, q);
	return res;
}

void Scheme::addPolyAndEqual(Ciphertext& cipher, ZZ* poly, long logp) {
	ZZ q = ring.qvec[cipher.logq];
	ring.addAndEqual(cipher.bxy, poly, q);
}

Ciphertext Scheme::sub(Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZ q = ring.qvec[cipher1.logq];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.sub(axy, cipher1.axy, cipher2.axy, q);
	ring.sub(bxy, cipher1.bxy, cipher2.bxy, q);

	return Ciphertext(axy, bxy, cipher1.logp, cipher1.logq, cipher1.Nx, cipher1.Ny, cipher1.nx, cipher1.ny, cipher1.isComplex);
}

void Scheme::subAndEqual(Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZ q = ring.qvec[cipher1.logq];

	ring.subAndEqual(cipher1.axy, cipher2.axy, q);
	ring.subAndEqual(cipher1.bxy, cipher2.bxy, q);
}

void Scheme::subAndEqual2(Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZ q = ring.qvec[cipher1.logq];

	ring.subAndEqual2(cipher1.axy, cipher2.axy, q);
	ring.subAndEqual2(cipher1.bxy, cipher2.bxy, q);
}

Ciphertext Scheme::imult(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.multByMonomial(axy, cipher.axy, ring.Nxh, 0);
	ring.multByMonomial(bxy, cipher.bxy, ring.Nxh, 0);

	return Ciphertext(axy, bxy, cipher.logp, cipher.logq, cipher.Nx, cipher.Ny, cipher.nx, cipher.ny, cipher.isComplex);
}

Ciphertext Scheme::idiv(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.multByMonomial(axy, cipher.axy, 3 * ring.Nxh, 0);
	ring.multByMonomial(bxy, cipher.bxy, 3 * ring.Nxh, 0);

	return Ciphertext(axy, bxy, cipher.logp, cipher.logq, cipher.Nx, cipher.Ny, cipher.nx, cipher.ny, cipher.isComplex);
}

void Scheme::imultAndEqual(Ciphertext& cipher) {
	ring.multByMonomialAndEqual(cipher.axy, ring.Nxh, 0);
	ring.multByMonomialAndEqual(cipher.bxy, ring.Nxh, 0);
}

void Scheme::idivAndEqual(Ciphertext& cipher) {
	ring.multByMonomialAndEqual(cipher.axy, 3 * ring.Nxh, 0);
	ring.multByMonomialAndEqual(cipher.bxy, 3 * ring.Nxh, 0);
}
Ciphertext Scheme::mult(Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZ q = ring.qvec[cipher1.logq];
	ZZ qQ = ring.qvec[cipher1.logq + ring.logQ];

	ZZ* ab1xy = new ZZ[ring.N];
	ZZ* ab2xy = new ZZ[ring.N];
	ZZ* aaxy = new ZZ[ring.N];
	ZZ* bbxy = new ZZ[ring.N];
	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.add(ab1xy, cipher1.axy, cipher1.bxy, q);
	ring.add(ab2xy, cipher2.axy, cipher2.bxy, q);
	ring.multAndEqual(ab1xy, ab2xy, q);

	ring.mult(bbxy, cipher1.bxy, cipher2.bxy, q);
	ring.mult(aaxy, cipher1.axy, cipher2.axy, q);

	Key key = keyMap.at(MULTIPLICATION);

	ring.mult(axy, aaxy, key.axy, qQ);
	ring.mult(bxy, aaxy, key.bxy, qQ);

	ring.rightShiftAndEqual(axy, ring.logQ);
	ring.rightShiftAndEqual(bxy, ring.logQ);

	ring.addAndEqual(axy, ab1xy, q);
	ring.subAndEqual(axy, bbxy, q);
	ring.subAndEqual(axy, aaxy, q);
	ring.addAndEqual(bxy, bbxy, q);

	delete[] ab1xy;
	delete[] ab2xy;
	delete[] aaxy;
	delete[] bbxy;

	return Ciphertext(axy, bxy, cipher1.logp + cipher2.logp, cipher1.logq, cipher1.Nx, cipher1.Ny, cipher1.nx, cipher1.ny, cipher1.isComplex);
}

void Scheme::multAndEqual(Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZ q = ring.qvec[cipher1.logq];
	ZZ qQ = ring.qvec[cipher1.logq + ring.logQ];

	ZZ* ab1xy = new ZZ[ring.N];
	ZZ* ab2xy = new ZZ[ring.N];
	ZZ* aaxy = new ZZ[ring.N];
	ZZ* bbxy = new ZZ[ring.N];

	ring.add(ab1xy, cipher1.axy, cipher1.bxy, q);
	ring.add(ab2xy, cipher2.axy, cipher2.bxy, q);
	ring.multAndEqual(ab1xy, ab2xy, q);

	ring.mult(aaxy, cipher1.axy, cipher2.axy, q);
	ring.mult(bbxy, cipher1.bxy, cipher2.bxy, q);

	Key key = keyMap.at(MULTIPLICATION);

	ring.mult(cipher1.axy, aaxy, key.axy, qQ);
	ring.mult(cipher1.bxy, aaxy, key.bxy, qQ);

	ring.rightShiftAndEqual(cipher1.axy, ring.logQ);
	ring.rightShiftAndEqual(cipher1.bxy, ring.logQ);

	ring.addAndEqual(cipher1.axy, ab1xy, q);
	ring.subAndEqual(cipher1.axy, bbxy, q);
	ring.subAndEqual(cipher1.axy, aaxy, q);
	ring.addAndEqual(cipher1.bxy, bbxy, q);

	cipher1.logp += cipher2.logp;
	delete[] ab1xy;
	delete[] ab2xy;
	delete[] aaxy;
	delete[] bbxy;
}

Ciphertext Scheme::square(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ qQ = ring.qvec[cipher.logq + ring.logQ];

	ZZ* abxy = new ZZ[ring.N];
	ZZ* aaxy = new ZZ[ring.N];
	ZZ* bbxy = new ZZ[ring.N];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.square(bbxy, cipher.bxy, q);
	ring.mult(abxy, cipher.axy, cipher.bxy, q);
	ring.addAndEqual(abxy, abxy, q);
	ring.square(aaxy, cipher.axy, q);

	Key key = keyMap.at(MULTIPLICATION);

	ring.mult(axy, aaxy, key.axy, qQ);
	ring.mult(bxy, aaxy, key.bxy, qQ);

	ring.rightShiftAndEqual(axy, ring.logQ);
	ring.rightShiftAndEqual(bxy, ring.logQ);

	ring.addAndEqual(axy, abxy, q);
	ring.addAndEqual(bxy, bbxy, q);

	delete[] abxy;
	delete[] aaxy;
	delete[] bbxy;


	return Ciphertext(axy, bxy, 2 * cipher.logp, cipher.logq, cipher.Nx, cipher.Ny, cipher.nx, cipher.ny, cipher.isComplex);
}

void Scheme::squareAndEqual(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ qQ = ring.qvec[cipher.logq + ring.logQ];

	ZZ* abxy = new ZZ[ring.N];
	ZZ* aaxy = new ZZ[ring.N];
	ZZ* bbxy = new ZZ[ring.N];

	ring.square(bbxy, cipher.bxy, q);
	ring.mult(abxy, cipher.bxy, cipher.axy, q);
	ring.addAndEqual(abxy, abxy, q);
	ring.square(aaxy, cipher.axy, q);

	Key key = keyMap.at(MULTIPLICATION);
	ring.mult(cipher.axy, aaxy, key.axy, qQ);
	ring.mult(cipher.bxy, aaxy, key.bxy, qQ);

	ring.rightShiftAndEqual(cipher.axy, ring.logQ);
	ring.rightShiftAndEqual(cipher.bxy, ring.logQ);

	ring.addAndEqual(cipher.axy, abxy, q);
	ring.addAndEqual(cipher.bxy, bbxy, q);

	cipher.logp *= 2;

	delete[] abxy;
	delete[] aaxy;
	delete[] bbxy;
}

Ciphertext Scheme::multByConst(Ciphertext& cipher, RR& cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ZZ cnstZZ = EvaluatorUtils::scaleUpToZZ(cnst, logp);

	ring.multByConst(axy, cipher.axy, cnstZZ, q);
	ring.multByConst(bxy, cipher.bxy, cnstZZ, q);

	return Ciphertext(axy, bxy, cipher.logp + logp, cipher.logq, cipher.Nx, cipher.Ny, cipher.nx, cipher.ny, cipher.isComplex);
}


Ciphertext Scheme::multByConst(Ciphertext& cipher, double cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ZZ cnstZZ = EvaluatorUtils::scaleUpToZZ(cnst, logp);

	ring.multByConst(axy, cipher.axy, cnstZZ, q);
	ring.multByConst(bxy, cipher.bxy, cnstZZ, q);

	return Ciphertext(axy, bxy, cipher.logp + logp, cipher.logq, cipher.Nx, cipher.Ny, cipher.nx, cipher.ny, cipher.isComplex);
}

void Scheme::multByConstAndEqual(Ciphertext& cipher, RR& cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ cnstZZ = EvaluatorUtils::scaleUpToZZ(cnst, logp);

	ring.multByConstAndEqual(cipher.axy, cnstZZ, q);
	ring.multByConstAndEqual(cipher.bxy, cnstZZ, q);

	cipher.logp += logp;
}

void Scheme::multByConstAndEqual(Ciphertext& cipher, double cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ cnstZZ = EvaluatorUtils::scaleUpToZZ(cnst, logp);

	ring.multByConstAndEqual(cipher.axy, cnstZZ, q);
	ring.multByConstAndEqual(cipher.bxy, cnstZZ, q);

	cipher.logp += logp;
}

Ciphertext Scheme::multByXPoly(Ciphertext& cipher, ZZ* xpoly, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.multXpoly(axy, cipher.axy, xpoly, q);
	ring.multXpoly(bxy, cipher.bxy, xpoly, q);

	return Ciphertext(axy, bxy, cipher.logp + logp, cipher.logq, cipher.Nx, cipher.Ny, cipher.nx, cipher.ny, cipher.isComplex);
}

void Scheme::multByXPolyAndEqual(Ciphertext& cipher, ZZ* xpoly, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ring.multXpolyAndEqual(cipher.axy, xpoly, q);
	ring.multXpolyAndEqual(cipher.bxy, xpoly, q);

	cipher.logp += logp;
}

Ciphertext Scheme::multByYPoly(Ciphertext& cipher, ZZ* ypoly, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.multYpoly(axy, cipher.axy, ypoly, q);
	ring.multYpoly(bxy, cipher.bxy, ypoly, q);

	return Ciphertext(axy, bxy, cipher.logp + logp, cipher.logq, cipher.Nx, cipher.Ny, cipher.nx, cipher.ny, cipher.isComplex);
}

void Scheme::multByYPolyAndEqual(Ciphertext& cipher, ZZ* ypoly, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ring.multYpolyAndEqual(cipher.axy, ypoly, q);
	ring.multYpolyAndEqual(cipher.bxy, ypoly, q);

	cipher.logp += logp;
}

Ciphertext Scheme::multByPoly(Ciphertext& cipher, ZZ* poly, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.mult(axy, cipher.axy, poly, q);
	ring.mult(bxy, cipher.bxy, poly, q);

	return Ciphertext(axy, bxy, cipher.logp + logp, cipher.logq, cipher.Nx, cipher.Ny, cipher.nx, cipher.ny, cipher.isComplex);
}

void Scheme::multByPolyAndEqual(Ciphertext& cipher, ZZ* poly, long logp) {
	ZZ q = ring.qvec[cipher.logq];

	ring.multAndEqual(cipher.axy, poly, q);
	ring.multAndEqual(cipher.bxy, poly, q);

	cipher.logp += logp;
}

Ciphertext Scheme::multByMonomial(Ciphertext& cipher, const long dx, const long dy) {
	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.multByMonomial(axy, cipher.axy, dx, dy);
	ring.multByMonomial(bxy, cipher.bxy, dx, dy);

	return Ciphertext(axy, bxy, cipher.logp, cipher.logq, cipher.Nx, cipher.Ny, cipher.nx, cipher.ny, cipher.isComplex);
}

void Scheme::multByMonomialAndEqual(Ciphertext& cipher, const long dx, const long dy) {
	ring.multByMonomialAndEqual(cipher.axy, dx, dy);
	ring.multByMonomialAndEqual(cipher.bxy, dx, dy);
}

Ciphertext Scheme::multByPo2(Ciphertext& cipher, long degree) {
	ZZ q = ring.qvec[cipher.logq];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.leftShift(axy, cipher.axy, degree, q);
	ring.leftShift(bxy, cipher.bxy, degree, q);

	return Ciphertext(axy, bxy, cipher.logp, cipher.logq, cipher.Nx, cipher.Ny, cipher.nx, cipher.ny, cipher.isComplex);
}

void Scheme::multByPo2AndEqual(Ciphertext& cipher, long degree) {
	ZZ q = ring.qvec[cipher.logq];

	ring.leftShiftAndEqual(cipher.axy, degree, q);
	ring.leftShiftAndEqual(cipher.bxy, degree, q);
}

void Scheme::doubleAndEqual(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];

	ring.doubleAndEqual(cipher.axy, q);
	ring.doubleAndEqual(cipher.bxy, q);
}

Ciphertext Scheme::divByPo2(Ciphertext& cipher, long logd) {
	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.rightShift(axy, cipher.axy, logd);
	ring.rightShift(bxy, cipher.bxy, logd);

	return Ciphertext(axy, bxy, cipher.logp, cipher.logq - logd, cipher.Nx, cipher.Ny, cipher.nx, cipher.ny, cipher.isComplex);
}

void Scheme::divByPo2AndEqual(Ciphertext& cipher, long logd) {
	ring.rightShiftAndEqual(cipher.axy, logd);
	ring.rightShiftAndEqual(cipher.bxy, logd);
	cipher.logq -= logd;
}


//----------------------------------------------------------------------------------
//   RESCALING & MODULUS DOWN
//----------------------------------------------------------------------------------


Ciphertext Scheme::reScaleBy(Ciphertext& cipher, long dlogq) {
	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.rightShift(axy, cipher.axy, dlogq);
	ring.rightShift(bxy, cipher.bxy, dlogq);

	return Ciphertext(axy, bxy, cipher.logp - dlogq, cipher.logq - dlogq, cipher.Nx, cipher.Ny, cipher.nx, cipher.ny, cipher.isComplex);
}

Ciphertext Scheme::reScaleTo(Ciphertext& cipher, long logq) {
	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	long dlogq = cipher.logq - logq;
	ring.rightShift(axy, cipher.axy, dlogq);
	ring.rightShift(bxy, cipher.bxy, dlogq);

	return Ciphertext(axy, bxy, cipher.logp - dlogq, logq, cipher.Nx, cipher.Ny, cipher.nx, cipher.ny, cipher.isComplex);
}

void Scheme::reScaleByAndEqual(Ciphertext& cipher, long dlogq) {
	ring.rightShiftAndEqual(cipher.axy, dlogq);
	ring.rightShiftAndEqual(cipher.bxy, dlogq);
	cipher.logq -= dlogq;
	cipher.logp -= dlogq;
}

void Scheme::reScaleToAndEqual(Ciphertext& cipher, long logq) {
	long dlogq = cipher.logq - logq;
	ring.rightShiftAndEqual(cipher.axy, dlogq);
	ring.rightShiftAndEqual(cipher.bxy, dlogq);
	cipher.logq = logq;
	cipher.logp -= dlogq;
}

Ciphertext Scheme::modDownBy(Ciphertext& cipher, long dlogq) {
	long logq = cipher.logq - dlogq;
	ZZ q = ring.qvec[logq];
	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.mod(axy, cipher.axy, q);
	ring.mod(bxy, cipher.bxy, q);
	return Ciphertext(axy, bxy, cipher.logp, logq, cipher.Nx, cipher.Ny, cipher.nx, cipher.ny, cipher.isComplex);
}

void Scheme::modDownByAndEqual(Ciphertext& cipher, long dlogq) {
	ZZ q = ring.qvec[cipher.logq - dlogq];

	ring.modAndEqual(cipher.axy, q);
	ring.modAndEqual(cipher.bxy, q);

	cipher.logq -= dlogq;
}

Ciphertext Scheme::modDownTo(Ciphertext& cipher, long logq) {
	ZZ q = ring.qvec[logq];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ring.mod(axy, cipher.axy, q);
	ring.mod(bxy, cipher.bxy, q);

	return Ciphertext(axy, bxy, cipher.logp, logq, cipher.Nx, cipher.Ny, cipher.nx, cipher.ny, cipher.isComplex);
}

void Scheme::modDownToAndEqual(Ciphertext& cipher, long newlogq) {
	ZZ q = ring.qvec[newlogq];

	ring.modAndEqual(cipher.axy, q);
	ring.modAndEqual(cipher.bxy, q);

	cipher.logq = newlogq;
}


//----------------------------------------------------------------------------------
//   ROTATIONS & CONJUGATIONS & TRANSPOSITION
//----------------------------------------------------------------------------------


Ciphertext Scheme::leftRotateFast(Ciphertext& cipher, long rx, long ry) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ qQ = ring.qvec[cipher.logq + ring.logQ];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];
	ZZ* rxbxy = new ZZ[ring.N];

	ring.leftRotate(rxbxy, cipher.bxy, rx, ry);
	ring.leftRotate(bxy, cipher.axy, rx, ry);

	Key key = leftRotKeyMap.at({rx, ry});

	ring.mult(axy, bxy, key.axy, qQ);
	ring.multAndEqual(bxy, key.bxy, qQ);

	ring.rightShiftAndEqual(axy, ring.logQ);
	ring.rightShiftAndEqual(bxy, ring.logQ);

	ring.addAndEqual(bxy, rxbxy, q);

	delete[] rxbxy;
	return Ciphertext(axy, bxy, cipher.logp, cipher.logq, cipher.Nx, cipher.Ny, cipher.nx, cipher.ny, cipher.isComplex);
}

Ciphertext Scheme::rightRotateFast(Ciphertext& cipher, long rx, long ry) {
	Ciphertext res = cipher;
	rightRotateFastAndEqual(res, rx, ry);
	return res;
}

void Scheme::leftRotateFastAndEqual(Ciphertext& cipher, long rx, long ry) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ qQ = ring.qvec[cipher.logq + ring.logQ];

	ZZ* rxbxy = new ZZ[ring.N];

	ring.leftRotate(rxbxy, cipher.bxy, rx, ry);
	ring.leftRotate(cipher.bxy, cipher.axy, rx, ry);

	Key key = leftRotKeyMap.at({rx, ry});

	ring.mult(cipher.axy, cipher.bxy, key.axy, qQ);
	ring.multAndEqual(cipher.bxy, key.bxy, qQ);

	ring.rightShiftAndEqual(cipher.axy, ring.logQ);
	ring.rightShiftAndEqual(cipher.bxy, ring.logQ);

	ring.addAndEqual(cipher.bxy, rxbxy, q);

	delete[] rxbxy;
}

void Scheme::rightRotateFastAndEqual(Ciphertext& cipher, long rx, long ry) {
	long rrx = rx == 0 ? 0 : ring.Nxh - rx;
	long rry = ry == 0 ? 0 : ring.Ny - ry;
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


void Scheme::leftRotateAndEqual(Ciphertext& cipher, long rx, long ry) {
	long rxRem = rx % cipher.nx;
	long logrxRem = log2((double)rxRem) + 1;
	long ipow;
	for (long i = 0; i < logrxRem; ++i) {
		if(bit(rxRem, i)) {
			ipow = 1 << i;
			leftRotateFastAndEqual(cipher, ipow, 0);
		}
	}
	long ryRem = ry % cipher.ny;
	long logryRem = log2((double)ryRem) + 1;
	for (long i = 0; i < logryRem; ++i) {
		if(bit(ryRem, i)) {
			ipow = 1 << i;
			leftRotateFastAndEqual(cipher, 0, ipow);
		}
	}
}

void Scheme::rightRotateAndEqual(Ciphertext& cipher, long rx, long ry) {
	long rxRem = rx % cipher.nx;
	long logrxRem = log2((double)rxRem) + 1;
	long ipow;
	for (long i = 0; i < logrxRem; ++i) {
		if(bit(rxRem, i)) {
			ipow = 1 << i;
			rightRotateFastAndEqual(cipher, ipow, 0);
		}
	}
	long ryRem = ry % cipher.ny;
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

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];

	ZZ* conjbxy = new ZZ[ring.N];

	ring.conjugate(conjbxy, cipher.bxy);
	ring.conjugate(bxy, cipher.axy);

	Key key = keyMap.at(CONJUGATION);

	ring.mult(axy, bxy, key.axy, qQ);
	ring.multAndEqual(bxy, key.bxy, qQ);

	ring.rightShiftAndEqual(axy, ring.logQ);
	ring.rightShiftAndEqual(bxy, ring.logQ);

	ring.addAndEqual(bxy, conjbxy, q);

	delete[] conjbxy;

	return Ciphertext(axy, bxy, cipher.logp, cipher.logq, cipher.Nx, cipher.Ny, cipher.nx, cipher.ny, cipher.isComplex);
}

void Scheme::conjugateAndEqual(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ qQ = ring.qvec[cipher.logq + ring.logQ];

	ZZ* conjbxy = new ZZ[ring.N];

	ring.conjugate(conjbxy, cipher.bxy);
	ring.conjugate(cipher.bxy, cipher.axy);

	Key key = keyMap.at(CONJUGATION);

	ring.mult(cipher.axy, cipher.bxy, key.axy, qQ);
	ring.multAndEqual(cipher.bxy, key.bxy, qQ);

	ring.rightShiftAndEqual(cipher.axy, ring.logQ);
	ring.rightShiftAndEqual(cipher.bxy, ring.logQ);

	ring.addAndEqual(cipher.bxy, conjbxy, q);

	delete[] conjbxy;
}

Ciphertext Scheme::transpose(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ qQ = ring.qvec[cipher.logq + ring.logQ];

	ZZ* axy = new ZZ[ring.N];
	ZZ* bxy = new ZZ[ring.N];
	ZZ* trbxy = new ZZ[ring.N];

	ring.transpose(trbxy, cipher.bxy);
	ring.transpose(bxy, cipher.axy);

	Key key = keyMap.at(TRANSPOSITION);

	ring.mult(axy, bxy, key.axy, qQ);
	ring.multAndEqual(bxy, key.bxy, qQ);

	ring.rightShiftAndEqual(axy, ring.logQ);
	ring.rightShiftAndEqual(bxy, ring.logQ);

	ring.addAndEqual(bxy, trbxy, q);

	delete[] trbxy;

	return Ciphertext(axy, bxy, cipher.logp, cipher.logq, cipher.Nx, cipher.Ny, cipher.ny, cipher.nx, cipher.isComplex);
}

void Scheme::transposeAndEqual(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ qQ = ring.qvec[cipher.logq + ring.logQ];

	ZZ* trbxy = new ZZ[ring.N];

	ring.transpose(trbxy, cipher.bxy);
	ring.transpose(cipher.bxy, cipher.axy);

	Key key = keyMap.at(TRANSPOSITION);

	ring.mult(cipher.axy, cipher.bxy, key.axy, qQ);
	ring.multAndEqual(cipher.bxy, key.bxy, qQ);

	ring.rightShiftAndEqual(cipher.axy, ring.logQ);
	ring.rightShiftAndEqual(cipher.bxy, ring.logQ);

	ring.addAndEqual(cipher.bxy, trbxy, q);

	swap(cipher.nx, cipher.ny);

	delete[] trbxy;
}
