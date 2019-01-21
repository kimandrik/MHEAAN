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
#include "SerializationUtils.h"

Scheme::Scheme(SecretKey& secretKey, Ring& ring, bool isSerialized) : ring(ring), isSerialized(isSerialized) {
	addEncKey(secretKey);
	addMultKey(secretKey);
};


//----------------------------------------------------------------------------------
//   KEYS GENERATION
//----------------------------------------------------------------------------------


void Scheme::addEncKey(SecretKey& secretKey) {
	ZZ ax[N], bx[N];

	ring.sampleRLWE(ax, bx, secretKey.sx, logQQ);

	Key* key = new Key();

	ring.toNTT(key->rax, ax, nprimes);
	ring.toNTT(key->rbx, bx, nprimes);

	if(isSerialized) {
		string path = "serkey/ENCRYPTION.txt";
		SerializationUtils::writeKey(*key, path);
		serKeyMap.insert(pair<long, string>(ENCRYPTION, path));
		delete key;
	} else {
		keyMap.insert(pair<long, Key&>(ENCRYPTION, *key));
	}
}

void Scheme::addMultKey(SecretKey& secretKey) {
	ZZ sx2[N], ax[N], bx[N];

	ring.sampleRLWE(ax, bx, secretKey.sx, logQQ);

	ring.square(sx2, secretKey.sx, 1, Q);
	ring.leftShiftAndEqual(sx2, logQ, QQ);
	ring.addAndEqual(bx, sx2, QQ);

	Key* key = new Key();
	ring.toNTT(key->rax, ax, nprimes);
	ring.toNTT(key->rbx, bx, nprimes);

	if(isSerialized) {
		string path = "serkey/MULTIPLICATION.txt";
		SerializationUtils::writeKey(*key, path);
		serKeyMap.insert(pair<long, string>(MULTIPLICATION, path));
		delete key;
	} else {
		keyMap.insert(pair<long, Key&>(MULTIPLICATION, *key));
	}
}

void Scheme::addConjKey(SecretKey& secretKey) {
	ZZ sxcnj[N], ax[N], bx[N];

	ring.sampleRLWE(ax, bx, secretKey.sx, logQQ);

	ring.conjugate(sxcnj, secretKey.sx);
	ring.leftShiftAndEqual(sxcnj, logQ, QQ);
	ring.addAndEqual(bx, sxcnj, QQ);

	Key* key = new Key();
	ring.toNTT(key->rax, ax, nprimes);
	ring.toNTT(key->rbx, bx, nprimes);

	if(isSerialized) {
		string path = "serkey/CONJUGATION.txt";
		SerializationUtils::writeKey(*key, path);
		serKeyMap.insert(pair<long, string>(CONJUGATION, path));
		delete key;
	} else {
		keyMap.insert(pair<long, Key&>(CONJUGATION, *key));
	}
}

void Scheme::addLeftRotKey(SecretKey& secretKey, long r0, long r1) {
	ZZ sxrot[N], ax[N], bx[N];

	ring.sampleRLWE(ax, bx, secretKey.sx, logQQ);
	ring.leftRotate(sxrot, secretKey.sx, r0, r1);
	ring.leftShiftAndEqual(sxrot, logQ, QQ);
	ring.addAndEqual(bx, sxrot, QQ);

	Key* key = new Key();
	ring.toNTT(key->rax, ax, nprimes);
	ring.toNTT(key->rbx, bx, nprimes);

	if(isSerialized) {
		string path = "serkey/ROTATION_" + to_string(r0) + "_" + to_string(r1) + ".txt";
		SerializationUtils::writeKey(*key, path);
		serLeftRotKeyMap.insert(pair<pair<long, long>, string>({r0, r1}, path));
		delete key;
	} else {
		leftRotKeyMap.insert(pair<pair<long, long>, Key&>({r0, r1}, *key));
	}
}

void Scheme::addLeftX0RotKeys(SecretKey& secretKey) {
	for (long i = 1; i < N0h; i <<= 1) {
		if(leftRotKeyMap.find({i, 0}) == leftRotKeyMap.end() && serLeftRotKeyMap.find({i, 0}) == serLeftRotKeyMap.end()) {
			addLeftRotKey(secretKey, i, 0);
		}
	}
}

void Scheme::addLeftX1RotKeys(SecretKey& secretKey) {
	for (long i = 1; i < N1; i <<=1) {
		if(leftRotKeyMap.find({0, i}) == leftRotKeyMap.end() && serLeftRotKeyMap.find({0, i}) == serLeftRotKeyMap.end()) {
			addLeftRotKey(secretKey, 0, i);
		}
	}
}

void Scheme::addRightX0RotKeys(SecretKey& secretKey) {
	for (long i = 1; i < N0h; i <<=1) {
		long idx = N0h - i;
		if(leftRotKeyMap.find({idx, 0}) == leftRotKeyMap.end() && serLeftRotKeyMap.find({idx, 0}) == serLeftRotKeyMap.end()) {
			addLeftRotKey(secretKey, idx, 0);
		}
	}
}

void Scheme::addRightX1RotKeys(SecretKey& secretKey) {
	for (long i = 1; i < N1; i<<=1) {
		long idx = N1 - i;
		if(leftRotKeyMap.find({0, idx}) == leftRotKeyMap.end() && serLeftRotKeyMap.find({0, idx}) == serLeftRotKeyMap.end()) {
			addLeftRotKey(secretKey, 0, idx);
		}
	}
}

void Scheme::addBootContext(long logn0, long logn1, long logp) {
	if (bootContextMap.find({logn0, logn1}) == bootContextMap.end()) {
		long n0 = 1 << logn0;
		long logk0 = logn0 >> 1;
		long k0 = 1 << logk0;

		uint64_t** rpVec = new uint64_t*[n0];
		uint64_t** rpInvVec = new uint64_t*[n0];
		uint64_t* rp1 = NULL;
		uint64_t* rp2 = NULL;

		long* bndVec = new long[n0];
		long* bndInvVec = new long[n0];
		long bnd1 = 0;
		long bnd2 = 0;

		long np;
		complex<double>* pvals = new complex<double>[n0];
		ZZ pVec[N0];

		long gap0 = N0h >> logn0;
		long deg;
		for (long ki = 0; ki < n0; ki += k0) {
			for (long pos = ki; pos < ki + k0; ++pos) {
				for (long i = 0; i < n0 - pos; ++i) {
					deg = ((M0 - ring.gM0Pows[i + pos]) * i * gap0) % M0;
					pvals[i] = ring.ksiM0Pows[deg];
				}
				for (long i = n0 - pos; i < n0; ++i) {
					deg = ((M0 - ring.gM0Pows[i + pos - n0]) * i * gap0) % M0;
					pvals[i] = ring.ksiM0Pows[deg];
				}
				EvaluatorUtils::rightRotateAndEqual(pvals, n0, 1, ki, 0);
				ring.IEMBX0(pvals, n0);
				for (long i = 0, jd = N0h, id = 0; i < n0; ++i, jd += gap0, id += gap0) {
					pVec[id] = EvaluatorUtils::scaleUpToZZ(pvals[i].real(), logp);
					pVec[jd] = EvaluatorUtils::scaleUpToZZ(pvals[i].imag(), logp);
				}
				bndVec[pos] = ring.MaxBits(pVec, N0);
				np = ceil((logQ + bndVec[pos] + logN0 + 3)/(double)pbnd);
				rpVec[pos] = new uint64_t[np << logN0];
				ring.toNTTX0(rpVec[pos], pVec, np);
			}
		}

		for (long ki = 0; ki < n0; ki += k0) {
			for (long pos = ki; pos < ki + k0; ++pos) {
				for (long i = 0; i < n0 - pos; ++i) {
					deg = (ring.gM0Pows[i] * (i + pos) * gap0) % M0;
					pvals[i] = ring.ksiM0Pows[deg];
				}
				for (long i = n0 - pos; i < n0; ++i) {
					deg = (ring.gM0Pows[i] * (i + pos - n0) * gap0) % M0;
					pvals[i] = ring.ksiM0Pows[deg];
				}
				EvaluatorUtils::rightRotateAndEqual(pvals, n0, 1, ki, 0);
				ring.IEMBX0(pvals, n0);
				for (long i = 0, jd = N0h, id = 0; i < n0; ++i, jd += gap0, id += gap0) {
					pVec[id] = EvaluatorUtils::scaleUpToZZ(pvals[i].real(), logp);
					pVec[jd] = EvaluatorUtils::scaleUpToZZ(pvals[i].imag(), logp);
				}
				bndInvVec[pos] = ring.MaxBits(pVec, N0);
				np = ceil((logQ + bndInvVec[pos] + logN0 + 3)/(double)pbnd);
				rpInvVec[pos] = new uint64_t[np << logN0];
				ring.toNTTX0(rpInvVec[pos], pVec, np);
			}
		}

		delete[] pvals;

		BootContext* bootContext = new BootContext(rpVec, rpInvVec, rp1, rp2, bndVec, bndInvVec, bnd1, bnd2, logp);
		bootContextMap.insert(pair<pair<long, long>, BootContext&>({logn0, logn1}, *bootContext));
	}
}

void Scheme::addSqrMatContext(long logn, long logp) {
	if (sqrMatContextMap.find(logn) == sqrMatContextMap.end()) {
		long n = (1 << logn);

		Plaintext* msgvec = new Plaintext[n];
		double* tmp = new double[n * n]();
		for (long i = 0; i < n; ++i) {
			for (long j = 0; j < n; ++j) {
				tmp[j + (((j + n - i) % n) * n)] = 1.0;
			}
			encode(msgvec[i], tmp, n, n, logp);

			for (long j = 0; j < n; ++j) {
				tmp[j + (((j + n - i) % n) * n)] = 0.0;
			}
		}
		delete[] tmp;

		SqrMatContext* sqrMatContext = new SqrMatContext(msgvec);
		sqrMatContextMap.insert(pair<long, SqrMatContext&>(logn, *sqrMatContext));
	}
}

void Scheme::addBootKey(SecretKey& secretKey, long logn0, long logn1, long logp) {
	addBootContext(logn0, logn1, logp);

	addConjKey(secretKey);
	addLeftX0RotKeys(secretKey);
	addLeftX1RotKeys(secretKey);

	long logn0h = logn0 / 2;
	long k0 = 1 << logn0h;
	long m0 = 1 << (logn0 - logn0h);

	for (long i = 1; i < k0; ++i) {
		if(leftRotKeyMap.find({i,0}) == leftRotKeyMap.end() && serLeftRotKeyMap.find({i, 0}) == serLeftRotKeyMap.end()) {
			addLeftRotKey(secretKey, i, 0);
		}
	}

	for (long i = 1; i < m0; ++i) {
		long idx = i * k0;
		if(leftRotKeyMap.find({idx, 0}) == leftRotKeyMap.end() && serLeftRotKeyMap.find({idx, 0}) == serLeftRotKeyMap.end()) {
			addLeftRotKey(secretKey, idx, 0);
		}
	}

	long logn1h = logn1 / 2;
	long k1 = 1 << logn1h;
	long m1 = 1 << (logn1 - logn1h);

	for (long i = 1; i < k1; ++i) {
		if(leftRotKeyMap.find({0,i}) == leftRotKeyMap.end() && serLeftRotKeyMap.find({0, i}) == serLeftRotKeyMap.end()) {
			addLeftRotKey(secretKey, 0, i);
		}
	}

	for (long i = 1; i < m1; ++i) {
		long idx = i * k1;
		if(leftRotKeyMap.find({0, idx}) == leftRotKeyMap.end() && serLeftRotKeyMap.find({0, idx}) == serLeftRotKeyMap.end()) {
			addLeftRotKey(secretKey, 0, idx);
		}
	}
}

void Scheme::addSqrMatKeys(SecretKey& secretKey, long logn, long logp) {
	addSqrMatContext(logn, logp);
	long n = (1 << logn);
	for (long i = 1; i < n; ++i) {
		long idx = N0h - i;
		if(leftRotKeyMap.find({idx, 0}) == leftRotKeyMap.end() && serLeftRotKeyMap.find({idx, 0}) == serLeftRotKeyMap.end()) {
			addLeftRotKey(secretKey, idx, 0);
		}
	}
	addLeftX1RotKeys(secretKey);
}

void Scheme::addTransposeKeys(SecretKey& secretKey, long logn, long logp) {
	addSqrMatContext(logn, logp);
	long n = (1 << logn);
	for (long i = 1; i < n; ++i) {
		long i1 = N1 - i;
		if(leftRotKeyMap.find({i, i1}) == leftRotKeyMap.end() && serLeftRotKeyMap.find({i, i1}) == serLeftRotKeyMap.end()) {
			addLeftRotKey(secretKey, i, i1);
		}
	}
}


//----------------------------------------------------------------------------------
//   ENCRYPTION & DECRYPTION
//----------------------------------------------------------------------------------

void Scheme::encode(Plaintext& msg, complex<double>* vals, long n0, long n1, long logp) {
	ring.encode(msg.mx, vals, n0, n1, logp);
	msg.n0 = n0;
	msg.n1 = n1;
	msg.logp = logp;
}

void Scheme::encode(Plaintext& msg, double* vals, long n0, long n1, long logp) {
	ring.encode(msg.mx, vals, n0, n1, logp);
	msg.n0 = n0;
	msg.n1 = n1;
	msg.logp = logp;
}

void Scheme::rlwe(Ciphertext& res, long logq) {
	ZZ qQ = ring.qvec[logq + logQ];
	long* vx = new long[N];

	Key& key = isSerialized ? SerializationUtils::readKey(serKeyMap.at(ENCRYPTION)) : keyMap.at(ENCRYPTION);
	long np = ceil((1 + logQQ + logN + 3)/(double)pbnd);
	ring.sampleZO(vx);

	ring.multNTT(res.ax, vx, key.rax, np, qQ);
	ring.addGauss(res.ax, qQ);
	ring.rightShiftAndEqual(res.ax, logQ);

	ring.multNTT(res.bx, vx, key.rbx, np, qQ);
	ring.addGauss(res.bx, qQ);
	ring.rightShiftAndEqual(res.bx, logQ);
	delete[] vx;

	if(isSerialized) delete &key;
}

void Scheme::encrypt(Ciphertext& res, complex<double>* vals, long n0, long n1, long logp, long logq) {
	ZZ q = ring.qvec[logq];
	Plaintext msg;
	encryptZeros(res, n0, n1, logp, logq);
	encode(msg, vals, n0, n1, logp);
	addAndEqual(res, msg);
}

void Scheme::encrypt(Ciphertext& res, double* vals, long n0, long n1, long logp, long logq) {
	ZZ q = ring.qvec[logq];
	Plaintext msg;
	encryptZeros(res, n0, n1, logp, logq);
	encode(msg, vals, n0, n1, logp);
	addAndEqual(res, msg);
}

void Scheme::encryptSingle(Ciphertext& res, complex<double> val, long logp, long logq) {
	encryptZeros(res, 1, 1, logp, logq);
}

void Scheme::encryptSingle(Ciphertext& res, double val, long logp, long logq) {
	encryptZeros(res, 1, 1, logp, logq);
}

void Scheme::encryptZeros(Ciphertext& res, long n0, long n1, long logp, long logq) {
	rlwe(res, logq);
	res.n0 = n0;
	res.n1 = n1;
	res.logp = logp;
	res.logq = logq;
}

void Scheme::decryptMsg(Plaintext& msg, Ciphertext& cipher, SecretKey& secretKey) {
	ZZ q = ring.qvec[cipher.logq];
	long np = ceil((1 + cipher.logq + logN + 3)/(double)pbnd);
	ring.mult(msg.mx, cipher.ax, secretKey.sx, np, q);
	ring.addAndEqual(msg.mx, cipher.bx, q);
}

complex<double>* Scheme::decrypt(SecretKey& secretKey, Ciphertext& cipher) {
	Plaintext msg;
	decryptMsg(msg, cipher, secretKey);
	return ring.decode(msg.mx, cipher.n0, cipher.n1, cipher.logp, cipher.logq);
}

complex<double> Scheme::decryptSingle(SecretKey& secretKey, Ciphertext& cipher) {
	complex<double> res;
	Plaintext msg;
	decryptMsg(msg, cipher, secretKey);
	return res;
}


//----------------------------------------------------------------------------------
//   HOMOMORPHIC OPERATIONS
//----------------------------------------------------------------------------------


void Scheme::negate(Ciphertext& res, Ciphertext& cipher) {
	res.copyParams(cipher);
	ring.negate(res.ax, cipher.ax);
	ring.negate(res.bx, cipher.bx);
}

void Scheme::negateAndEqual(Ciphertext& cipher) {
	ring.negateAndEqual(cipher.ax);
	ring.negateAndEqual(cipher.bx);
}

void Scheme::add(Ciphertext& res, Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZ q = ring.qvec[cipher1.logq];
	res.copyParams(cipher1);
	ring.add(res.ax, cipher1.ax, cipher2.ax, q);
	ring.add(res.bx, cipher1.bx, cipher2.bx, q);
}

void Scheme::addAndEqual(Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZ q = ring.qvec[cipher1.logq];
	ring.addAndEqual(cipher1.ax, cipher2.ax, q);
	ring.addAndEqual(cipher1.bx, cipher2.bx, q);
}

void Scheme::addConst(Ciphertext& res, Ciphertext& cipher, double cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];
	res.copy(cipher);
	ZZ cnstZZ = logp < 0 ? -EvaluatorUtils::scaleUpToZZ(cnst, cipher.logp) : -EvaluatorUtils::scaleUpToZZ(cnst, logp);
	for (long i = 0; i < N; i+=N0) {
		AddMod(res.bx[i], res.bx[i], cnstZZ, q);
	}
}

void Scheme::addConst(Ciphertext& res, Ciphertext& cipher, RR& cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];
	res.copy(cipher);
	ZZ cnstZZ = logp < 0 ? -EvaluatorUtils::scaleUpToZZ(cnst, cipher.logp) : -EvaluatorUtils::scaleUpToZZ(cnst, logp);
	for (long i = 0; i < N; i+=N0) {
		AddMod(res.bx[i], res.bx[i], cnstZZ, q);
	}
}

void Scheme::addConstAndEqual(Ciphertext& cipher, double cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ cnstZZ = logp < 0 ? -EvaluatorUtils::scaleUpToZZ(cnst, cipher.logp) : -EvaluatorUtils::scaleUpToZZ(cnst, logp);
	for (long i = 0; i < N; i+=N0) {
		AddMod(cipher.bx[i], cipher.bx[i], cnstZZ, q);
	}
}

void Scheme::addConstAndEqual(Ciphertext& cipher, RR& cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ cnstZZ = logp < 0 ? -EvaluatorUtils::scaleUpToZZ(cnst, cipher.logp) : -EvaluatorUtils::scaleUpToZZ(cnst, logp);
	for (long i = 0; i < N; i+=N0) {
		AddMod(cipher.bx[i], cipher.bx[i], cnstZZ, q);
	}
}

void Scheme::add(Ciphertext& res, Ciphertext& cipher, Plaintext& msg) {
	ZZ q = ring.qvec[cipher.logq];
	res.copy(cipher);
	ring.addAndEqual(res.bx, msg.mx, q);
}

void Scheme::addAndEqual(Ciphertext& cipher, Plaintext& msg) {
	ZZ q = ring.qvec[cipher.logq];
	ring.addAndEqual(cipher.bx, msg.mx, q);
}

void Scheme::sub(Ciphertext& res, Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZ q = ring.qvec[cipher1.logq];
	res.copyParams(cipher1);
	ring.sub(res.ax, cipher1.ax, cipher2.ax, q);
	ring.sub(res.bx, cipher1.bx, cipher2.bx, q);
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

void Scheme::imult(Ciphertext& res, Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	res.copyParams(cipher);
	ring.multByMonomial(res.ax, cipher.ax, N0h, 0, q);
	ring.multByMonomial(res.bx, cipher.bx, N0h, 0, q);
}

void Scheme::idiv(Ciphertext& res, Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	res.copyParams(cipher);
	ring.multByMonomial(res.ax, cipher.ax, 3 * N0h, 0, q);
	ring.multByMonomial(res.bx, cipher.bx, 3 * N0h, 0, q);
}

void Scheme::imultAndEqual(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	ring.multByMonomialAndEqual(cipher.ax, N0h, 0, q);
	ring.multByMonomialAndEqual(cipher.bx, N0h, 0, q);
}

void Scheme::idivAndEqual(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	ring.multByMonomialAndEqual(cipher.ax, 3 * N0h, 0, q);
	ring.multByMonomialAndEqual(cipher.bx, 3 * N0h, 0, q);
}

void Scheme::mult(Ciphertext& res, Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZ q = ring.qvec[cipher1.logq]; // 2^1200
	ZZ qQ = ring.qvec[cipher1.logq + logQ]; // 2^2400

	long np = ceil((2 + cipher1.logq + cipher2.logq + logN + 3)/(double)pbnd);
	uint64_t* ra1 = new uint64_t[np << logN];
	uint64_t* rb1 = new uint64_t[np << logN];
	uint64_t* ra2 = new uint64_t[np << logN];
	uint64_t* rb2 = new uint64_t[np << logN];

	ring.toNTT(ra1, cipher1.ax, np);
	ring.toNTT(rb1, cipher1.bx, np);
	ring.toNTT(ra2, cipher2.ax, np);
	ring.toNTT(rb2, cipher2.bx, np);

	ZZ aax[N], bbx[N], abx[N];
	ring.multDNTT(aax, ra1, ra2, np, q);
	ring.multDNTT(bbx, rb1, rb2, np, q);

	ring.addNTTAndEqual(ra1, rb1, np);
	ring.addNTTAndEqual(ra2, rb2, np);
	ring.multDNTT(abx, ra1, ra2, np, q);
	delete[] ra1; delete[] ra2; delete[] rb1; delete[] rb2;

	np = ceil((cipher1.logq + logQQ + logN + 3)/(double)pbnd);
	uint64_t* raa = new uint64_t[np << logN];
	ring.toNTT(raa, aax, np);
	Key& key = isSerialized ? SerializationUtils::readKey(serKeyMap.at(MULTIPLICATION)) : keyMap.at(MULTIPLICATION);
	res.copyParams(cipher1);
	res.logp += cipher2.logp;
	ring.multDNTT(res.ax, raa, key.rax, np, qQ);
	ring.multDNTT(res.bx, raa, key.rbx, np, qQ);
	delete[] raa;
	if(isSerialized) delete &key;

	ring.rightShiftAndEqual(res.ax, logQ);
	ring.rightShiftAndEqual(res.bx, logQ);

	ring.addAndEqual(res.ax, abx, q);
	ring.subAndEqual(res.ax, bbx, q);
	ring.subAndEqual(res.ax, aax, q);
	ring.addAndEqual(res.bx, bbx, q);
}

void Scheme::multAndEqual(Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZ q = ring.qvec[cipher1.logq];
	ZZ qQ = ring.qvec[cipher1.logq + logQ];

	long np = ceil((2 + cipher1.logq + cipher2.logq + logN + 3)/(double)pbnd);
	uint64_t* ra1 = new uint64_t[np << logN];
	uint64_t* rb1 = new uint64_t[np << logN];
	uint64_t* ra2 = new uint64_t[np << logN];
	uint64_t* rb2 = new uint64_t[np << logN];
	ZZ* aax = new ZZ[N];
	ZZ* bbx = new ZZ[N];
	ZZ* abx = new ZZ[N];

	ring.toNTT(ra1, cipher1.ax, np);
	ring.toNTT(rb1, cipher1.bx, np);
	ring.toNTT(ra2, cipher2.ax, np);
	ring.toNTT(rb2, cipher2.bx, np);

	ring.multDNTT(aax, ra1, ra2, np, q);
	ring.multDNTT(bbx, rb1, rb2, np, q);

	ring.addNTTAndEqual(ra1, rb1, np);
	ring.addNTTAndEqual(ra2, rb2, np);
	ring.multDNTT(abx, ra1, ra2, np, q);
	delete[] ra1; delete[] ra2; delete[] rb1; delete[] rb2;

	np = ceil((cipher1.logq + logQQ + logN + 3)/(double)pbnd);
	uint64_t* raa = new uint64_t[np << logN];
	ring.toNTT(raa, aax, np);
	Key& key = isSerialized ? SerializationUtils::readKey(serKeyMap.at(MULTIPLICATION)) : keyMap.at(MULTIPLICATION);
	ring.multDNTT(cipher1.ax, raa, key.rax, np, qQ);
	ring.multDNTT(cipher1.bx, raa, key.rbx, np, qQ);
	delete[] raa;
	if(isSerialized) delete &key;

	ring.rightShiftAndEqual(cipher1.ax, logQ);
	ring.rightShiftAndEqual(cipher1.bx, logQ);

	ring.addAndEqual(cipher1.ax, abx, q);
	ring.subAndEqual(cipher1.ax, bbx, q);
	ring.subAndEqual(cipher1.ax, aax, q);
	ring.addAndEqual(cipher1.bx, bbx, q);
	delete[] aax; delete[] bbx; delete[] abx;
	cipher1.logp += cipher2.logp;
}

void Scheme::square(Ciphertext& res, Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ qQ = ring.qvec[cipher.logq + logQ];

	long np = ceil((2 * cipher.logq + logN + 3)/(double)pbnd);

	uint64_t* ra = new uint64_t[np << logN];
	uint64_t* rb = new uint64_t[np << logN];
	ring.toNTT(ra, cipher.ax, np);
	ring.toNTT(rb, cipher.bx, np);

	ZZ aax[N], bbx[N], abx[N];
	ring.squareNTT(bbx, rb, np, q);
	ring.squareNTT(aax, ra, np, q);
	ring.multDNTT(abx, ra, rb, np, q);
	ring.leftShiftAndEqual(abx, 1, q);
	delete[] ra; delete[] rb;
	res.copyParams(cipher);
	res.logp *= 2;
	np = ceil((cipher.logq + logQQ + logN + 3)/(double)pbnd);
	uint64_t* raa = new uint64_t[np << logN];
	ring.toNTT(raa, aax, np);
	Key& key = isSerialized ? SerializationUtils::readKey(serKeyMap.at(MULTIPLICATION)) : keyMap.at(MULTIPLICATION);
	ring.multDNTT(res.ax, raa, key.rax, np, qQ);
	ring.multDNTT(res.bx, raa, key.rbx, np, qQ);
	delete[] raa;
	if(isSerialized) delete &key;

	ring.rightShiftAndEqual(res.ax, logQ);
	ring.rightShiftAndEqual(res.bx, logQ);

	ring.addAndEqual(res.ax, abx, q);
	ring.addAndEqual(res.bx, bbx, q);
}

void Scheme::squareAndEqual(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ qQ = ring.qvec[cipher.logq + logQ];
	long np = ceil((2 * cipher.logq + logN + 3)/(double)pbnd);
	uint64_t* ra = new uint64_t[np << logN];
	uint64_t* rb = new uint64_t[np << logN];
	ring.toNTT(ra, cipher.ax, np);
	ring.toNTT(rb, cipher.bx, np);

	ZZ aax[N], bbx[N], abx[N];
	ring.squareNTT(bbx, rb, np, q);
	ring.multDNTT(abx, ra, rb, np, q);
	ring.leftShiftAndEqual(abx, 1, q);
	ring.squareNTT(aax, ra, np, q);
	delete[] ra; delete[] rb;

	np = ceil((cipher.logq + logQQ + logN + 3)/(double)pbnd);
	uint64_t* raa = new uint64_t[np << logN];
	ring.toNTT(raa, aax, np);
	Key& key = isSerialized ? SerializationUtils::readKey(serKeyMap.at(MULTIPLICATION)) : keyMap.at(MULTIPLICATION);
	ring.multDNTT(cipher.ax, raa, key.rax, np, qQ);
	ring.multDNTT(cipher.bx, raa, key.rbx, np, qQ);
	delete[] raa;
	if(isSerialized) delete &key;

	ring.rightShiftAndEqual(cipher.ax, logQ);
	ring.rightShiftAndEqual(cipher.bx, logQ);
	ring.addAndEqual(cipher.ax, abx, q);
	ring.addAndEqual(cipher.bx, bbx, q);
	cipher.logp *= 2;
}

void Scheme::multConst(Ciphertext& res, Ciphertext& cipher, RR& cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ cnstZZ = EvaluatorUtils::scaleUpToZZ(cnst, logp);
	res.copyParams(cipher);
	ring.multByConst(res.ax, cipher.ax, cnstZZ, q);
	ring.multByConst(res.bx, cipher.bx, cnstZZ, q);
	res.logp += logp;
}


void Scheme::multConst(Ciphertext& res, Ciphertext& cipher, double cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ cnstZZ = EvaluatorUtils::scaleUpToZZ(cnst, logp);
	res.copyParams(cipher);
	ring.multByConst(res.ax, cipher.ax, cnstZZ, q);
	ring.multByConst(res.bx, cipher.bx, cnstZZ, q);
	res.logp += logp;
}

void Scheme::multConst(Ciphertext& res, Ciphertext& cipher, complex<double> cnst, long logp) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ axi[N], bxi[N];
	ring.multByMonomial(axi, cipher.ax, N0h, 0, q);
	ring.multByMonomial(bxi, cipher.bx, N0h, 0, q);
	ZZ cnstrZZ = EvaluatorUtils::scaleUpToZZ(cnst.real(), logp);
	ZZ cnstiZZ = EvaluatorUtils::scaleUpToZZ(cnst.imag(), logp);
	res.copyParams(cipher);
	ring.multByConst(res.ax, cipher.ax, cnstrZZ, q);
	ring.multByConst(res.bx, cipher.bx, cnstrZZ, q);
	ring.multByConstAndEqual(axi, cnstiZZ, q);
	ring.multByConstAndEqual(bxi, cnstiZZ, q);
	ring.addAndEqual(res.ax, axi, q);
	ring.addAndEqual(res.bx, bxi, q);
	res.logp += logp;
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
	ZZ axi[N], bxi[N];
	ring.multByMonomial(axi, cipher.ax, N0h, 0, q);
	ring.multByMonomial(bxi, cipher.bx, N0h, 0, q);
	ring.multByConstAndEqual(axi, cnstiZZ, q);
	ring.multByConstAndEqual(bxi, cnstiZZ, q);
	ring.multByConstAndEqual(cipher.ax, cnstrZZ, q);
	ring.multByConstAndEqual(cipher.bx, cnstrZZ, q);
	ring.addAndEqual(cipher.ax, axi, q);
	ring.addAndEqual(cipher.bx, bxi, q);
	cipher.logp += logp;
}

void Scheme::multPolyX0(Ciphertext& res, Ciphertext& cipher, ZZ* poly, long logp) {
	ZZ q = ring.qvec[cipher.logq];
	res.copy(cipher);
	long bnd = ring.MaxBits(poly, N0);
	long np = ceil((cipher.logq + bnd + logN0 + 3)/(double)pbnd);
	uint64_t* rpoly = new uint64_t[np << logN0];
	ring.toNTTX0(rpoly, poly, np);
	ring.multNTTX0AndEqual(res.ax, rpoly, np, q);
	ring.multNTTX0AndEqual(res.bx, rpoly, np, q);
	delete[] rpoly;
	res.logp += logp;
}

void Scheme::multPolyX0AndEqual(Ciphertext& cipher, ZZ* poly, long logp) {
	ZZ q = ring.qvec[cipher.logq];
	long bnd = ring.MaxBits(poly, N0);
	long np = ceil((cipher.logq + bnd + logN0 + 3)/(double)pbnd);
	uint64_t* rpoly = new uint64_t[np << logN0];
	ring.toNTTX0(rpoly, poly, np);
	ring.multNTTX0AndEqual(cipher.ax, rpoly, np, q);
	ring.multNTTX0AndEqual(cipher.bx, rpoly, np, q);
	delete[] rpoly;
	cipher.logp += logp;
}

void Scheme::multPolyNTTX0(Ciphertext& res, Ciphertext& cipher, uint64_t* rpoly, long bnd, long logp) {
	ZZ q = ring.qvec[cipher.logq];
	res.copy(cipher);
	long np = ceil((cipher.logq + bnd + logN0 + 3)/(double)pbnd);
	ring.multNTTX0AndEqual(res.ax, rpoly, np, q);
	ring.multNTTX0AndEqual(res.bx, rpoly, np, q);
	res.logp += logp;
}

void Scheme::multPolyNTTX0AndEqual(Ciphertext& cipher, uint64_t* rpoly, long bnd, long logp) {
	ZZ q = ring.qvec[cipher.logq];
	long np = ceil((cipher.logq + bnd + logN0 + 3)/(double)pbnd);
	ring.multNTTX0AndEqual(cipher.ax, rpoly, np, q);
	ring.multNTTX0AndEqual(cipher.bx, rpoly, np, q);
	cipher.logp += logp;
}

void Scheme::multPolyX1(Ciphertext& res, Ciphertext& cipher, ZZ* rpoly, ZZ* ipoly, long logp) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ axi[N], bxi[N];
	ring.multByMonomial(axi, cipher.ax, N0h, 0, q);
	ring.multByMonomial(bxi, cipher.bx, N0h, 0, q);
	long bnd = ring.MaxBits(ipoly, N1);
	long np = ceil((cipher.logq + bnd + logN1 + 3)/(double)pbnd);
	uint64_t* ripoly = new uint64_t[np << logN];
	ring.toNTT(ripoly, ipoly, np);
	ring.multNTTX1AndEqual(axi, ripoly, np, q);
	ring.multNTTX1AndEqual(bxi, ripoly, np, q);
	delete[] ripoly;
	res.copy(cipher);
	bnd = ring.MaxBits(rpoly, N1);
	np = ceil((cipher.logq + bnd + logN1 + 3)/(double)pbnd);
	uint64_t* rrpoly = new uint64_t[np << logN];
	ring.toNTT(rrpoly, rpoly, np);
	ring.multNTTX1AndEqual(res.ax, rrpoly, np, q);
	ring.multNTTX1AndEqual(res.bx, rrpoly, np, q);
	delete[] rrpoly;
	ring.addAndEqual(res.ax, axi, q);
	ring.addAndEqual(res.bx, bxi, q);
	res.logp += logp;
}

void Scheme::multPolyX1AndEqual(Ciphertext& cipher, ZZ* rpoly, ZZ* ipoly, long logp) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ axi[N], bxi[N];
	ring.multByMonomial(axi, cipher.ax, N0h, 0, q);
	ring.multByMonomial(bxi, cipher.bx, N0h, 0, q);
	long bnd = ring.MaxBits(rpoly, N1);
	long np = ceil((cipher.logq + bnd + logN1 + 3)/(double)pbnd);
	uint64_t* rrpoly = new uint64_t[np << logN];
	ring.toNTT(rrpoly, rpoly, np);
	ring.multNTTX1AndEqual(cipher.ax, rrpoly, np, q);
	ring.multNTTX1AndEqual(cipher.bx, rrpoly, np, q);
	delete[] rrpoly;
	bnd = ring.MaxBits(ipoly, N1);
	np = ceil((cipher.logq + bnd + logN1 + 3)/(double)pbnd);
	uint64_t* ripoly = new uint64_t[np << logN];
	ring.toNTT(ripoly, ipoly, np);
	ring.multNTTX1AndEqual(axi, ripoly, np, q);
	ring.multNTTX1AndEqual(bxi, ripoly, np, q);
	delete[] ripoly;
	ring.addAndEqual(cipher.ax, axi, q);
	ring.addAndEqual(cipher.bx, bxi, q);
	cipher.logp += logp;
}

void Scheme::mult(Ciphertext& res, Ciphertext& cipher, Plaintext& msg) {
	ZZ q = ring.qvec[cipher.logq];

	long bnd = ring.MaxBits(msg.mx, N);
	long np = ceil((cipher.logq + bnd + logN + 3)/(double)pbnd);
	uint64_t* rpoly = new uint64_t[np << logN];
	ring.toNTT(rpoly, msg.mx, np);
	res.copy(cipher);
	ring.multNTTAndEqual(res.ax, rpoly, np, q);
	ring.multNTTAndEqual(res.bx, rpoly, np, q);
	delete[] rpoly;
	res.logp += msg.logp;
}

void Scheme::multAndEqual(Ciphertext& cipher, Plaintext& msg) {
	ZZ q = ring.qvec[cipher.logq];
	long bnd = ring.MaxBits(msg.mx, N);
	long np = ceil((cipher.logq + bnd + logN + 3)/(double)pbnd);
	uint64_t* rpoly = new uint64_t[np << logN];
	ring.toNTT(rpoly, msg.mx, np);
	ring.multNTTAndEqual(cipher.ax, rpoly, np, q);
	ring.multNTTAndEqual(cipher.bx, rpoly, np, q);
	delete[] rpoly;
	cipher.logp += msg.logp;
}


void Scheme::multPolyNTT(Ciphertext& res, Ciphertext& cipher, uint64_t* rpoly, long bnd, long logp) {
	ZZ q = ring.qvec[cipher.logq];
	res.copy(cipher);
	long np = ceil((cipher.logq + bnd + logN + 3)/(double)pbnd);
	ring.multNTTAndEqual(res.ax, rpoly, np, q);
	ring.multNTTAndEqual(res.bx, rpoly, np, q);
	res.logp += logp;
}

void Scheme::multPolyNTTAndEqual(Ciphertext& cipher, uint64_t* rpoly, long bnd, long logp) {
	ZZ q = ring.qvec[cipher.logq];
	long np = ceil((cipher.logq + bnd + logN + 3)/(double)pbnd);
	ring.multNTTAndEqual(cipher.ax, rpoly, np, q);
	ring.multNTTAndEqual(cipher.bx, rpoly, np, q);
	cipher.logp += logp;
}

void Scheme::multByMonomial(Ciphertext& res, Ciphertext& cipher, const long d0, const long d1) {
	ZZ q = ring.qvec[cipher.logq];
	res.copyParams(cipher);
	ring.multByMonomial(res.ax, cipher.ax, d0, d1, q);
	ring.multByMonomial(res.bx, cipher.bx, d0, d1, q);
}

void Scheme::multByMonomialAndEqual(Ciphertext& cipher, const long d0, const long d1) {
	ZZ q = ring.qvec[cipher.logq];
	ring.multByMonomialAndEqual(cipher.ax, d0, d1, q);
	ring.multByMonomialAndEqual(cipher.bx, d0, d1, q);
}

void Scheme::multPo2(Ciphertext& res, Ciphertext& cipher, long bits) {
	ZZ q = ring.qvec[cipher.logq];
	res.copyParams(cipher);
	ring.leftShift(res.ax, cipher.ax, bits, q);
	ring.leftShift(res.bx, cipher.bx, bits, q);
}

void Scheme::multPo2AndEqual(Ciphertext& cipher, long bits) {
	ZZ q = ring.qvec[cipher.logq];
	ring.leftShiftAndEqual(cipher.ax, bits, q);
	ring.leftShiftAndEqual(cipher.bx, bits, q);
}

void Scheme::divPo2(Ciphertext& res, Ciphertext& cipher, long logd) {
	res.copyParams(cipher);
	ring.rightShift(res.ax, cipher.ax, logd);
	ring.rightShift(res.bx, cipher.bx, logd);
	res.logq -= logd;
}

void Scheme::divPo2AndEqual(Ciphertext& cipher, long logd) {
	ring.rightShiftAndEqual(cipher.ax, logd);
	ring.rightShiftAndEqual(cipher.bx, logd);
	cipher.logq -= logd;
}


//----------------------------------------------------------------------------------
//   RESCALING & MODULUS DOWN
//----------------------------------------------------------------------------------


void Scheme::reScaleBy(Ciphertext& res, Ciphertext& cipher, long dlogq) {
	res.copyParams(cipher);
	ring.rightShift(res.ax, cipher.ax, dlogq);
	ring.rightShift(res.bx, cipher.bx, dlogq);
	res.logq -= dlogq;
	res.logp -= dlogq;
}

void Scheme::reScaleTo(Ciphertext& res, Ciphertext& cipher, long logq) {
	long dlogq = cipher.logq - logq;
	reScaleBy(res, cipher, dlogq);
}

void Scheme::reScaleByAndEqual(Ciphertext& cipher, long dlogq) {
	ring.rightShiftAndEqual(cipher.ax, dlogq);
	ring.rightShiftAndEqual(cipher.bx, dlogq);
	cipher.logq -= dlogq;
	cipher.logp -= dlogq;
}

void Scheme::reScaleToAndEqual(Ciphertext& cipher, long logq) {
	long dlogq = cipher.logq - logq;
	reScaleByAndEqual(cipher, dlogq);
}

void Scheme::modDownBy(Ciphertext& res, Ciphertext& cipher, long dlogq) {
	long logq = cipher.logq - dlogq;
	modDownTo(res, cipher, logq);
}

void Scheme::modDownByAndEqual(Ciphertext& cipher, long dlogq) {
	long logq = cipher.logq - dlogq;
	modDownToAndEqual(cipher, logq);
}

void Scheme::modDownTo(Ciphertext& res, Ciphertext& cipher, long logq) {
	ZZ q = ring.qvec[logq];
	res.copyParams(cipher);
	ring.mod(res.ax, cipher.ax, q);
	ring.mod(res.bx, cipher.bx, q);
	res.logq = logq;
}

void Scheme::modDownToAndEqual(Ciphertext& cipher, long logq) {
	ZZ q = ring.qvec[logq];
	ring.modAndEqual(cipher.ax, q);
	ring.modAndEqual(cipher.bx, q);
	cipher.logq = logq;
}


//----------------------------------------------------------------------------------
//   ROTATIONS & CONJUGATIONS
//----------------------------------------------------------------------------------


void Scheme::leftRotate(Ciphertext& res, Ciphertext& cipher, long r0, long r1) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ qQ = ring.qvec[cipher.logq + logQ];

	ZZ axrot[N], bxrot[N];

	ring.leftRotate(axrot, cipher.ax, r0, r1);
	ring.leftRotate(bxrot, cipher.bx, r0, r1);
	res.copyParams(cipher);
	long np = ceil((cipher.logq + logQQ + logN + 3)/(double)pbnd);
	uint64_t* rarot = new uint64_t[np << logN];
	ring.toNTT(rarot, axrot, np);
	Key& key = isSerialized ? SerializationUtils::readKey(serLeftRotKeyMap.at({r0, r1})) : leftRotKeyMap.at({r0, r1});
	ring.multDNTT(res.ax, rarot, key.rax, np, qQ);
	ring.multDNTT(res.bx, rarot, key.rbx, np, qQ);
	delete[] rarot;
	if(isSerialized) delete &key;

	ring.rightShiftAndEqual(res.ax, logQ);
	ring.rightShiftAndEqual(res.bx, logQ);

	ring.addAndEqual(res.bx, bxrot, q);
}

void Scheme::rightRotate(Ciphertext& res, Ciphertext& cipher, long r0, long r1) {
	long rr0 = r0 == 0 ? 0 : N0h - r0;
	long rr1 = r1 == 0 ? 0 : N1 - r1;
	leftRotate(res, cipher, rr0, rr1);
}

void Scheme::leftRotateAndEqual(Ciphertext& cipher, long r0, long r1) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ qQ = ring.qvec[cipher.logq + logQ];

	ZZ axrot[N], bxrot[N];

	ring.leftRotate(axrot, cipher.ax, r0, r1);
	ring.leftRotate(bxrot, cipher.bx, r0, r1);

	Key& key = isSerialized ? SerializationUtils::readKey(serLeftRotKeyMap.at({r0, r1})) : leftRotKeyMap.at({r0, r1});

	long np = ceil((cipher.logq + logQQ + logN + 3)/(double)pbnd);
	uint64_t* rarot = new uint64_t[np << logN];
	ring.toNTT(rarot, axrot, np);
	ring.multDNTT(cipher.bx, rarot, key.rbx, np, qQ);
	ring.multDNTT(cipher.ax, rarot, key.rax, np, qQ);
	delete[] rarot;
	if(isSerialized) delete &key;

	ring.rightShiftAndEqual(cipher.ax, logQ);
	ring.rightShiftAndEqual(cipher.bx, logQ);

	ring.addAndEqual(cipher.bx, bxrot, q);
}

void Scheme::rightRotateAndEqual(Ciphertext& cipher, long r0, long r1) {
	long rr0 = r0 == 0 ? 0 : N0h - r0;
	long rr1 = r1 == 0 ? 0 : N1 - r1;
	leftRotateAndEqual(cipher, rr0, rr1);
}

void Scheme::conjugate(Ciphertext& res, Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ qQ = ring.qvec[cipher.logq + logQ];

	ZZ axcnj[N], bxcnj[N];
	ring.conjugate(axcnj, cipher.ax);
	ring.conjugate(bxcnj, cipher.bx);

	long np = ceil((cipher.logq + logQQ + logN + 3)/(double)pbnd);
	uint64_t* racnj = new uint64_t[np << logN];
	ring.toNTT(racnj, axcnj, np);
	res.copyParams(cipher);
	Key& key = isSerialized ? SerializationUtils::readKey(serKeyMap.at(CONJUGATION)) : keyMap.at(CONJUGATION);
	ring.multDNTT(res.ax, racnj, key.rax, np, qQ);
	ring.multDNTT(res.bx, racnj, key.rbx, np, qQ);
	delete[] racnj;
	if(isSerialized) delete &key;

	ring.rightShiftAndEqual(res.ax, logQ);
	ring.rightShiftAndEqual(res.bx, logQ);

	ring.addAndEqual(res.bx, bxcnj, q);
}

void Scheme::conjugateAndEqual(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ qQ = ring.qvec[cipher.logq + logQ];

	ZZ axcnj[N], bxcnj[N];
	ring.conjugate(axcnj, cipher.ax);
	ring.conjugate(bxcnj, cipher.bx);

	Key& key = isSerialized ? SerializationUtils::readKey(serKeyMap.at(CONJUGATION)) : keyMap.at(CONJUGATION);

	long np = ceil((cipher.logq + logQQ + logN + 3)/(double)pbnd);
	uint64_t* racnj = new uint64_t[np << logN];
	ring.toNTT(racnj, axcnj, np);
	ring.multDNTT(cipher.ax, racnj, key.rax, np, qQ);
	ring.multDNTT(cipher.bx, racnj, key.rbx, np, qQ);
	delete[] racnj;
	if(isSerialized) delete &key;

	ring.rightShiftAndEqual(cipher.ax, logQ);
	ring.rightShiftAndEqual(cipher.bx, logQ);

	ring.addAndEqual(cipher.bx, bxcnj, q);
}


//----------------------------------------------------------------------------------
//   BOOTSTRAPPING
//----------------------------------------------------------------------------------


void Scheme::normalizeAndEqual(Ciphertext& cipher) {
	ZZ q = ring.qvec[cipher.logq];
	ZZ qh = ring.qvec[cipher.logq - 1];
	for (long i = 0; i < N; ++i) {
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
	mutex m;

	long n0 = cipher.n0;
	long n1 = cipher.n1;

	long logn0 = log2(n0);
	long logn1 = log2(n1);

	long logk0 = logn0 / 2;
	long k0 = 1 << logk0;

	Ciphertext* rotvec = new Ciphertext[k0];
	NTL_EXEC_RANGE(k0, first, last);
	for (long j = first; j < last; ++j) {
		rotvec[j].copy(cipher);
		if(j > 0) leftRotateAndEqual(rotvec[j], j, 0);
	}
	NTL_EXEC_RANGE_END;

	BootContext& bootContext = ring.bootContextMap.at({logn0, logn1});
	cipher.free();
	cipher.logp += bootContext.logp;

	Ciphertext aux;
	aux.copyParams(cipher);
	for (long ki = 0; ki < n0; ki += k0) {
		NTL_EXEC_RANGE(k0, first, last);
		for (long j = first; j < last; ++j) {
			Ciphertext tmp(rotvec[j]);
			multPolyNTTX0AndEqual(tmp, bootContext.rpxVec[j + ki], bootContext.bndVec[j + ki], bootContext.logp);
			m.lock();
			addAndEqual(aux, tmp);
			m.unlock();
		}
		NTL_EXEC_RANGE_END;
		if(ki > 0) leftRotateAndEqual(aux, ki, 0);
		addAndEqual(cipher, aux);
		aux.free();
	}

	reScaleByAndEqual(cipher, bootContext.logp);

	delete[] rotvec;
}

void Scheme::coeffToSlotX1AndEqual(Ciphertext& cipher) {
	mutex m;

	long n0 = cipher.n0;
	long n1 = cipher.n1;

	long logn0 = log2(n0);
	long logn1 = log2(n1);

	long logk1 = logn1 / 2;
	long k1 = 1 << logk1;

	Ciphertext rot(cipher);
	Ciphertext aux;
	for (long i = 1; i < N1; i <<= 1) {
		leftRotate(aux, rot, 0, i);
		addAndEqual(rot, aux);
	}
	aux.free();

	Ciphertext* rotvec = new Ciphertext[k1];
	NTL_EXEC_RANGE(k1, first, last);
	for (long j = first; j < last; ++j) {
		rotvec[j].copy(cipher);
		if(j > 0) leftRotateAndEqual(rotvec[j], 0, j);
	}
	NTL_EXEC_RANGE_END;

	BootContext& bootContext = ring.bootContextMap.at({logn0, logn1});
	cipher.free();
	cipher.logp += bootContext.logp;

	aux.copyParams(cipher);
	for (long ki = 0; ki < n1; ki += k1) {
		NTL_EXEC_RANGE(k1, first, last);
		for (long j = first; j < last; ++j) {
			Ciphertext tmp(rotvec[j]);
			complex<double> cnst = conj(ring.dftM1Pows[logn1][j + ki]) * (double)n1/(double)M1;
			multConstAndEqual(tmp, cnst, bootContext.logp);
			m.lock();
			addAndEqual(aux, tmp);
			m.unlock();
		}
		NTL_EXEC_RANGE_END;

		if(ki > 0) leftRotateAndEqual(aux, 0, ki);
		addAndEqual(cipher, aux);
		aux.free();
	}

	multConstAndEqual(rot, (double)n1/(double)M1, bootContext.logp);
	subAndEqual(cipher, rot);
	reScaleByAndEqual(cipher, bootContext.logp);

	delete[] rotvec;
}

void Scheme::coeffToSlotAndEqual(Ciphertext& cipher) {
	coeffToSlotX1AndEqual(cipher);
	coeffToSlotX0AndEqual(cipher);
}

void Scheme::slotToCoeffX0AndEqual(Ciphertext& cipher) {
	mutex m;

	long n0 = cipher.n0;
	long n1 = cipher.n1;

	long logn0 = log2(n0);
	long logn1 = log2(n1);

	long logk0 = logn0 / 2;
	long k0 = 1 << logk0;

	Ciphertext* rotvec = new Ciphertext[k0];
	NTL_EXEC_RANGE(k0, first, last);
	for (long j = first; j < last; ++j) {
		rotvec[j].copy(cipher);
		if(j > 0) leftRotateAndEqual(rotvec[j], j, 0);
	}
	NTL_EXEC_RANGE_END;

	BootContext& bootContext = ring.bootContextMap.at({logn0, logn1});
	cipher.free();
	cipher.logp += bootContext.logp;

	Ciphertext aux;
	aux.copyParams(cipher);
	for (long ki = 0; ki < n0; ki+=k0) {
		NTL_EXEC_RANGE(k0, first, last);
		for (long j = first; j < last; ++j) {
			Ciphertext tmp(rotvec[j]);
			multPolyNTTX0AndEqual(tmp, bootContext.rpxInvVec[j + ki], bootContext.bndInvVec[j + ki], bootContext.logp);
			m.lock();
			addAndEqual(aux, tmp);
			m.unlock();
		}
		NTL_EXEC_RANGE_END;
		if(ki > 0) leftRotateAndEqual(aux, ki, 0);
		addAndEqual(cipher, aux);
		aux.free();
	}

	reScaleByAndEqual(cipher, bootContext.logp);

	delete[] rotvec;
}

void Scheme::slotToCoeffX1AndEqual(Ciphertext& cipher) {
	mutex m;

	long n0 = cipher.n0;
	long n1 = cipher.n1;

	long logn0 = log2(n0);
	long logn1 = log2(n1);

	long logk1 = logn1 / 2;
	long k1 = 1 << logk1;

	Ciphertext* rotvec = new Ciphertext[k1];
	NTL_EXEC_RANGE(k1, first, last);
	for (long j = first; j < last; ++j) {
		rotvec[j].copy(cipher);
		if(j > 0) leftRotateAndEqual(rotvec[j], 0, j);
	}
	NTL_EXEC_RANGE_END;

	BootContext& bootContext = ring.bootContextMap.at({logn0, logn1});
	cipher.free();
	cipher.logp += bootContext.logp;

	Ciphertext aux;
	aux.copyParams(cipher);
	for (long ki = 0; ki < n1; ki+=k1) {
		NTL_EXEC_RANGE(k1, first, last);
		for (long j = first; j < last; ++j) {
			complex<double> cnst = ring.dftM1Pows[logn1][n1-j-ki];
			Ciphertext tmp(rotvec[j]);
			multConstAndEqual(tmp, cnst, bootContext.logp);
			m.lock();
			addAndEqual(aux, tmp);
			m.unlock();
		}
		NTL_EXEC_RANGE_END;

		if(ki > 0) leftRotateAndEqual(aux, 0, ki);
		addAndEqual(cipher, aux);
		aux.free();
	}

	reScaleByAndEqual(cipher, bootContext.logp);

	delete[] rotvec;
}

void Scheme::slotToCoeffAndEqual(Ciphertext& cipher) {
	slotToCoeffX0AndEqual(cipher);
	slotToCoeffX1AndEqual(cipher);
}

void Scheme::exp2piAndEqual(Ciphertext& cipher, long logp) {
	Ciphertext cipher2(cipher);
	squareAndEqual(cipher2);
	reScaleByAndEqual(cipher2, logp);

	Ciphertext cipher4(cipher2);
	squareAndEqual(cipher4);
	reScaleByAndEqual(cipher4, logp);

	RR c = 1/(2*Pi);
	Ciphertext cipher01(cipher);
	addConstAndEqual(cipher01, c, logp);

	c = 2*Pi;
	multConstAndEqual(cipher01, c, logp);
	reScaleByAndEqual(cipher01, logp);

	c = 3/(2*Pi);
	Ciphertext cipher23(cipher);
	addConstAndEqual(cipher23, c, logp);

	c = 4*Pi*Pi*Pi/3;
	multConstAndEqual(cipher23, c, logp);
	reScaleByAndEqual(cipher23, logp);

	multAndEqual(cipher23, cipher2);
	reScaleByAndEqual(cipher23, logp);

	addAndEqual(cipher23, cipher01);

	c = 5/(2*Pi);
	Ciphertext cipher45(cipher);
	addConstAndEqual(cipher45, c, logp);

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

void Scheme::removeIPartAndEqual(Ciphertext& cipher, long logT, long logI) {
	long logn0 = log2(cipher.n0);
	long logn1 = log2(cipher.n1);

	BootContext& bootContext = ring.bootContextMap.at({logn0, logn1});

	Ciphertext aux(cipher);
	Ciphertext cimag(cipher);
	conjugateAndEqual(aux);
	subAndEqual(cimag, aux);
	addAndEqual(cipher, aux);
	imultAndEqual(cipher);

	divPo2AndEqual(cipher, logT + logn0 + logn1 + 1);
	divPo2AndEqual(cimag, logT + logn0 + logn1 + 1);

	exp2piAndEqual(cipher, bootContext.logp);
	exp2piAndEqual(cimag, bootContext.logp);

	for (long i = 0; i < logI + logT; ++i) {
		squareAndEqual(cipher);
		squareAndEqual(cimag);
		reScaleByAndEqual(cipher, bootContext.logp);
		reScaleByAndEqual(cimag, bootContext.logp);
	}
	conjugate(aux, cimag);
	subAndEqual(cimag, aux);
	conjugate(aux, cipher);
	subAndEqual(cipher, aux);
	imultAndEqual(cipher);
	subAndEqual2(cimag, cipher);

	RR c = 0.25/Pi;
	multConstAndEqual(cipher, c, bootContext.logp);
	reScaleByAndEqual(cipher, bootContext.logp + logI + logN1 - logn1);
}

void Scheme::bootstrapX0AndEqual(Ciphertext& cipher, long logq, long logQ, long logT, long logI) {
	long logn0 = log2(cipher.n0);
	long logp = cipher.logp;

	modDownToAndEqual(cipher, logq);
	normalizeAndEqual(cipher);

	cipher.logq = logQ;
	cipher.logp = logq + logI;
	Ciphertext rot;
	for (long i = logn0; i < logN0h; ++i) {
		leftRotate(rot, cipher, (1 << i), 0);
		addAndEqual(cipher, rot);
	}

	coeffToSlotX0AndEqual(cipher);
	divPo2AndEqual(cipher, logN0h);
	removeIPartAndEqual(cipher, logT, logI);
	slotToCoeffX0AndEqual(cipher);
	cipher.logp = logp;
}

void Scheme::bootstrapX1AndEqual(Ciphertext& cipher, long logq, long logQ, long logT, long logI) {
	long logn1 = log2(cipher.n1);
	long logp = cipher.logp;

	modDownToAndEqual(cipher, logq);
	normalizeAndEqual(cipher);

	cipher.logq = logQ;
	cipher.logp = logq + logI;
	Ciphertext rot;
	for (long i = logn1; i < logN1; ++i) {
		leftRotate(rot, cipher, 0, (1 << i));
		addAndEqual(cipher, rot);
	}

	coeffToSlotX1AndEqual(cipher);
	divPo2AndEqual(cipher, logN1);
	removeIPartAndEqual(cipher, logT, logI);
	slotToCoeffX1AndEqual(cipher);
	cipher.logp = logp;
}

void Scheme::bootstrapAndEqual(Ciphertext& cipher, long logq, long logQ, long logT, long logI) {
	long logn0 = log2(cipher.n0);
	long logn1 = log2(cipher.n1);
	long logp = cipher.logp;

	modDownToAndEqual(cipher, logq);
	normalizeAndEqual(cipher);

	cipher.logq = logQ;
	cipher.logp = logq + logI;
	Ciphertext rot;
	for (long i = logn0; i < logN0h; ++i) {
		leftRotate(rot, cipher, (1 << i), 0);
		addAndEqual(cipher, rot);
	}
	for (long i = logn1; i < logN1; ++i) {
		leftRotate(rot, cipher, 0, (1 << i));
		addAndEqual(cipher, rot);
	}

	coeffToSlotAndEqual(cipher);
	divPo2AndEqual(cipher, logN0h + logN1);
	removeIPartAndEqual(cipher, logT, logI);
	slotToCoeffAndEqual(cipher);
	cipher.logp = logp;
}
