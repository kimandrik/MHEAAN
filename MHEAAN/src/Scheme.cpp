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

Scheme::Scheme(SecretKey* secretKey, Ring* ring, bool isSerialized) : ring(ring), isSerialized(isSerialized) {
	addEncKey(secretKey);
	addMultKey(secretKey);
};


//----------------------------------------------------------------------------------
//   KEYS GENERATION
//----------------------------------------------------------------------------------


void Scheme::addEncKey(SecretKey* secretKey) {

	long np = ceil((1 + ring->logQQ + ring->logN + 3)/(double)ring->pbnd);
	ZZ* ax = ring->sampleUniform(ring->logQQ);
	ZZ* bx = ring->mult(secretKey->sx, ax, np, ring->QQ);
	ZZ* ex = ring->sampleGauss();
	ring->subAndEqual2(ex, bx, ring->QQ);
	delete[] ex;

	np = ceil((2 * ring->logQQ + ring->logN + 3)/(double)ring->pbnd);
	uint64_t* rax = ring->toNTT(ax, np);
	uint64_t* rbx = ring->toNTT(bx, np);
	delete[] ax;
	delete[] bx;

	Key* key = new Key(rax, rbx, ring->N, np);
	if(isSerialized) {
		string path = "serkey/ENCRYPTION.txt";
		SerializationUtils::writeKey(key, path);
		serKeyMap.insert(pair<long, string>(ENCRYPTION, path));
		delete key;
	} else {
		keyMap.insert(pair<long, Key*>(ENCRYPTION, key));
	}
}

void Scheme::addMultKey(SecretKey* secretKey) {
	long np = ceil((1 + 1 + ring->logN + 3)/(double)ring->pbnd);
	ZZ* sx2 = ring->square(secretKey->sx, np, ring->Q);
	ring->leftShiftAndEqual(sx2, ring->logQ, ring->QQ);
	ZZ* ex = ring->sampleGauss();

	ring->addAndEqual(ex, sx2, ring->QQ);
	delete[] sx2;

	np = ceil((1 + ring->logQQ + ring->logN + 3)/(double)ring->pbnd);
	ZZ* ax = ring->sampleUniform(ring->logQQ);
	ZZ* bx = ring->mult(secretKey->sx, ax, np, ring->QQ);
	ring->subAndEqual2(ex, bx, ring->QQ);
	delete[] ex;

	np = ceil((2 * ring->logQQ + ring->logN + 3)/(double)ring->pbnd);
	uint64_t* rax = ring->toNTT(ax, np);
	uint64_t* rbx = ring->toNTT(bx, np);
	delete[] ax; delete[] bx;

	Key* key = new Key(rax, rbx, ring->N, np);
	if(isSerialized) {
		string path = "serkey/MULTIPLICATION.txt";
		SerializationUtils::writeKey(key, path);
		serKeyMap.insert(pair<long, string>(MULTIPLICATION, path));
		delete key;
	} else {
		keyMap.insert(pair<long, Key*>(MULTIPLICATION, key));
	}
}

void Scheme::addConjKey(SecretKey* secretKey) {
	ZZ* sxcnj = ring->conjugate(secretKey->sx);
	ring->leftShiftAndEqual(sxcnj, ring->logQ, ring->QQ);
	ZZ* ex = ring->sampleGauss();
	ring->addAndEqual(ex, sxcnj, ring->QQ);
	delete[] sxcnj;

	ZZ* ax = ring->sampleUniform(ring->logQQ);
	long np = ceil((1 + ring->logQQ + ring->logN + 3)/(double)ring->pbnd);
	ZZ* bx = ring->mult(secretKey->sx, ax, np, ring->QQ);
	ring->subAndEqual2(ex, bx, ring->QQ);
	delete[] ex;

	np = ceil((2 * ring->logQQ + ring->logN + 3)/(double)ring->pbnd);
	uint64_t* rax = ring->toNTT(ax, np);
	uint64_t* rbx = ring->toNTT(bx, np);
	delete[] ax;
	delete[] bx;

	Key* key = new Key(rax, rbx, ring->N, np);
	if(isSerialized) {
		string path = "serkey/CONJUGATION.txt";
		SerializationUtils::writeKey(key, path);
		serKeyMap.insert(pair<long, string>(CONJUGATION, path));
		delete key;
	} else {
		keyMap.insert(pair<long, Key*>(CONJUGATION, key));
	}
}

void Scheme::addLeftRotKey(SecretKey* secretKey, long r0, long r1) {
	ZZ* sxrot = ring->leftRotate(secretKey->sx, r0, r1);
	ring->leftShiftAndEqual(sxrot, ring->logQ, ring->QQ);
	ZZ* ex = ring->sampleGauss();
	ring->addAndEqual(ex, sxrot, ring->QQ);
	delete[] sxrot;

	long np = ceil((1 + ring->logQQ + ring->logN + 3)/(double)ring->pbnd);
	ZZ* ax = ring->sampleUniform(ring->logQQ);
	ZZ* bx = ring->mult(secretKey->sx, ax, np, ring->QQ);
	ring->subAndEqual2(ex, bx, ring->QQ);
	delete[] ex;

	np = ceil((2 * ring->logQQ + ring->logN + 3)/(double)ring->pbnd);
	uint64_t* rax = ring->toNTT(ax, np);
	uint64_t* rbx = ring->toNTT(bx, np);
	delete[] ax;
	delete[] bx;

	Key* key = new Key(rax, rbx, ring->N, np);
	if(isSerialized) {
		string path = "serkey/ROTATION_" + to_string(r0) + "_" + to_string(r1) + ".txt";
		SerializationUtils::writeKey(key, path);
		serLeftRotKeyMap.insert(pair<pair<long, long>, string>({r0, r1}, path));
		delete key;
	} else {
		leftRotKeyMap.insert(pair<pair<long, long>, Key*>({r0, r1}, key));
	}
}

void Scheme::addLeftX0RotKeys(SecretKey* secretKey) {
	for (long i = 1; i < ring->N0h; i <<= 1) {
		if(leftRotKeyMap.find({i, 0}) == leftRotKeyMap.end() && serLeftRotKeyMap.find({i, 0}) == serLeftRotKeyMap.end()) {
			addLeftRotKey(secretKey, i, 0);
		}
	}
}

void Scheme::addLeftX1RotKeys(SecretKey* secretKey) {
	for (long i = 1; i < ring->N1; i <<=1) {
		if(leftRotKeyMap.find({0, i}) == leftRotKeyMap.end() && serLeftRotKeyMap.find({0, i}) == serLeftRotKeyMap.end()) {
			addLeftRotKey(secretKey, 0, i);
		}
	}
}

void Scheme::addRightX0RotKeys(SecretKey* secretKey) {
	for (long i = 1; i < ring->N0h; i <<=1) {
		long idx = ring->N0h - i;
		if(leftRotKeyMap.find({idx, 0}) == leftRotKeyMap.end() && serLeftRotKeyMap.find({idx, 0}) == serLeftRotKeyMap.end()) {
			addLeftRotKey(secretKey, idx, 0);
		}
	}
}

void Scheme::addRightX1RotKeys(SecretKey* secretKey) {
	for (long i = 1; i < ring->N1; i<<=1) {
		long idx = ring->N1 - i;
		if(leftRotKeyMap.find({0, idx}) == leftRotKeyMap.end() && serLeftRotKeyMap.find({0, idx}) == serLeftRotKeyMap.end()) {
			addLeftRotKey(secretKey, 0, idx);
		}
	}
}

void Scheme::addBootKey(SecretKey* secretKey, long logn0, long logn1, long logp) {
	ring->addBootContext(logn0, logn1, logp);

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

void Scheme::addSqrMatKeys(SecretKey* secretKey, long logn, long logp) {
	ring->addSqrMatContext(logn, logp);
	long n = (1 << logn);
	for (long i = 1; i < n; ++i) {
		long idx = ring->N0h - i;
		if(leftRotKeyMap.find({idx, 0}) == leftRotKeyMap.end() && serLeftRotKeyMap.find({idx, 0}) == serLeftRotKeyMap.end()) {
			addLeftRotKey(secretKey, idx, 0);
		}
	}
	addLeftX1RotKeys(secretKey);
}


//----------------------------------------------------------------------------------
//   ENCODING & DECODING
//----------------------------------------------------------------------------------


Plaintext* Scheme::encode(complex<double>* vals, long n0, long n1, long logp, long logq) {
	ZZ* mx = ring->encode(vals, n0, n1, logp + ring->logQ);
	return new Plaintext(mx, logp, logq, ring->N0, ring->N1, n0, n1);
}

Plaintext* Scheme::encode(double* vals, long n0, long n1, long logp, long logq) {
	ZZ* mx = ring->encode(vals, n0, n1, logp + ring->logQ);
	return new Plaintext(mx, logp, logq, ring->N0, ring->N1, n0, n1);
}

Plaintext* Scheme::encodeSingle(complex<double> val, long logp, long logq) {
	ZZ* mx = new ZZ[ring->N];

	mx[0] = EvaluatorUtils::scaleUpToZZ(val.real(), logp + ring->logQ);
	mx[ring->Nh] = EvaluatorUtils::scaleUpToZZ(val.imag(), logp + ring->logQ);

	return new Plaintext(mx, logp, logq, ring->N0, ring->N1, 1, 1);
}

Plaintext* Scheme::encodeSingle(double val, long logp, long logq) {
	ZZ* mx = new ZZ[ring->N];
	mx[0] = EvaluatorUtils::scaleUpToZZ(val, logp + ring->logQ);
	return new Plaintext(mx, logp, logq, ring->N0, ring->N1, 1, 1);
}

complex<double>* Scheme::decode(Plaintext* msg) {
	complex<double>* vals = ring->decode(msg->mx, msg->n0, msg->n1, msg->logp, msg->logq);
	return vals;
}

complex<double> Scheme::decodeSingle(Plaintext* msg) {
	complex<double> res;
	ZZ q = ring->qvec[msg->logq];
	ZZ qh = ring->qvec[msg->logq - 1];
	ZZ tmp = msg->mx[0];
	while(tmp < 0) tmp += q;
	while(tmp > qh) tmp -= q;
	res.real(EvaluatorUtils::scaleDownToReal(tmp, msg->logp));

	tmp = msg->mx[ring->Nh];
	while(tmp < 0) tmp += q;
	while(tmp > qh) tmp -= q;
	res.imag(EvaluatorUtils::scaleDownToReal(tmp, msg->logp));

	return res;
}


//----------------------------------------------------------------------------------
//   ENCRYPTION & DECRYPTION
//----------------------------------------------------------------------------------


Ciphertext* Scheme::encryptMsg(Plaintext* msg) {
	ZZ qQ = ring->qvec[msg->logq + ring->logQ];

	Key* key = isSerialized ? SerializationUtils::readKey(serKeyMap.at(ENCRYPTION)) : keyMap.at(ENCRYPTION);
	long np = ceil((1 + ring->logQQ + ring->logN + 3)/(double)ring->pbnd);

	ZZ* vx = ring->sampleZO();
	ZZ* ax = ring->multNTT(vx, key->rax, np, qQ);
	ZZ* ex = ring->sampleGauss();
	ring->addAndEqual(ax, ex, qQ);

	ZZ* bx = ring->multNTT(vx, key->rbx, np, qQ);
	ring->sampleGauss(ex);
	ring->addAndEqual(bx, ex, qQ);
	delete[] ex;
	delete[] vx;

	ring->addAndEqual(bx, msg->mx, qQ);
	ring->rightShiftAndEqual(ax, ring->logQ);
	ring->rightShiftAndEqual(bx, ring->logQ);

	if(isSerialized) delete key;

	return new Ciphertext(ax, bx, msg->logp, msg->logq, msg->N0, msg->N1, msg->n0, msg->n1);
}

Ciphertext* Scheme::encrypt(complex<double>* vals, long n0, long n1, long logp, long logq) {
	Plaintext* msg = encode(vals, n0, n1, logp, logq);
	Ciphertext* res = encryptMsg(msg);
	delete msg;
	return res;
}

Ciphertext* Scheme::encrypt(double* vals, long n0, long n1, long logp, long logq) {
	Plaintext* msg = encode(vals, n0, n1, logp, logq);
	Ciphertext* res = encryptMsg(msg);
	delete msg;
	return res;
}

Ciphertext* Scheme::encryptSingle(complex<double> val, long logp, long logq) {
	Plaintext* msg = encodeSingle(val, logp, logq);
	Ciphertext* res = encryptMsg(msg);
	delete msg;
	return res;
}

Ciphertext* Scheme::encryptSingle(double val, long logp, long logq) {
	Plaintext* msg = encodeSingle(val, logp, logq);
	Ciphertext* res = encryptMsg(msg);
	delete msg;
	return res;
}

Ciphertext* Scheme::encryptZeros(long n0, long n1, long logp, long logq) {
	Ciphertext* czeros = encryptSingle(0.0, logp, logq);
	czeros->n0 = n0;
	czeros->n1 = n1;
	return czeros;
}

Plaintext* Scheme::decryptMsg(SecretKey* secretKey, Ciphertext* cipher) {
	ZZ q = ring->qvec[cipher->logq];

	long np = ceil((1 + cipher->logq + ring->logN + 3)/(double)ring->pbnd);
	ZZ* mx = ring->mult(cipher->ax, secretKey->sx, np, q);
	ring->addAndEqual(mx, cipher->bx, q);

	return new Plaintext(mx, cipher->logp, cipher->logq, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}

complex<double>* Scheme::decrypt(SecretKey* secretKey, Ciphertext* cipher) {
	Plaintext* msg = decryptMsg(secretKey, cipher);
	complex<double>* res = decode(msg);
	delete msg;
	return res;
}

complex<double> Scheme::decryptSingle(SecretKey* secretKey, Ciphertext* cipher) {
	Plaintext* msg = decryptMsg(secretKey, cipher);
	complex<double> res = decodeSingle(msg);
	delete msg;
	return res;

}


//----------------------------------------------------------------------------------
//   HOMOMORPHIC OPERATIONS
//----------------------------------------------------------------------------------


Ciphertext* Scheme::negate(Ciphertext* cipher) {
	ZZ* ax = ring->negate(cipher->ax);
	ZZ* bx = ring->negate(cipher->bx);

	return new Ciphertext(ax, bx, cipher->logp, cipher->logq, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}

void Scheme::negateAndEqual(Ciphertext* cipher) {
	ring->negateAndEqual(cipher->ax);
	ring->negateAndEqual(cipher->bx);
}

Ciphertext* Scheme::add(Ciphertext* cipher1, Ciphertext* cipher2) {
	ZZ q = ring->qvec[cipher1->logq];

	ZZ* ax = ring->add(cipher1->ax, cipher2->ax, q);
	ZZ* bx = ring->add(cipher1->bx, cipher2->bx, q);

	return new Ciphertext(ax, bx, cipher1->logp, cipher1->logq, cipher1->N0, cipher1->N1, cipher1->n0, cipher1->n1);
}

void Scheme::addAndEqual(Ciphertext* cipher1, Ciphertext* cipher2) {
	ZZ q = ring->qvec[cipher1->logq];

	ring->addAndEqual(cipher1->ax, cipher2->ax, q);
	ring->addAndEqual(cipher1->bx, cipher2->bx, q);
}

Ciphertext* Scheme::addConst(Ciphertext* cipher, double cnst, long logp) {
	ZZ q = ring->qvec[cipher->logq];

	Ciphertext* res = new Ciphertext(cipher);
	ZZ cnstZZ = logp < 0 ? EvaluatorUtils::scaleUpToZZ(cnst, cipher->logp) : EvaluatorUtils::scaleUpToZZ(cnst, logp);
	for (long i = 0; i < ring->N1; ++i) {
		AddMod(res->bx[i*ring->N0], res->bx[i*ring->N0], -cnstZZ, q);

	}
	return res;
}

Ciphertext* Scheme::addConst(Ciphertext* cipher, RR& cnst, long logp) {
	ZZ q = ring->qvec[cipher->logq];

	Ciphertext* res = new Ciphertext(cipher);
	ZZ cnstZZ = logp < 0 ? -EvaluatorUtils::scaleUpToZZ(cnst, cipher->logp) : -EvaluatorUtils::scaleUpToZZ(cnst, logp);
	for (long i = 0; i < ring->N; i+=ring->N0) {
		AddMod(res->bx[i], res->bx[i], cnstZZ, q);
	}
	return res;
}

void Scheme::addConstAndEqual(Ciphertext* cipher, double cnst, long logp) {
	ZZ q = ring->qvec[cipher->logq];

	ZZ cnstZZ = logp < 0 ? -EvaluatorUtils::scaleUpToZZ(cnst, cipher->logp) : -EvaluatorUtils::scaleUpToZZ(cnst, logp);
	for (long i = 0; i < ring->N; i+=ring->N0) {
		AddMod(cipher->bx[i], cipher->bx[i], cnstZZ, q);
	}
}

void Scheme::addConstAndEqual(Ciphertext* cipher, RR& cnst, long logp) {
	ZZ q = ring->qvec[cipher->logq];

	ZZ cnstZZ = logp < 0 ? -EvaluatorUtils::scaleUpToZZ(cnst, cipher->logp) : -EvaluatorUtils::scaleUpToZZ(cnst, logp);
	for (long i = 0; i < ring->N; i+=ring->N0) {
		AddMod(cipher->bx[i], cipher->bx[i], cnstZZ, q);
	}
}

Ciphertext* Scheme::addPoly(Ciphertext* cipher, ZZ* poly, long logp) {
	ZZ q = ring->qvec[cipher->logq];

	Ciphertext* res = new Ciphertext(cipher);
	ring->addAndEqual(res->bx, poly, q);
	return res;
}

void Scheme::addPolyAndEqual(Ciphertext* cipher, ZZ* poly, long logp) {
	ZZ q = ring->qvec[cipher->logq];
	ring->addAndEqual(cipher->bx, poly, q);
}

Ciphertext* Scheme::sub(Ciphertext* cipher1, Ciphertext* cipher2) {
	ZZ q = ring->qvec[cipher1->logq];

	ZZ* ax = ring->sub(cipher1->ax, cipher2->ax, q);
	ZZ* bx = ring->sub(cipher1->bx, cipher2->bx, q);

	return new Ciphertext(ax, bx, cipher1->logp, cipher1->logq, cipher1->N0, cipher1->N1, cipher1->n0, cipher1->n1);
}

void Scheme::subAndEqual(Ciphertext* cipher1, Ciphertext* cipher2) {
	ZZ q = ring->qvec[cipher1->logq];

	ring->subAndEqual(cipher1->ax, cipher2->ax, q);
	ring->subAndEqual(cipher1->bx, cipher2->bx, q);
}

void Scheme::subAndEqual2(Ciphertext* cipher1, Ciphertext* cipher2) {
	ZZ q = ring->qvec[cipher1->logq];

	ring->subAndEqual2(cipher1->ax, cipher2->ax, q);
	ring->subAndEqual2(cipher1->bx, cipher2->bx, q);
}

Ciphertext* Scheme::imult(Ciphertext* cipher) {
	ZZ q = ring->qvec[cipher->logq];

	ZZ* ax = ring->multByMonomial(cipher->ax, ring->N0h, 0, q);
	ZZ* bx = ring->multByMonomial(cipher->bx, ring->N0h, 0, q);

	return new Ciphertext(ax, bx, cipher->logp, cipher->logq, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}

Ciphertext* Scheme::idiv(Ciphertext* cipher) {
	ZZ q = ring->qvec[cipher->logq];

	ZZ* ax = ring->multByMonomial(cipher->ax, 3 * ring->N0h, 0, q);
	ZZ* bx = ring->multByMonomial(cipher->bx, 3 * ring->N0h, 0, q);

	return new Ciphertext(ax, bx, cipher->logp, cipher->logq, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}

void Scheme::imultAndEqual(Ciphertext* cipher) {
	ZZ q = ring->qvec[cipher->logq];

	ring->multByMonomialAndEqual(cipher->ax, ring->N0h, 0, q);
	ring->multByMonomialAndEqual(cipher->bx, ring->N0h, 0, q);
}

void Scheme::idivAndEqual(Ciphertext* cipher) {
	ZZ q = ring->qvec[cipher->logq];

	ring->multByMonomialAndEqual(cipher->ax, 3 * ring->N0h, 0, q);
	ring->multByMonomialAndEqual(cipher->bx, 3 * ring->N0h, 0, q);
}

Ciphertext* Scheme::mult(Ciphertext* cipher1, Ciphertext* cipher2) {
	ZZ q = ring->qvec[cipher1->logq];
	ZZ qQ = ring->qvec[cipher1->logq + ring->logQ];

	long np = ceil((2 + cipher1->logq + cipher2->logq + ring->logN + 3)/(double)ring->pbnd);
	uint64_t* ra1 = ring->toNTT(cipher1->ax, np);
	uint64_t* rb1 = ring->toNTT(cipher1->bx, np);
	uint64_t* ra2 = ring->toNTT(cipher2->ax, np);
	uint64_t* rb2 = ring->toNTT(cipher2->bx, np);

	ZZ* aax = ring->multDNTT(ra1, ra2, np, q);
	ZZ* bbx = ring->multDNTT(rb1, rb2, np, q);

	ring->addNTTAndEqual(ra1, rb1, np);
	ring->addNTTAndEqual(ra2, rb2, np);
	ZZ* abx = ring->multDNTT(ra1, ra2, np, q);
	delete[] ra1; delete[] ra2; delete[] rb1; delete[] rb2;

	np = ceil((cipher1->logq + ring->logQQ + ring->logN + 3)/(double)ring->pbnd);
	uint64_t* raa = ring->toNTT(aax, np);
	Key* key = isSerialized ? SerializationUtils::readKey(serKeyMap.at(MULTIPLICATION)) : keyMap.at(MULTIPLICATION);
	ZZ* ax = ring->multDNTT(raa, key->rax, np, qQ);
	ZZ* bx = ring->multDNTT(raa, key->rbx, np, qQ);
	delete[] raa;
	if(isSerialized) delete key;

	ring->rightShiftAndEqual(ax, ring->logQ);
	ring->rightShiftAndEqual(bx, ring->logQ);

	ring->addAndEqual(ax, abx, q);
	ring->subAndEqual(ax, bbx, q);
	ring->subAndEqual(ax, aax, q);
	ring->addAndEqual(bx, bbx, q);
	delete[] abx; delete[] aax; delete[] bbx;

	return new Ciphertext(ax, bx, cipher1->logp + cipher2->logp, cipher1->logq, cipher1->N0, cipher1->N1, cipher1->n0, cipher1->n1);
}

void Scheme::multAndEqual(Ciphertext* cipher1, Ciphertext* cipher2) {
	ZZ q = ring->qvec[cipher1->logq];
	ZZ qQ = ring->qvec[cipher1->logq + ring->logQ];

	long np = ceil((2 + cipher1->logq + cipher2->logq + ring->logN + 3)/(double)ring->pbnd);
	uint64_t* ra1 = ring->toNTT(cipher1->ax, np);
	uint64_t* rb1 = ring->toNTT(cipher1->bx, np);
	uint64_t* ra2 = ring->toNTT(cipher2->ax, np);
	uint64_t* rb2 = ring->toNTT(cipher2->bx, np);

	ZZ* aax = ring->multDNTT(ra1, ra2, np, q);
	ZZ* bbx = ring->multDNTT(rb1, rb2, np, q);

	np = ceil((cipher1->logq + ring->logQQ + ring->logN + 3)/(double)ring->pbnd);
	uint64_t* raa = ring->toNTT(aax, np);
	Key* key = isSerialized ? SerializationUtils::readKey(serKeyMap.at(MULTIPLICATION)) : keyMap.at(MULTIPLICATION);
	ring->multDNTT(cipher1->ax, raa, key->rax, np, qQ);
	ring->multDNTT(cipher1->bx, raa, key->rbx, np, qQ);
	delete[] raa;
	if(isSerialized) delete key;

	ring->rightShiftAndEqual(cipher1->ax, ring->logQ);
	ring->rightShiftAndEqual(cipher1->bx, ring->logQ);

	ring->addNTTAndEqual(ra1, rb1, np);
	ring->addNTTAndEqual(ra2, rb2, np);
	delete[] rb1; delete[] rb2;

	ZZ* abx = ring->multDNTT(ra1, ra2, np, q);
	delete[] ra1; delete[] ra2;

	ring->addAndEqual(cipher1->ax, abx, q);
	ring->subAndEqual(cipher1->ax, bbx, q);
	ring->subAndEqual(cipher1->ax, aax, q);
	ring->addAndEqual(cipher1->bx, bbx, q);
	delete[] abx; delete[] aax; delete[] bbx;

	cipher1->logp += cipher2->logp;
}

Ciphertext* Scheme::square(Ciphertext* cipher) {
	ZZ q = ring->qvec[cipher->logq];
	ZZ qQ = ring->qvec[cipher->logq + ring->logQ];

	long np = ceil((2 * cipher->logq + ring->logN + 3)/(double)ring->pbnd);
	uint64_t* ra = ring->toNTT(cipher->ax, np);
	uint64_t* rb = ring->toNTT(cipher->bx, np);

	ZZ* bbx = ring->squareNTT(rb, np, q);
	ZZ* aax = ring->squareNTT(ra, np, q);
	ZZ* abx = ring->multDNTT(ra, rb, np, q);
	ring->leftShiftAndEqual(abx, 1, q);
	delete[] ra; delete[] rb;

	np = ceil((cipher->logq + ring->logQQ + ring->logN + 3)/(double)ring->pbnd);
	uint64_t* raa = ring->toNTT(aax, np);
	delete[] aax;
	Key* key = isSerialized ? SerializationUtils::readKey(serKeyMap.at(MULTIPLICATION)) : keyMap.at(MULTIPLICATION);
	ZZ* ax = ring->multDNTT(raa, key->rax, np, qQ);
	ZZ* bx = ring->multDNTT(raa, key->rbx, np, qQ);
	delete[] raa;
	if(isSerialized) delete key;

	ring->rightShiftAndEqual(ax, ring->logQ);
	ring->rightShiftAndEqual(bx, ring->logQ);

	ring->addAndEqual(ax, abx, q);
	ring->addAndEqual(bx, bbx, q);
	delete[] abx; delete[] bbx;

	return new Ciphertext(ax, bx, 2 * cipher->logp, cipher->logq, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}

void Scheme::squareAndEqual(Ciphertext* cipher) {
	ZZ q = ring->qvec[cipher->logq];
	ZZ qQ = ring->qvec[cipher->logq + ring->logQ];

	long np = ceil((2 * cipher->logq + ring->logN + 3)/(double)ring->pbnd);
	uint64_t* ra = ring->toNTT(cipher->ax, np);
	uint64_t* rb = ring->toNTT(cipher->bx, np);

	ZZ* bbx = ring->squareNTT(rb, np, q);
	ZZ* abx = ring->multDNTT(ra, rb, np, q);
	ring->leftShiftAndEqual(abx, 1, q);
	ZZ* aax = ring->squareNTT(ra, np, q);
	delete[] ra; delete[] rb;

	np = ceil((cipher->logq + ring->logQQ + ring->logN + 3)/(double)ring->pbnd);
	uint64_t* raa = ring->toNTT(aax, np);
	delete[] aax;
	Key* key = isSerialized ? SerializationUtils::readKey(serKeyMap.at(MULTIPLICATION)) : keyMap.at(MULTIPLICATION);
	ring->multDNTT(cipher->ax, raa, key->rax, np, qQ);
	ring->multDNTT(cipher->bx, raa, key->rbx, np, qQ);
	delete[] raa;
	if(isSerialized) delete key;

	ring->rightShiftAndEqual(cipher->ax, ring->logQ);
	ring->rightShiftAndEqual(cipher->bx, ring->logQ);

	ring->addAndEqual(cipher->ax, abx, q);
	ring->addAndEqual(cipher->bx, bbx, q);
	delete[] abx; delete[] bbx;

	cipher->logp *= 2;
}

Ciphertext* Scheme::multConst(Ciphertext* cipher, RR& cnst, long logp) {
	ZZ q = ring->qvec[cipher->logq];

	ZZ cnstZZ = EvaluatorUtils::scaleUpToZZ(cnst, logp);

	ZZ* ax = ring->multByConst(cipher->ax, cnstZZ, q);
	ZZ* bx = ring->multByConst(cipher->bx, cnstZZ, q);

	return new Ciphertext(ax, bx, cipher->logp + logp, cipher->logq, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}


Ciphertext* Scheme::multConst(Ciphertext* cipher, double cnst, long logp) {
	ZZ q = ring->qvec[cipher->logq];

	ZZ cnstZZ = EvaluatorUtils::scaleUpToZZ(cnst, logp);

	ZZ* ax = ring->multByConst(cipher->ax, cnstZZ, q);
	ZZ* bx = ring->multByConst(cipher->bx, cnstZZ, q);

	return new Ciphertext(ax, bx, cipher->logp + logp, cipher->logq, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}

Ciphertext* Scheme::multConst(Ciphertext* cipher, complex<double> cnst, long logp) {
	ZZ q = ring->qvec[cipher->logq];

	ZZ* axi = ring->multByMonomial(cipher->ax, ring->N0h, 0, q);
	ZZ* bxi = ring->multByMonomial(cipher->bx, ring->N0h, 0, q);

	ZZ cnstrZZ = EvaluatorUtils::scaleUpToZZ(cnst.real(), logp);
	ZZ cnstiZZ = EvaluatorUtils::scaleUpToZZ(cnst.imag(), logp);

	ZZ* ax = ring->multByConst(cipher->ax, cnstrZZ, q);
	ZZ* bx = ring->multByConst(cipher->bx, cnstrZZ, q);

	ring->multByConstAndEqual(axi, cnstiZZ, q);
	ring->multByConstAndEqual(bxi, cnstiZZ, q);

	ring->addAndEqual(ax, axi, q);
	ring->addAndEqual(bx, bxi, q);
	delete[] axi; delete bxi;

	return new Ciphertext(ax, bx, cipher->logp + logp, cipher->logq, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}

void Scheme::multConstAndEqual(Ciphertext* cipher, RR& cnst, long logp) {
	ZZ q = ring->qvec[cipher->logq];

	ZZ cnstZZ = EvaluatorUtils::scaleUpToZZ(cnst, logp);

	ring->multByConstAndEqual(cipher->ax, cnstZZ, q);
	ring->multByConstAndEqual(cipher->bx, cnstZZ, q);

	cipher->logp += logp;
}

void Scheme::multConstAndEqual(Ciphertext* cipher, double cnst, long logp) {
	ZZ q = ring->qvec[cipher->logq];

	ZZ cnstZZ = EvaluatorUtils::scaleUpToZZ(cnst, logp);

	ring->multByConstAndEqual(cipher->ax, cnstZZ, q);
	ring->multByConstAndEqual(cipher->bx, cnstZZ, q);

	cipher->logp += logp;
}

void Scheme::multConstAndEqual(Ciphertext* cipher, complex<double> cnst, long logp) {
	ZZ q = ring->qvec[cipher->logq];


	ZZ cnstrZZ = EvaluatorUtils::scaleUpToZZ(cnst.real(), logp);
	ZZ cnstiZZ = EvaluatorUtils::scaleUpToZZ(cnst.imag(), logp);

	ZZ* axi = ring->multByMonomial(cipher->ax, ring->N0h, 0, q);
	ZZ* bxi = ring->multByMonomial(cipher->bx, ring->N0h, 0, q);

	ring->multByConstAndEqual(axi, cnstiZZ, q);
	ring->multByConstAndEqual(bxi, cnstiZZ, q);

	ring->multByConstAndEqual(cipher->ax, cnstrZZ, q);
	ring->multByConstAndEqual(cipher->bx, cnstrZZ, q);

	ring->addAndEqual(cipher->ax, axi, q);
	ring->addAndEqual(cipher->bx, bxi, q);
	delete[] axi; delete[] bxi;
	cipher->logp += logp;
}


Ciphertext* Scheme::multPolyX0(Ciphertext* cipher, ZZ* poly, long logp) {
	ZZ q = ring->qvec[cipher->logq];

	long bnd = ring->MaxBits(poly, ring->N0);
	long np = ceil((cipher->logq + bnd + ring->logN0 + 3)/(double)ring->pbnd);
	uint64_t* rpoly = ring->toNTTX0(poly, np);
	ZZ* ax = ring->multNTTX0(cipher->ax, rpoly, np, q);
	ZZ* bx = ring->multNTTX0(cipher->bx, rpoly, np, q);
	delete[] rpoly;

	return new Ciphertext(ax, bx, cipher->logp + logp, cipher->logq, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}

void Scheme::multPolyX0AndEqual(Ciphertext* cipher, ZZ* poly, long logp) {
	ZZ q = ring->qvec[cipher->logq];

	long bnd = ring->MaxBits(poly, ring->N0);
	long np = ceil((cipher->logq + bnd + ring->logN0 + 3)/(double)ring->pbnd);
	uint64_t* rpoly = ring->toNTTX0(poly, np);
	ring->multNTTX0AndEqual(cipher->ax, rpoly, np, q);
	ring->multNTTX0AndEqual(cipher->bx, rpoly, np, q);
	delete[] rpoly;
	cipher->logp += logp;
}

Ciphertext* Scheme::multPolyNTTX0(Ciphertext* cipher, uint64_t* rpoly, long bnd, long logp) {
	ZZ q = ring->qvec[cipher->logq];

	long np = ceil((cipher->logq + bnd + ring->logN0 + 3)/(double)ring->pbnd);
	ZZ* ax = ring->multNTTX0(cipher->ax, rpoly, np, q);
	ZZ* bx = ring->multNTTX0(cipher->bx, rpoly, np, q);

	return new Ciphertext(ax, bx, cipher->logp + logp, cipher->logq, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}

void Scheme::multPolyNTTX0AndEqual(Ciphertext* cipher, uint64_t* rpoly, long bnd, long logp) {
	ZZ q = ring->qvec[cipher->logq];

	long np = ceil((cipher->logq + bnd + ring->logN0 + 3)/(double)ring->pbnd);
	ring->multNTTX0AndEqual(cipher->ax, rpoly, np, q);
	ring->multNTTX0AndEqual(cipher->bx, rpoly, np, q);
	cipher->logp += logp;
}

Ciphertext* Scheme::multPolyX1(Ciphertext* cipher, ZZ* rpoly, ZZ* ipoly, long logp) {
	ZZ q = ring->qvec[cipher->logq];

	ZZ* axi = ring->multByMonomial(cipher->ax, ring->N0h, 0, q);
	ZZ* bxi = ring->multByMonomial(cipher->bx, ring->N0h, 0, q);

	long bnd = ring->MaxBits(ipoly, ring->N1);
	long np = ceil((cipher->logq + bnd + ring->logN1 + 3)/(double)ring->pbnd);
	uint64_t* ripoly = ring->toNTT(ipoly, np);
	ring->multNTTX1AndEqual(axi, ripoly, np, q);
	ring->multNTTX1AndEqual(bxi, ripoly, np, q);
	delete[] ripoly;

	bnd = ring->MaxBits(rpoly, ring->N1);
	np = ceil((cipher->logq + bnd + ring->logN1 + 3)/(double)ring->pbnd);
	uint64_t* rrpoly = ring->toNTT(rpoly, np);
	ZZ* axr = ring->multNTTX1(cipher->ax, rrpoly, np, q);
	ZZ* bxr = ring->multNTTX1(cipher->bx, rrpoly, np, q);
	delete[] rrpoly;

	ring->addAndEqual(axr, axi, q);
	ring->addAndEqual(bxr, bxi, q);
	delete[] axi; delete[] bxi;

	return new Ciphertext(axr, bxr, cipher->logp + logp, cipher->logq, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}

void Scheme::multPolyX1AndEqual(Ciphertext* cipher, ZZ* rpoly, ZZ* ipoly, long logp) {
	ZZ q = ring->qvec[cipher->logq];

	ZZ* axi = ring->multByMonomial(cipher->ax, ring->N0h, 0, q);
	ZZ* bxi = ring->multByMonomial(cipher->bx, ring->N0h, 0, q);

	long bnd = ring->MaxBits(rpoly, ring->N1);
	long np = ceil((cipher->logq + bnd + ring->logN1 + 3)/(double)ring->pbnd);
	uint64_t* rrpoly = ring->toNTT(rpoly, np);
	ring->multNTTX1AndEqual(cipher->ax, rrpoly, np, q);
	ring->multNTTX1AndEqual(cipher->bx, rrpoly, np, q);
	delete[] rrpoly;

	bnd = ring->MaxBits(ipoly, ring->N1);
	np = ceil((cipher->logq + bnd + ring->logN1 + 3)/(double)ring->pbnd);
	uint64_t* ripoly = ring->toNTT(ipoly, np);
	ring->multNTTX1AndEqual(axi, ripoly, np, q);
	ring->multNTTX1AndEqual(bxi, ripoly, np, q);
	delete[] ripoly;

	ring->addAndEqual(cipher->ax, axi, q);
	ring->addAndEqual(cipher->bx, bxi, q);
	delete[] axi; delete[] bxi;

	cipher->logp += logp;
}

Ciphertext* Scheme::multPoly(Ciphertext* cipher, ZZ* poly, long logp) {
	ZZ q = ring->qvec[cipher->logq];

	long bnd = ring->MaxBits(poly, ring->N);
	long np = ceil((cipher->logq + logp + ring->logN + 3)/(double)ring->pbnd);
	uint64_t* rpoly = ring->toNTT(poly, np);
	ZZ* ax = ring->multNTT(cipher->ax, rpoly, np, q);
	ZZ* bx = ring->multNTT(cipher->bx, rpoly, np, q);
	delete[] rpoly;

	return new Ciphertext(ax, bx, cipher->logp + logp, cipher->logq, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}

void Scheme::multPolyAndEqual(Ciphertext* cipher, ZZ* poly, long logp) {
	ZZ q = ring->qvec[cipher->logq];

	long bnd = ring->MaxBits(poly, ring->N);
	long np = ceil((cipher->logq + bnd + ring->logN + 3)/(double)ring->pbnd);
	uint64_t* rpoly = ring->toNTT(poly, np);
	ring->multNTTAndEqual(cipher->ax, rpoly, np, q);
	ring->multNTTAndEqual(cipher->bx, rpoly, np, q);
	delete[] rpoly;
	cipher->logp += logp;
}


Ciphertext Scheme::multPolyNTT(Ciphertext* cipher, uint64_t* rpoly, long bnd, long logp) {
	ZZ q = ring->qvec[cipher->logq];

	long np = ceil((cipher->logq + bnd + ring->logN + 3)/(double)ring->pbnd);
	ZZ* ax = ring->multNTT(cipher->ax, rpoly, np, q);
	ZZ* bx = ring->multNTT(cipher->bx, rpoly, np, q);

	return new Ciphertext(ax, bx, cipher->logp + logp, cipher->logq, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}

void Scheme::multPolyNTTAndEqual(Ciphertext* cipher, uint64_t* rpoly, long bnd, long logp) {
	ZZ q = ring->qvec[cipher->logq];

	long np = ceil((cipher->logq + bnd + ring->logN + 3)/(double)ring->pbnd);
	ring->multNTTAndEqual(cipher->ax, rpoly, np, q);
	ring->multNTTAndEqual(cipher->bx, rpoly, np, q);
	cipher->logp += logp;
}

Ciphertext* Scheme::multByMonomial(Ciphertext* cipher, const long d0, const long d1) {
	ZZ q = ring->qvec[cipher->logq];

	ZZ* ax = ring->multByMonomial(cipher->ax, d0, d1, q);
	ZZ* bx = ring->multByMonomial(cipher->bx, d0, d1, q);

	return new Ciphertext(ax, bx, cipher->logp, cipher->logq, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}

void Scheme::multByMonomialAndEqual(Ciphertext* cipher, const long d0, const long d1) {
	ZZ q = ring->qvec[cipher->logq];

	ring->multByMonomialAndEqual(cipher->ax, d0, d1, q);
	ring->multByMonomialAndEqual(cipher->bx, d0, d1, q);
}

Ciphertext* Scheme::multPo2(Ciphertext* cipher, long bits) {
	ZZ q = ring->qvec[cipher->logq];

	ZZ* ax = ring->leftShift(cipher->ax, bits, q);
	ZZ* bx = ring->leftShift(cipher->bx, bits, q);

	return new Ciphertext(ax, bx, cipher->logp, cipher->logq, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}

void Scheme::multPo2AndEqual(Ciphertext* cipher, long bits) {
	ZZ q = ring->qvec[cipher->logq];

	ring->leftShiftAndEqual(cipher->ax, bits, q);
	ring->leftShiftAndEqual(cipher->bx, bits, q);
}

Ciphertext* Scheme::divPo2(Ciphertext* cipher, long logd) {
	ZZ* ax = ring->rightShift(cipher->ax, logd);
	ZZ* bx = ring->rightShift(cipher->bx, logd);

	return new Ciphertext(ax, bx, cipher->logp, cipher->logq - logd, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}

void Scheme::divPo2AndEqual(Ciphertext* cipher, long logd) {
	ring->rightShiftAndEqual(cipher->ax, logd);
	ring->rightShiftAndEqual(cipher->bx, logd);
	cipher->logq -= logd;
}


//----------------------------------------------------------------------------------
//   RESCALING & MODULUS DOWN
//----------------------------------------------------------------------------------


Ciphertext* Scheme::reScaleBy(Ciphertext* cipher, long dlogq) {
	ZZ* ax = ring->rightShift(cipher->ax, dlogq);
	ZZ* bx = ring->rightShift(cipher->bx, dlogq);

	return new Ciphertext(ax, bx, cipher->logp - dlogq, cipher->logq - dlogq, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}

Ciphertext* Scheme::reScaleTo(Ciphertext* cipher, long logq) {
	long dlogq = cipher->logq - logq;
	ZZ* ax = ring->rightShift(cipher->ax, dlogq);
	ZZ* bx = ring->rightShift(cipher->bx, dlogq);

	return new Ciphertext(ax, bx, cipher->logp - dlogq, logq, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}

void Scheme::reScaleByAndEqual(Ciphertext* cipher, long dlogq) {
	ring->rightShiftAndEqual(cipher->ax, dlogq);
	ring->rightShiftAndEqual(cipher->bx, dlogq);
	cipher->logq -= dlogq;
	cipher->logp -= dlogq;
}

void Scheme::reScaleToAndEqual(Ciphertext* cipher, long logq) {
	long dlogq = cipher->logq - logq;
	ring->rightShiftAndEqual(cipher->ax, dlogq);
	ring->rightShiftAndEqual(cipher->bx, dlogq);
	cipher->logq = logq;
	cipher->logp -= dlogq;
}

Ciphertext* Scheme::modDownBy(Ciphertext* cipher, long dlogq) {
	long logq = cipher->logq - dlogq;
	ZZ q = ring->qvec[logq];

	ZZ* ax = ring->mod(cipher->ax, q);
	ZZ* bx = ring->mod(cipher->bx, q);

	return new Ciphertext(ax, bx, cipher->logp, logq, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}

void Scheme::modDownByAndEqual(Ciphertext* cipher, long dlogq) {
	ZZ q = ring->qvec[cipher->logq - dlogq];

	ring->modAndEqual(cipher->ax, q);
	ring->modAndEqual(cipher->bx, q);

	cipher->logq -= dlogq;
}

Ciphertext* Scheme::modDownTo(Ciphertext* cipher, long logq) {
	ZZ q = ring->qvec[logq];

	ZZ* ax = ring->mod(cipher->ax, q);
	ZZ* bx = ring->mod(cipher->bx, q);

	return new Ciphertext(ax, bx, cipher->logp, logq, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}

void Scheme::modDownToAndEqual(Ciphertext* cipher, long newlogq) {
	ZZ q = ring->qvec[newlogq];

	ring->modAndEqual(cipher->ax, q);
	ring->modAndEqual(cipher->bx, q);

	cipher->logq = newlogq;
}


//----------------------------------------------------------------------------------
//   ROTATIONS & CONJUGATIONS
//----------------------------------------------------------------------------------


Ciphertext* Scheme::leftRotateFast(Ciphertext* cipher, long r0, long r1) {
	ZZ q = ring->qvec[cipher->logq];
	ZZ qQ = ring->qvec[cipher->logq + ring->logQ];

	ZZ* bxrot = ring->leftRotate(cipher->bx, r0, r1);
	ZZ* axrot = ring->leftRotate(cipher->ax, r0, r1);

	long np = ceil((cipher->logq + ring->logQQ + ring->logN + 3)/(double)ring->pbnd);
	uint64_t* rarot = ring->toNTT(axrot, np);
	Key* key = isSerialized ? SerializationUtils::readKey(serLeftRotKeyMap.at({r0, r1})) : leftRotKeyMap.at({r0, r1});
	ZZ* ax = ring->multDNTT(rarot, key->rax, np, qQ);
	ZZ* bx = ring->multDNTT(rarot, key->rbx, np, qQ);
	delete[] rarot;
	if(isSerialized) delete key;

	ring->rightShiftAndEqual(ax, ring->logQ);
	ring->rightShiftAndEqual(bx, ring->logQ);

	ring->addAndEqual(bx, bxrot, q);
	delete[] bxrot; delete[] axrot;

	return new Ciphertext(ax, bx, cipher->logp, cipher->logq, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}

Ciphertext* Scheme::rightRotateFast(Ciphertext* cipher, long rx, long ry) {
	Ciphertext* res = new Ciphertext(cipher);
	rightRotateFastAndEqual(res, rx, ry);
	return res;
}

void Scheme::leftRotateFastAndEqual(Ciphertext* cipher, long r0, long r1) {
	ZZ q = ring->qvec[cipher->logq];
	ZZ qQ = ring->qvec[cipher->logq + ring->logQ];

	ZZ* axrot = ring->leftRotate(cipher->ax, r0, r1);
	ZZ* bxrot = ring->leftRotate(cipher->bx, r0, r1);

	Key* key = isSerialized ? SerializationUtils::readKey(serLeftRotKeyMap.at({r0, r1})) : leftRotKeyMap.at({r0, r1});

	long np = ceil((cipher->logq + ring->logQQ + ring->logN + 3)/(double)ring->pbnd);
	uint64_t* rarot = ring->toNTT(axrot, np);
	ring->multDNTT(cipher->bx, rarot, key->rbx, np, qQ);
	ring->multDNTT(cipher->ax, rarot, key->rax, np, qQ);
	delete[] rarot;
	if(isSerialized) delete key;

	ring->rightShiftAndEqual(cipher->ax, ring->logQ);
	ring->rightShiftAndEqual(cipher->bx, ring->logQ);

	ring->addAndEqual(cipher->bx, bxrot, q);
	delete[] axrot; delete[] bxrot;
}

void Scheme::rightRotateFastAndEqual(Ciphertext* cipher, long r0, long r1) {
	long rr0 = r0 == 0 ? 0 : ring->N0h - r0;
	long rr1 = r1 == 0 ? 0 : ring->N1 - r1;
	leftRotateFastAndEqual(cipher, rr0, rr1);
}

Ciphertext* Scheme::leftRotate(Ciphertext* cipher, long r0, long r1) {
	Ciphertext* res = new Ciphertext(cipher);
	leftRotateAndEqual(res, r0, r1);
	return res;
}

Ciphertext* Scheme::rightRotate(Ciphertext* cipher, long rx, long ry) {
	Ciphertext* res = new Ciphertext(cipher);
	rightRotateAndEqual(res, rx, ry);
	return res;
}


void Scheme::leftRotateAndEqual(Ciphertext* cipher, long r0, long r1) {
	r0 %= cipher->n0;
	long logr0 = log2((double)r0) + 1;
	long ipow;
	for (long i = 0; i < logr0; ++i) {
		if(bit(r0, i)) {
			ipow = 1 << i;
			leftRotateFastAndEqual(cipher, ipow, 0);
		}
	}
	r1 %= cipher->n1;
	long logr1 = log2((double)r1) + 1;
	for (long i = 0; i < logr1; ++i) {
		if(bit(r1, i)) {
			ipow = 1 << i;
			leftRotateFastAndEqual(cipher, 0, ipow);
		}
	}
}

void Scheme::rightRotateAndEqual(Ciphertext* cipher, long r0, long r1) {
	r0 %= cipher->n0;
	long logr0 = log2((double)r0) + 1;
	long ipow;
	for (long i = 0; i < logr0; ++i) {
		if(bit(r0, i)) {
			ipow = 1 << i;
			rightRotateFastAndEqual(cipher, ipow, 0);
		}
	}
	r1 %= cipher->n1;
	long logryRem = log2((double)r1) + 1;
	for (long i = 0; i < logryRem; ++i) {
		if(bit(r1, i)) {
			ipow = 1 << i;
			rightRotateFastAndEqual(cipher, 0, ipow);
		}
	}
}

Ciphertext* Scheme::conjugate(Ciphertext* cipher) {
	ZZ q = ring->qvec[cipher->logq];
	ZZ qQ = ring->qvec[cipher->logq + ring->logQ];

	ZZ* bxcnj = ring->conjugate(cipher->bx);
	ZZ* axcnj = ring->conjugate(cipher->ax);

	long np = ceil((cipher->logq + ring->logQQ + ring->logN + 3)/(double)ring->pbnd);
	uint64_t* racnj = ring->toNTT(axcnj, np);
	Key* key = isSerialized ? SerializationUtils::readKey(serKeyMap.at(CONJUGATION)) : keyMap.at(CONJUGATION);
	ZZ* ax = ring->multDNTT(racnj, key->rax, np, qQ);
	ZZ* bx = ring->multDNTT(racnj, key->rbx, np, qQ);
	delete[] racnj;
	if(isSerialized) delete key;

	ring->rightShiftAndEqual(ax, ring->logQ);
	ring->rightShiftAndEqual(bx, ring->logQ);

	ring->addAndEqual(bx, bxcnj, q);
	delete[] axcnj; delete[] bxcnj;

	return new Ciphertext(ax, bx, cipher->logp, cipher->logq, cipher->N0, cipher->N1, cipher->n0, cipher->n1);
}

void Scheme::conjugateAndEqual(Ciphertext* cipher) {
	ZZ q = ring->qvec[cipher->logq];
	ZZ qQ = ring->qvec[cipher->logq + ring->logQ];

	ZZ* axcnj = ring->conjugate(cipher->ax);
	ZZ* bxcnj = ring->conjugate(cipher->bx);

	Key* key = isSerialized ? SerializationUtils::readKey(serKeyMap.at(CONJUGATION)) : keyMap.at(CONJUGATION);

	long np = ceil((cipher->logq + ring->logQQ + ring->logN + 3)/(double)ring->pbnd);
	uint64_t* racnj = ring->toNTT(axcnj, np);
	ring->multDNTT(cipher->ax, racnj, key->rax, np, qQ);
	ring->multDNTT(cipher->bx, racnj, key->rbx, np, qQ);
	delete[] racnj;
	if(isSerialized) delete key;

	ring->rightShiftAndEqual(cipher->ax, ring->logQ);
	ring->rightShiftAndEqual(cipher->bx, ring->logQ);

	ring->addAndEqual(cipher->bx, bxcnj, q);
	delete[] axcnj; delete[] bxcnj;
}


//----------------------------------------------------------------------------------
//   BOOTSTRAPPING
//----------------------------------------------------------------------------------


void Scheme::normalizeAndEqual(Ciphertext* cipher) {
	ZZ q = ring->qvec[cipher->logq];
	ZZ qh = ring->qvec[cipher->logq - 1];
	for (long i = 0; i < ring->N; ++i) {
		if (cipher->ax[i] > qh) {
			cipher->ax[i] -= q;
		} else if (cipher->ax[i] < -qh) {
			cipher->ax[i] += q;
		}
		if (cipher->bx[i] > qh) {
			cipher->bx[i] -= q;
		} else if (cipher->bx[i] < -qh) {
			cipher->bx[i] += q;
		}
	}
}

void Scheme::coeffToSlotX0AndEqual(Ciphertext*& cipher) {
	long n0 = cipher->n0;
	long n1 = cipher->n1;

	long logn0 = log2(n0);
	long logn1 = log2(n1);

	long logk0 = logn0 / 2;
	long k0 = 1 << logk0;

	Ciphertext** rotvec = new Ciphertext*[k0];

	NTL_EXEC_RANGE(k0, first, last);
	for (long j = first; j < last; ++j) {
		rotvec[j] = (j == 0) ? new Ciphertext(cipher) : leftRotateFast(cipher, j, 0);
	}
	NTL_EXEC_RANGE_END;

	BootContext* bootContext = ring->bootContextMap.at({logn0, logn1});

	Ciphertext** tmpvec = new Ciphertext*[k0];

	NTL_EXEC_RANGE(k0, first, last);
	for (long j = first; j < last; ++j) {
		tmpvec[j] = multPolyNTTX0(rotvec[j], bootContext->rpxVec[j], bootContext->bndVec[j], bootContext->logp);
	}
	NTL_EXEC_RANGE_END;

	for (long j = 1; j < k0; ++j) {
		addAndEqual(tmpvec[0], tmpvec[j]);
	}

	delete cipher;
	cipher = new Ciphertext(tmpvec[0]);

	for (long ki = k0; ki < n0; ki += k0) {
		NTL_EXEC_RANGE(k0, first, last);
		for (long j = first; j < last; ++j) {
			delete tmpvec[j];
			tmpvec[j] = multPolyNTTX0(rotvec[j], bootContext->rpxVec[j + ki], bootContext->bndVec[j + ki], bootContext->logp);
		}
		NTL_EXEC_RANGE_END;
		for (long j = 1; j < k0; ++j) {
			addAndEqual(tmpvec[0], tmpvec[j]);
		}
		leftRotateFastAndEqual(tmpvec[0], ki, 0);
		addAndEqual(cipher, tmpvec[0]);
	}
	reScaleByAndEqual(cipher, bootContext->logp);
	for (long i = 0; i < k0; ++i) {
		delete rotvec[i];
	}
	delete[] rotvec;
	for (long i = 0; i < k0; ++i) {
		delete tmpvec[i];
	}
	delete[] tmpvec;
}

void Scheme::coeffToSlotX1AndEqual(Ciphertext*& cipher) {
	long n0 = cipher->n0;
	long n1 = cipher->n1;
	long logn0 = log2(n0);
	long logn1 = log2(n1);

	long logk1 = logn1 / 2;
	long k1 = 1 << logk1;

	BootContext* bootContext = ring->bootContextMap.at({logn0, logn1});

	Ciphertext* rot = new Ciphertext(cipher);
	for (long i = 1; i < ring->N1; i <<= 1) {
		Ciphertext* tmp = leftRotateFast(rot, 0, i);
		addAndEqual(rot, tmp);
		delete tmp;
	}

	Ciphertext** rotvec = new Ciphertext*[k1];

	NTL_EXEC_RANGE(k1, first, last);
	for (long j = first; j < last; ++j) {
		rotvec[j] = (j == 0) ? new Ciphertext(cipher) : leftRotateFast(cipher, 0, j);
	}
	NTL_EXEC_RANGE_END;

	Ciphertext** tmpvec = new Ciphertext*[k1];

	NTL_EXEC_RANGE(k1, first, last);
	for (long j = first; j < last; ++j) {
		complex<double> cnst = conj(ring->dftM1Pows[logn1][j]) * (double)n1/(double)ring->M1;
		tmpvec[j] = multConst(rotvec[j], cnst, bootContext->logp);
	}
	NTL_EXEC_RANGE_END;

	for (long j = 1; j < k1; ++j) {
		addAndEqual(tmpvec[0], tmpvec[j]);
	}
	delete cipher;
	cipher = new Ciphertext(tmpvec[0]);

	for (long ki = k1; ki < n1; ki += k1) {
		NTL_EXEC_RANGE(k1, first, last);
		for (long j = first; j < last; ++j) {
			complex<double> cnst = conj(ring->dftM1Pows[logn1][j + ki]) * (double)n1/(double)ring->M1;
			delete tmpvec[j];
			tmpvec[j] = multConst(rotvec[j], cnst, bootContext->logp);
		}
		NTL_EXEC_RANGE_END;
		for (long j = 1; j < k1; ++j) {
			addAndEqual(tmpvec[0], tmpvec[j]);
		}
		leftRotateFastAndEqual(tmpvec[0], 0, ki);
		addAndEqual(cipher, tmpvec[0]);
	}

	multConstAndEqual(rot, (double)n1/(double)ring->M1, bootContext->logp);
	subAndEqual(cipher, rot);
	reScaleByAndEqual(cipher, bootContext->logp);
	for (long i = 0; i < k1; ++i) {
		delete rotvec[i];
	}
	delete[] rotvec;
	for (long i = 0; i < k1; ++i) {
		delete tmpvec[i];
	}
	delete[] tmpvec;
	delete rot;
}

void Scheme::coeffToSlotAndEqual(Ciphertext*& cipher) {
	coeffToSlotX1AndEqual(cipher);
	coeffToSlotX0AndEqual(cipher);
}

void Scheme::slotToCoeffX0AndEqual(Ciphertext*& cipher) {
	long n0 = cipher->n0;
	long n1 = cipher->n1;
	long logn0 = log2(n0);
	long logn1 = log2(n1);

	long logk0 = logn0 / 2;
	long k0 = 1 << logk0;

	Ciphertext** rotvec = new Ciphertext*[k0];

	NTL_EXEC_RANGE(k0, first, last);
	for (long j = first; j < last; ++j) {
		rotvec[j] = (j == 0) ? new Ciphertext(cipher) : leftRotateFast(cipher, j, 0);
	}
	NTL_EXEC_RANGE_END;

	BootContext* bootContext = ring->bootContextMap.at({logn0, logn1});

	Ciphertext** tmpvec = new Ciphertext*[k0];

	NTL_EXEC_RANGE(k0, first, last);
	for (long j = first; j < last; ++j) {
		tmpvec[j] = multPolyNTTX0(rotvec[j], bootContext->rpxInvVec[j], bootContext->bndInvVec[j], bootContext->logp);
	}
	NTL_EXEC_RANGE_END;

	for (long j = 1; j < k0; ++j) {
		addAndEqual(tmpvec[0], tmpvec[j]);
	}
	delete cipher;
	cipher = new Ciphertext(tmpvec[0]);
	for (long ki = k0; ki < n0; ki+=k0) {
		NTL_EXEC_RANGE(k0, first, last);
		for (long j = first; j < last; ++j) {
			delete tmpvec[j];
			tmpvec[j] = multPolyNTTX0(rotvec[j], bootContext->rpxInvVec[j + ki], bootContext->bndInvVec[j + ki], bootContext->logp);
		}
		NTL_EXEC_RANGE_END;

		for (long j = 1; j < k0; ++j) {
			addAndEqual(tmpvec[0], tmpvec[j]);
		}

		leftRotateFastAndEqual(tmpvec[0], ki, 0);
		addAndEqual(cipher, tmpvec[0]);
	}
	reScaleByAndEqual(cipher, bootContext->logp);
	for (long i = 0; i < k0; ++i) {
		delete rotvec[i];
	}
	delete[] rotvec;
	for (long i = 0; i < k0; ++i) {
		delete tmpvec[i];
	}
	delete[] tmpvec;
}

void Scheme::slotToCoeffX1AndEqual(Ciphertext*& cipher) {
	long n0 = cipher->n0;
	long n1 = cipher->n1;
	long logn0 = log2(n0);
	long logn1 = log2(n1);

	long logk1 = logn1 / 2;
	long k1 = 1 << logk1;

	Ciphertext** rotvec = new Ciphertext*[k1];

	NTL_EXEC_RANGE(k1, first, last);
	for (long j = first; j < last; ++j) {
		rotvec[j] = (j == 0) ? new Ciphertext(cipher) : leftRotateFast(cipher, 0, j);
	}
	NTL_EXEC_RANGE_END;

	BootContext* bootContext = ring->bootContextMap.at({logn0, logn1});

	Ciphertext** tmpvec = new Ciphertext*[k1];

	NTL_EXEC_RANGE(k1, first, last);
	for (long j = first; j < last; ++j) {
		complex<double> cnst = ring->dftM1Pows[logn1][n1-j];
		tmpvec[j] = multConst(rotvec[j], cnst, bootContext->logp);
	}
	NTL_EXEC_RANGE_END;

	for (long j = 1; j < k1; ++j) {
		addAndEqual(tmpvec[0], tmpvec[j]);
	}

	delete cipher;
	cipher = new Ciphertext(tmpvec[0]);

	for (long ki = k1; ki < n1; ki+=k1) {
		NTL_EXEC_RANGE(k1, first, last);
		for (long j = first; j < last; ++j) {
			complex<double> cnst = ring->dftM1Pows[logn1][n1-j-ki];
			delete tmpvec[j];
			tmpvec[j] = multConst(rotvec[j], cnst, bootContext->logp);
		}
		NTL_EXEC_RANGE_END;

		for (long j = 1; j < k1; ++j) {
			addAndEqual(tmpvec[0], tmpvec[j]);
		}

		leftRotateFastAndEqual(tmpvec[0], 0, ki);
		addAndEqual(cipher, tmpvec[0]);
	}

	reScaleByAndEqual(cipher, bootContext->logp);
	for (long i = 0; i < k1; ++i) {
		delete rotvec[i];
	}
	delete[] rotvec;
	for (long i = 0; i < k1; ++i) {
		delete tmpvec[i];
	}
	delete[] tmpvec;
}

void Scheme::slotToCoeffAndEqual(Ciphertext*& cipher) {
	slotToCoeffX0AndEqual(cipher);
	slotToCoeffX1AndEqual(cipher);
}

void Scheme::exp2piAndEqual(Ciphertext* cipher, long logp) {
	Ciphertext* cipher2 = square(cipher);
	reScaleByAndEqual(cipher2, logp);

	Ciphertext* cipher4 = square(cipher2);
	reScaleByAndEqual(cipher4, logp);

	RR c = 1/(2*Pi);
	Ciphertext* cipher01 = addConst(cipher, c, logp);

	c = 2*Pi;
	multConstAndEqual(cipher01, c, logp);
	reScaleByAndEqual(cipher01, logp);

	c = 3/(2*Pi);
	Ciphertext* cipher23 = addConst(cipher, c, logp);

	c = 4*Pi*Pi*Pi/3;
	multConstAndEqual(cipher23, c, logp);
	reScaleByAndEqual(cipher23, logp);

	multAndEqual(cipher23, cipher2);
	reScaleByAndEqual(cipher23, logp);

	addAndEqual(cipher23, cipher01);

	c = 5/(2*Pi);
	Ciphertext* cipher45 = addConst(cipher, c, logp);

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
	delete cipher2; delete cipher4; delete cipher01; delete cipher23; delete cipher45;
}

void Scheme::removeIPartAndEqual(Ciphertext* cipher, long logT, long logI) {
	long logn0 = log2(cipher->n0);
	long logn1 = log2(cipher->n1);

	BootContext* bootContext = ring->bootContextMap.at({logn0, logn1});

	Ciphertext* tmp = conjugate(cipher);
	Ciphertext* cimag = sub(cipher, tmp);
	addAndEqual(cipher, tmp);
	imultAndEqual(cipher);

	divPo2AndEqual(cipher, logT + logn0 + logn1 + 1);
	divPo2AndEqual(cimag, logT +  logn0 + logn1 + 1);

	exp2piAndEqual(cipher, bootContext->logp);
	exp2piAndEqual(cimag, bootContext->logp);

	for (long i = 0; i < logI + logT; ++i) {
		squareAndEqual(cipher);
		squareAndEqual(cimag);
		reScaleByAndEqual(cipher, bootContext->logp);
		reScaleByAndEqual(cimag, bootContext->logp);
	}

	delete tmp;
	tmp = conjugate(cimag);
	subAndEqual(cimag, tmp);
	delete tmp;
	tmp = conjugate(cipher);
	subAndEqual(cipher, tmp);
	imultAndEqual(cipher);
	subAndEqual2(cimag, cipher);

	RR c = 0.25/Pi;
	multConstAndEqual(cipher, c, bootContext->logp);
	reScaleByAndEqual(cipher, bootContext->logp + logI + ring->logN1-logn1);
	delete tmp; delete cimag;
}

void Scheme::bootstrapX0AndEqual(Ciphertext*& cipher, long logq, long logQ, long logT, long logI) {
	long logn0 = log2(cipher->n0);
	long logp = cipher->logp;

	modDownToAndEqual(cipher, logq);
	normalizeAndEqual(cipher);

	cipher->logq = logQ;
	cipher->logp = logq + logI;

	for (long i = logn0; i < ring->logN0h; ++i) {
		Ciphertext* rot = leftRotateFast(cipher, (1 << i), 0);
		addAndEqual(cipher, rot);
		delete rot;
	}
	coeffToSlotX0AndEqual(cipher);
	divPo2AndEqual(cipher, ring->logN0h);
	removeIPartAndEqual(cipher, logT, logI);
	slotToCoeffX0AndEqual(cipher);
	cipher->logp = logp;
}

void Scheme::bootstrapX1AndEqual(Ciphertext*& cipher, long logq, long logQ, long logT, long logI) {
	long logn1 = log2(cipher->n1);
	long logp = cipher->logp;

	modDownToAndEqual(cipher, logq);
	normalizeAndEqual(cipher);

	cipher->logq = logQ;
	cipher->logp = logq + logI;

	for (long i = logn1; i < ring->logN1; ++i) {
		Ciphertext* rot = leftRotateFast(cipher, 0, (1 << i));
		addAndEqual(cipher, rot);
		delete rot;
	}
	coeffToSlotX1AndEqual(cipher);
	divPo2AndEqual(cipher, ring->logN1);
	removeIPartAndEqual(cipher, logT, logI);
	slotToCoeffX1AndEqual(cipher);
	cipher->logp = logp;
}

void Scheme::bootstrapAndEqual(Ciphertext*& cipher, long logq, long logQ, long logT, long logI) {
	long logn0 = log2(cipher->n0);
	long logn1 = log2(cipher->n1);
	long logp = cipher->logp;

	modDownToAndEqual(cipher, logq);
	normalizeAndEqual(cipher);

	cipher->logq = logQ;
	cipher->logp = logq + logI;

	for (long i = logn0; i < ring->logN0h; ++i) {
		Ciphertext* rot = leftRotateFast(cipher, (1 << i), 0);
		addAndEqual(cipher, rot);
		delete rot;
	}
	for (long i = logn1; i < ring->logN1; ++i) {
		Ciphertext* rot = leftRotateFast(cipher, 0, (1 << i));
		addAndEqual(cipher, rot);
		delete rot;
	}
	coeffToSlotAndEqual(cipher);
	divPo2AndEqual(cipher, ring->logN0h + ring->logN1);
	removeIPartAndEqual(cipher, logT, logI);
	slotToCoeffAndEqual(cipher);
	cipher->logp = logp;
}
