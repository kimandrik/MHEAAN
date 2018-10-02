/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef MPHEAAN_MATRIXCONTEXT_H_
#define MPHEAAN_MATRIXCONTEXT_H_

#include <NTL/ZZ.h>

using namespace NTL;

class SqrMatContext {
public:

	ZZ** mvec;
	long* bndVec;
	long logp;

	SqrMatContext(ZZ** mvec, long logp);
};

#endif /* MATRIXCONTEXT_H_ */
