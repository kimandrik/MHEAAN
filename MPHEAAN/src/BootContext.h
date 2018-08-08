/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef MPHEAAN_BOOTCONTEXT_H_
#define MPHEAAN_BOOTCONTEXT_H_

#include <NTL/ZZ.h>

using namespace NTL;

class BootContext {
public:

	ZZ** pxVec;
	ZZ** pyrVec;
	ZZ** pyiVec;

	ZZ** pxInvVec;
	ZZ** pyrInvVec;
	ZZ** pyiInvVec;

	ZZ* p1;
	ZZ* p2;

	long logp;

	BootContext(ZZ** pxVec = NULL, ZZ** pyrVec = NULL, ZZ** pyiVec = NULL, ZZ** pxInvVec = NULL, ZZ** pyrInvVec = NULL, ZZ** pyiInvVec = NULL, ZZ* p1 = NULL, ZZ* p2 = NULL, long logp = 0);

};

#endif /* BOOTCONTEXT_H_ */
