/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#ifndef MHEAAN_SECRETKEY_H_
#define MHEAAN_SECRETKEY_H_

#include <NTL/ZZ.h>
#include "Ring2XY.h"

using namespace std;
using namespace NTL;

class SecretKey {
public:

	long Nx;
	long Ny;

	ZZ* sxy; ///< secret key

	SecretKey(Ring2XY& ring);

};

#endif
