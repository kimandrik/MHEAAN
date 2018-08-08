/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "BootContext.h"

BootContext::BootContext(ZZ** pxVec, ZZ** pyrVec, ZZ** pyiVec, ZZ** pxInvVec, ZZ** pyrInvVec, ZZ** pyiInvVec, ZZ* p1, ZZ* p2, long logp)
			: pxVec(pxVec), pyrVec(pyrVec), pyiVec(pyiVec), pxInvVec(pxInvVec), pyrInvVec(pyrInvVec), pyiInvVec(pyiInvVec), p1(p1), p2(p2), logp(logp) {
}
