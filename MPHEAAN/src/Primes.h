/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#ifndef MPHEAAN_PRIMES_H_
#define MPHEAAN_PRIMES_H_

#include <cstdint>

using namespace std;

/**
 * Primes are of form ...
 */
const uint64_t pPrimesVec[270] = {
		0x808000000a8a801,
		0x808000001c1c001,
		0x8080000020a0801,
		0x808000002424001,
		0x8080000026a6801,
		0x808000004848001,
		0x808000004c4c001,
		0x808000004e4e001,
		0x808000005b5b001,
		0x808000005d5d001,
		0x808000008302801,
		0x808000008a8a001,
		0x808000009393001,
		0x808000009797001,
		0x80800000b1b1001,
		0x80800000b4b4001,
		0x80800000b5b5001,
		0x80800000d3d3001,
		0x80800000f372801,
		0x808000011f1e001,
		0x808000012120001,
		0x80800001302f001,
		0x808000013433001,
		0x808000013d3c001,
		0x808000014948001,
		0x8080000153d2801,
		0x808000015a59001,
		0x808000016463001,
		0x8080000168e7801,
		0x808000016bea801,
		0x808000017574001,
		0x8080000177f6801,
		0x808000017cfb801,
		0x808000017f7e001,
		0x808000018786001,
		0x808000019c9b001,
		0x80800001a6a5001,
		0x80800001b432801,
		0x80800001b7b6001,
		0x80800001c1c0001,
		0x80800001cc4a801,
		0x80800001e15f801,
		0x80800001e2e1001,
		0x80800001e4e3001,
		0x80800001fe7c801,
		0x808000020280801,
		0x808000020583801,
		0x808000020a08001,
		0x808000020c0a001,
		0x808000020d8b801,
		0x8080000222a0801,
		0x8080000226a4801,
		0x808000024644001,
		0x8080000249c7801,
		0x808000024fcd801,
		0x808000025250001,
		0x8080000253d1801,
		0x8080000258d6801,
		0x808000025d5b001,
		0x80800002605e001,
		0x80800002615f001,
		0x808000026765001,
		0x808000026967001,
		0x808000027dfb801,
		0x808000027e7c001,
		0x808000028482001,
		0x808000029f1c801,
		0x80800002b1af001,
		0x80800002bd3a801,
		0x80800002d24f801,
		0x80800002d2d0001,
		0x80800002dd5a801,
		0x80800002e2e0001,
		0x80800002e764801,
		0x80800002ea67801,
		0x80800002f875801,
		0x80800002fe7b801,
		0x808000030603001,
		0x808000032ba8801,
		0x8080000337b4801,
		0x808000033a37001,
		0x80800003423f001,
		0x8080000349c6801,
		0x80800003524f001,
		0x8080000359d6801,
		0x80800003605d001,
		0x8080000361de801,
		0x8080000367e4801,
		0x808000036be8801,
		0x8080000376f3801,
		0x808000037e7b001,
		0x808000038f0b801,
		0x808000039c18801,
		0x80800003a8a5001,
		0x80800003ba36801,
		0x80800003bd39801,
		0x80800003c4c1001,
		0x80800003c844801,
		0x80800003dd59801,
		0x80800003e15d801,
		0x80800003ef6b801,
		0x80800003fb77801,
		0x80800003fcf9001,
		0x808000040703001,
		0x80800004100c001,
		0x80800004118d801,
		0x808000041915001,
		0x808000041b17001,
		0x808000042a26001,
		0x808000042ca8801,
		0x8080000431ad801,
		0x808000043935001,
		0x80800004423e001,
		0x80800004433f001,
		0x808000044642001,
		0x808000045fdb801,
		0x808000046f6b001,
		0x808000047e7a001,
		0x8080000480fc801,
		0x808000048601801,
		0x808000049b16801,
		0x80800004a39f001,
		0x80800004a823801,
		0x80800004a8a4001,
		0x80800004b7b3001,
		0x80800004bd38801,
		0x80800004da55801,
		0x80800004dfdb001,
		0x80800004e1dd001,
		0x80800004e762801,
		0x80800004e8e4001,
		0x80800004ef6a801,
		0x8080000502fe001,
		0x808000051c97801,
		0x80800005241f001,
		0x808000053eb9801,
		0x8080000541bc801,
		0x80800005615c001,
		0x808000056762001,
		0x808000057af5801,
		0x808000057e79001,
		0x8080000582fd801,
		0x808000058600801,
		0x808000058a85001,
		0x808000058d07801,
		0x80800005918c001,
		0x80800005948f001,
		0x808000059e18801,
		0x80800005a09b001,
		0x80800005a41e801,
		0x80800005aa24801,
		0x80800005b42e801,
		0x80800005da54801,
		0x80800005e1dc001,
		0x80800005f670801,
		0x80800005fcf7001,
		0x80800005fdf8001,
		0x80800006148e801,
		0x808000061610001,
		0x8080000628a2801,
		0x8080000631ab801,
		0x8080000634ae801,
		0x808000064e48001,
		0x8080000653cd801,
		0x808000065b55001,
		0x808000066de7801,
		0x808000067e78001,
		0x80800006827c001,
		0x808000069e17801,
		0x80800006a8a2001,
		0x80800006aa23801,
		0x80800006cfc9001,
		0x80800006d9d3001,
		0x80800006da53801,
		0x80800006db54801,
		0x80800006e760801,
		0x80800006faf4001,
		0x80800007047d801,
		0x80800007057e801,
		0x8080000705ff001,
		0x80800007150e001,
		0x808000071e17001,
		0x8080000728a1801,
		0x808000072b24001,
		0x808000072e27001,
		0x8080000732ab801,
		0x808000073ab3801,
		0x808000073c35001,
		0x8080000741ba801,
		0x80800007433c001,
		0x808000074fc8801,
		0x808000075a53001,
		0x8080000764dd801,
		0x808000077972001,
		0x808000077af3801,
		0x8080000780f9801,
		0x80800007857e001,
		0x80800007948d001,
		0x80800007960e801,
		0x80800007a099001,
		0x80800007a51d801,
		0x80800007b7b0001,
		0x80800007beb7001,
		0x80800007bf37801,
		0x80800007ce46801,
		0x80800007e45c801,
		0x80800007ea62801,
		0x80800007f36b801,
		0x80800007f66e801,
		0x80800007fdf6001,
		0x808000080c04001,
		0x808000081911001,
		0x808000082b23001,
		0x80800008524a001,
		0x8080000858d0801,
		0x8080000861d9801,
		0x808000086c64001,
		0x808000086de5801,
		0x80800008756d001,
		0x8080000877ef801,
		0x8080000885fd801,
		0x80800008877f001,
		0x808000088d04801,
		0x808000088f06801,
		0x80800008960d801,
		0x80800008a39b001,
		0x80800008a41b801,
		0x80800008beb6001,
		0x80800008cac2001,
		0x80800008dcd4001,
		0x80800008fe75801,
		0x80800009138a801,
		0x808000091c13001,
		0x808000092299801,
		0x80800009241b001,
		0x80800009281f001,
		0x808000092a21001,
		0x8080000935ac801,
		0x808000093930001,
		0x808000093d34001,
		0x808000093f36001,
		0x8080000950c7801,
		0x8080000955cc801,
		0x8080000967de801,
		0x8080000971e8801,
		0x80800009766d001,
		0x8080000986fd801,
		0x808000098c02801,
		0x808000098d84001,
		0x808000098e85001,
		0x808000099208801,
		0x80800009938a001,
		0x80800009a69d001,
		0x80800009aba2001,
		0x80800009b2a9001,
		0x80800009c238801,
		0x80800009d3ca001,
		0x80800009d5cc001,
		0x80800009e8df001,
		0x80800009f369801,
		0x80800009f56b801,
		0x80800009f7ee001,
		0x80800009fe74801,
		0x8080000a0d83801,
		0x8080000a148a801,
		0x8080000a178d801,
		0x8080000a1b11001,
		0x8080000a1d93801,
		0x8080000a259b801,
		0x8080000a281e001};

#endif /* PRIMES_H_ */
