#pragma once

#define MRK(x) (0xFF00 | (x))

#define BDCT_SOF 0xC0
#define ESDCT_SOF 0xC1
#define EPDCT_SOF 0xC2
#define DHT 0xC4
#define RST(m) (0xD0 | (m))
#define SOI 0xD8
#define EOI 0xD9
#define SOS 0xDA
#define DQT 0xDB
#define DNL 0xDC
#define DRI 0xDD
#define APP(n) (0xE0 | (n))
