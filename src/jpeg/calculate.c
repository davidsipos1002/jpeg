#include <jpeg/calculate.h>

#include <math.h>
#include <stdlib.h>
#include <mem.h>

static uint8_t zigzag[8][8] = 
{
    0, 1, 5, 6, 14, 15, 27, 28,
    2, 4, 7, 13, 16, 26, 29, 42,
    3, 8, 12, 17, 25, 30, 41, 43,
    9, 11, 18, 24, 31, 40, 44, 53,
    10, 19, 23, 32, 39, 45, 52, 54,
    20, 22, 33, 38, 46, 51, 55, 60,
    21, 34, 37, 47, 50, 56, 59, 61,
    35, 36, 48, 49, 57, 58, 62, 63
};

static uint8_t izigzag[] =
{
	// k = 0 i = 0 j = 0
	0x0,
	// k = 1 i = 0 j = 1
	0x1,
	// k = 2 i = 1 j = 0
	0x10,
	// k = 3 i = 2 j = 0
	0x20,
	// k = 4 i = 1 j = 1
	0x11,
	// k = 5 i = 0 j = 2
	0x2,
	// k = 6 i = 0 j = 3
	0x3,
	// k = 7 i = 1 j = 2
	0x12,
	// k = 8 i = 2 j = 1
	0x21,
	// k = 9 i = 3 j = 0
	0x30,
	// k = 10 i = 4 j = 0
	0x40,
	// k = 11 i = 3 j = 1
	0x31,
	// k = 12 i = 2 j = 2
	0x22,
	// k = 13 i = 1 j = 3
	0x13,
	// k = 14 i = 0 j = 4
	0x4,
	// k = 15 i = 0 j = 5
	0x5,
	// k = 16 i = 1 j = 4
	0x14,
	// k = 17 i = 2 j = 3
	0x23,
	// k = 18 i = 3 j = 2
	0x32,
	// k = 19 i = 4 j = 1
	0x41,
	// k = 20 i = 5 j = 0
	0x50,
	// k = 21 i = 6 j = 0
	0x60,
	// k = 22 i = 5 j = 1
	0x51,
	// k = 23 i = 4 j = 2
	0x42,
	// k = 24 i = 3 j = 3
	0x33,
	// k = 25 i = 2 j = 4
	0x24,
	// k = 26 i = 1 j = 5
	0x15,
	// k = 27 i = 0 j = 6
	0x6,
	// k = 28 i = 0 j = 7
	0x7,
	// k = 29 i = 1 j = 6
	0x16,
	// k = 30 i = 2 j = 5
	0x25,
	// k = 31 i = 3 j = 4
	0x34,
	// k = 32 i = 4 j = 3
	0x43,
	// k = 33 i = 5 j = 2
	0x52,
	// k = 34 i = 6 j = 1
	0x61,
	// k = 35 i = 7 j = 0
	0x70,
	// k = 36 i = 7 j = 1
	0x71,
	// k = 37 i = 6 j = 2
	0x62,
	// k = 38 i = 5 j = 3
	0x53,
	// k = 39 i = 4 j = 4
	0x44,
	// k = 40 i = 3 j = 5
	0x35,
	// k = 41 i = 2 j = 6
	0x26,
	// k = 42 i = 1 j = 7
	0x17,
	// k = 43 i = 2 j = 7
	0x27,
	// k = 44 i = 3 j = 6
	0x36,
	// k = 45 i = 4 j = 5
	0x45,
	// k = 46 i = 5 j = 4
	0x54,
	// k = 47 i = 6 j = 3
	0x63,
	// k = 48 i = 7 j = 2
	0x72,
	// k = 49 i = 7 j = 3
	0x73,
	// k = 50 i = 6 j = 4
	0x64,
	// k = 51 i = 5 j = 5
	0x55,
	// k = 52 i = 4 j = 6
	0x46,
	// k = 53 i = 3 j = 7
	0x37,
	// k = 54 i = 4 j = 7
	0x47,
	// k = 55 i = 5 j = 6
	0x56,
	// k = 56 i = 6 j = 5
	0x65,
	// k = 57 i = 7 j = 4
	0x74,
	// k = 58 i = 7 j = 5
	0x75,
	// k = 59 i = 6 j = 6
	0x66,
	// k = 60 i = 5 j = 7
	0x57,
	// k = 61 i = 6 j = 7
	0x67,
	// k = 62 i = 7 j = 6
	0x76,
	// k = 63 i = 7 j = 7
	0x77,
};

// These functions are just for initial functionality testing
// They will be rewritten to use ARM Neon and Intel AVX

// this function generates izigzag from the zigzag matrix
void gen_inverse_zigzag()
{
    printf("static uint8_t izigzag[] =\n");
    printf("{\n");
    for (int k = 0; k < 64; k++)
    {
        int pi, pj;
        for (int i = 0; i < 8; i++)
        {
            for (int j = 0; j < 8; j++)
            {
                if (zigzag[i][j] == k)
                {
                    pi = i;
                    pj = j;
                }
            }
        }
        printf("\t// k = %u i = %u j = %u\n", k, pi, pj);
        printf("\t0x%X,\n", (pi << 4) + pj);
    }
    printf("};\n");
}

// this function generates the above cos table
void gen_cos_table()
{
	printf("#ifdef __aarch64__\n");
	printf("static float16_t costable[8][8] = \n");
	printf("{\n");
	for (uint8_t i = 0; i < 8; i++)
	{
		for (uint8_t j = 0; j < 8; j++)
		{	
		}
	}
	printf("};\n");
	printf("#endif");
}

#ifdef __aarch64__
void dequantize(int16_t *b, int16_t *q)
{
	// we will do the multiplication using ARM Neon
	// we split the input in 8 integer blocks
	for (uint8_t i = 0; i < 8; i++)
	{
		int16_t *bp = b + i * 8;
		int16_t *qp = q + i * 8;
		// load the elements in of the 128-bit Q registers
		int16x8_t vb = vld1q_s16(bp);
		int16x8_t vq = vld1q_s16(qp);
		// do the multiplication
		int16x8_t vdq = vmulq_s16(vb, vq);
		// store the result
		vst1q_s16(bp, vdq);
	}
}
#else
void dequantize(int16_t *b, int16_t *q)
{
	for (uint8_t i = 0; i < 64; i++)
		b[i] *= q[i];
}
#endif

void unzigzag(int16_t *b, jpegf **mat)
{
	for (uint8_t i = 0; i < 64; i++)
	{
		uint8_t pos = izigzag[i];
		mat[pos >> 4][pos & 0xF] = b[i];
	}
}

static jpegf idct_helper(jpegf **mat, uint8_t ii, uint8_t jj)
{
	float ret = 0;
	for (uint8_t i = 0; i < 8; i++)
	{
		for (uint8_t j = 0; j < 8; j++)
		{
			jpegf ci = i == 0 ? (1 / sqrt(2)) : 1; 
			jpegf cj = j == 0 ? (1 / sqrt(2)) : 1;	
			jpegf cosi = cos((double) (2 * ii + 1) * i * M_PI / 16);
			jpegf cosj = cos((double) (2 * jj + 1) * j * M_PI / 16);
			ret += ci * cj * mat[i][j] * cosi * cosj;
		}
	}
	ret /= 4;
	return ret;
}

// defined in Annex A 3.3 of the standard
void idct(jpegf **mat)
{
	float res[8][8];
	for (uint8_t i = 0; i < 8; i++)
		for (uint8_t j = 0; j < 8; j++)
			res[i][j] = idct_helper(mat, i, j);
	for (uint8_t i = 0; i < 8; i++)
		for (uint8_t j = 0; j < 8; j++)
			mat[i][j] = res[i][j];
}

void undo_level_shift(jpegf **mat)
{
	for (uint8_t i = 0; i < 8; i++)
		for (uint8_t j = 0; j < 8; j++)
			mat[i][j] += 128;
}

static jpegf get_min(jpegf x, jpegf y)
{
	if (x < y)
		return x;
	return y;
}

static jpegf get_max(jpegf x, jpegf y)
{
	if (x < y)
		return y;
	return x;
}

// defined in JFIF standard
static jpegf get_round(jpegf x)
{
	return floor(x + 0.5);
}

// formulas taken from the JFIF standard
// R = Min(Max(0,Round(Y +1.402*(CR −128) )),255)
// G = Min(Max(0,Round(Y −(0.114*1.772*(CB −128)+0.299*1.402*(CR −128))/0.587)),255)
// B = Min(Max(0,Round(Y +1.772*(CB −128) )),255)
void convert_to_rgb(jpegf y, jpegf cb, jpegf cr, uint8_t *rp, uint8_t *gp, uint8_t *bp)
{
	jpegf r = get_min(get_max(0, get_round(y + 1.402 * (cr - 128))), 255); 	
	jpegf g = get_min(get_max(0, get_round(y - (0.114 * 1.772 * (cb - 128) + 0.299 * 1.402 * (cr - 128))/0.587)), 255);
	jpegf b = get_min(get_max(0, get_round(y + 1.772 * (cb - 128))), 255);
	*rp = (uint8_t) r;
	*gp = (uint8_t) g;
	*bp = (uint8_t) b;
}
