#include <jpeg/calculate.h>

#include <math.h>
#include <stdlib.h>
#include <mem.h>
#include <misc.h>

// lookup tables

static const uint8_t fzigzag[8][8] = 
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

static const uint8_t izigzag[] =
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

static const jpegf icostable[8][8] = 
{
	1.000000, 0.980957, 0.923828, 0.831543, 0.707031, 0.555664, 0.382568, 0.195068, 
	1.000000, 0.831543, 0.382568, -0.195068, -0.707031, -0.980957, -0.923828, -0.555664, 
	1.000000, 0.555664, -0.382568, -0.980957, -0.707031, 0.195068, 0.923828, 0.831543, 
	1.000000, 0.195068, -0.923828, -0.555664, 0.707031, 0.831543, -0.382568, -0.980957, 
	1.000000, -0.195068, -0.923828, 0.555664, 0.707031, -0.831543, -0.382568, 0.980957, 
	1.000000, -0.555664, -0.382568, 0.980957, -0.707031, -0.195068, 0.923828, -0.831543, 
	1.000000, -0.831543, 0.382568, 0.195068, -0.707031, 0.980957, -0.923828, 0.555664, 
	1.000000, -0.980957, 0.923828, -0.831543, 0.707031, -0.555664, 0.382568, -0.195068, 
};

static const jpegf fcostable[8][8] = 
{
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    0.980957, 0.831543, 0.555664, 0.195068, -0.195068, -0.555664, -0.831543, -0.980957, 
    0.923828, 0.382568, -0.382568, -0.923828, -0.923828, -0.382568, 0.382568, 0.923828, 
    0.831543, -0.195068, -0.980957, -0.555664, 0.555664, 0.980957, 0.195068, -0.831543, 
    0.707031, -0.707031, -0.707031, 0.707031, 0.707031, -0.707031, -0.707031, 0.707031, 
    0.555664, -0.980957, 0.195068, 0.831543, -0.831543, -0.195068, 0.980957, -0.555664, 
    0.382568, -0.923828, 0.923828, -0.382568, -0.382568, 0.923828, -0.923828, 0.382568, 
    0.195068, -0.555664, 0.831543, -0.980957, 0.980957, -0.831543, 0.555664, -0.195068, 
};

static jpegf c[8] = {0.7071067812, 1, 1, 1, 1, 1, 1, 1};

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
                if (fzigzag[i][j] == k)
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
	printf("static jpegf icostable[8][8] = \n");
	printf("{\n");
	for (uint8_t i = 0; i < 8; i++)
	{
		printf("\t");
		for (uint8_t j = 0; j < 8; j++)
			printf("%f, ", (jpegf) cos((double) (2 * i + 1) * j * M_PI / 16));
		printf("\n");
	}
	printf("};\n");

	printf("static jpegf fcostable[8][8] = \n");
	printf("{\n");
	for (uint8_t i = 0; i < 8; i++)
	{
		printf("\t");
		for (uint8_t j = 0; j < 8; j++)
			printf("%f, ", (jpegf) cos((double) (2 * j + 1) * i * M_PI / 16));
		printf("\n");
	}
	printf("};\n");
}

#ifdef USE_NEON
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

#ifdef USE_NEON
static jpegf idct_helper(jpegf **mat, uint8_t ii, uint8_t jj)
{
	// row normalization constants
	float16x8_t rowc;
	// cosine value for a row
	float16x8_t rowcos;
	// accumulated columns sums
	float16x8_t sumcol = vmovq_n_f16(0);
	// normalization constants for columns of a row
	float16x8_t colc;
	// cosine values for columns of a row
	float16x8_t colcos;
	// matrix row
	float16x8_t matrow;
	// partial result accumulator
	float16x8_t part;
	// final 4 values to be added see below
	float16_t r[4];

	// compute innersum for each row and accumulate the sums per columns in sumcol
	for (uint8_t i = 0; i < 8; i++)
	{
		// load the values
		rowc = vmovq_n_f16(c[i]);
		rowcos = vmovq_n_f16(icostable[ii][i]);
		colc = vld1q_f16(c);	
		colcos = vld1q_f16(icostable[jj]);
		matrow = vld1q_f16(mat[i]);

		// compute innersum elements
		part = vmulq_f16(rowc, rowcos);
		part = vmulq_f16(part, colc);
		part = vmulq_f16(part, colcos);
		part = vmulq_f16(part, matrow);
		
		// acummulate the innersum elements per column
		sumcol = vaddq_f16(sumcol, part);
	}
	
	// we need to add the 8 floats in sumcol 
	// and divide by 4 to get the final result
	
	// split the 8 floats into two groups of 4 floats
	float16x4_t hsumcol = vget_high_f16(sumcol);
	float16x4_t lsumcol = vget_low_f16(sumcol);
	float16x4_t sum = vadd_f16(lsumcol, hsumcol);
	vst1_f16(r, sum);
	float16_t finalsum = r[0] + r[1] + r[2] + r[3];
	return finalsum / 4;
}
#else
static jpegf idct_helper(jpegf **mat, uint8_t ii, uint8_t jj)
{
	float ret = 0;
	for (uint8_t i = 0; i < 8; i++)
	{
		for (uint8_t j = 0; j < 8; j++)
		{
			jpegf cosi = costable[ii][i];
			jpegf cosj = costable[jj][j];
			ret += c[i] * c[j] * mat[i][j] * cosi * cosj;
		}
	}
	ret /= 4;
	return ret;
}
#endif

// defined in Annex A 3.3 of the standard
void idct(jpegf **mat)
{
	jpegf res[8][8];
	for (uint8_t i = 0; i < 8; i++)
		for (uint8_t j = 0; j < 8; j++)
			res[i][j] = idct_helper(mat, i, j);
	for (uint8_t i = 0; i < 8; i++)
		for (uint8_t j = 0; j < 8; j++)
			mat[i][j] = res[i][j];
}

#ifdef USE_NEON
void undo_level_shift(jpegf **mat)
{
	// load constant 128 in of the 128-bit Q registers
	float16x8_t shift = vmovq_n_f16(128);
	for (uint8_t i = 0; i < 8;i ++)
	{
		float16_t *rowptr = mat[i];
		// load row to Q register
		float16x8_t rw = vld1q_f16(rowptr);
		// add the two registers
		float16x8_t shifted = vaddq_f16(rw, shift);
		// store result
		vst1q_f16(rowptr, shifted);
	}
}
#else
void undo_level_shift(jpegf **mat)
{
	for (uint8_t i = 0; i < 8; i++)
		for (uint8_t j = 0; j < 8; j++)
			mat[i][j] += 128;
}
#endif

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
#ifdef USE_NEON
void convert_to_rgb(jpegf y, jpegf cb, jpegf cr, uint8_t *rp, uint8_t *gp, uint8_t *bp)
{
	static const float16_t cbc[] = {0, -0.3441, 1.772, 0};
	static const float16_t crc[] = {1.402, -0.7141, 0, 0};

	float16x4_t cst = vmov_n_f16(128);
	float16x4_t vy = vmov_n_f16(y);
	float16x4_t vcb = vmov_n_f16(cb);
	float16x4_t vcr = vmov_n_f16(cr);
	float16x4_t vcbc = vld1_f16(cbc);
	float16x4_t vcrc = vld1_f16(crc);
	float16x4_t res = vmov_n_f16(0);
	float16_t ret[4];
	
	vcb = vsub_f16(vcb, cst);
	vcb = vmul_f16(vcb, vcbc);
	vcr = vsub_f16(vcr, cst);
	vcr = vmul_f16(vcr, vcrc);
	
	res = vadd_f16(res, vy);
	res = vadd_f16(res, vcb);
	res = vadd_f16(res, vcr);
	res = vrnd_f16(res);
	cst = vmov_n_f16(0);
	res = vmax_f16(res, cst);
	cst = vmov_n_f16(255);
	res = vmin_f16(res, cst); 

	vst1_f16(ret, res);
	*rp = ret[0];
	*gp = ret[1];
	*bp = ret[2];
}
#else
void convert_to_rgb(jpegf y, jpegf cb, jpegf cr, uint8_t *rp, uint8_t *gp, uint8_t *bp)
{
	jpegf r = get_min(get_max(0, get_round(y + 1.402 * (cr - 128))), 255); 	
	jpegf g = get_min(get_max(0, get_round(y - (0.114 * 1.772 * (cb - 128) + 0.299 * 1.402 * (cr - 128))/0.587)), 255);
	jpegf b = get_min(get_max(0, get_round(y + 1.772 * (cb - 128))), 255);
	*rp = (uint8_t) r;
	*gp = (uint8_t) g;
	*bp = (uint8_t) b;
}
#endif

// formulas taken from the JFIF standard
// Y = Min(Max(0,Round( 0.299*R +0.587*G +0.114*B )),255)
// CB = Min(Max(0,Round(( −0.299*R −0.587*G +0.886*B )/1.772 +128 )),255)
// CR = Min(Max(0,Round(( 0.701*R −0.587*G −0.114*B )/1.402 +128 )),255)
#ifdef USE_NEON
void convert_to_ycbcr(jpegf r, jpegf g, jpegf b, jpegf *yp, jpegf *cbp, jpegf *crp)
{
	static const float16_t cr[] = {0.299, -0.1687, 0.5, 0};
	static const float16_t cg[]	= {0.587, -0.3313, -0.4187, 0};
	static const float16_t cb[] = {0.114, 0.5, -0.0813, 0};
	static const float16_t sh[] = {0, 128, 128, 0};
	
	// load values
	float16x4_t vcr = vld1_f16(cr);
	float16x4_t vcg = vld1_f16(cg);
	float16x4_t vcb = vld1_f16(cb);
	float16x4_t vsh = vld1_f16(sh);
	float16x4_t vr = vmov_n_f16(r);
	float16x4_t vg = vmov_n_f16(g);
	float16x4_t vb = vmov_n_f16(b);
	
	// do the computation
	vr = vmul_f16(vcr, vr);
	vg = vmul_f16(vcg, vg);
	vb = vmul_f16(vcb, vb);
	
	float16x4_t s = vadd_f16(vr, vg);
	s = vadd_f16(s, vb);
	s = vadd_f16(s, vsh);
	s = vrnd_f16(s);
	float16x4_t cst = vmov_n_f16(0);
	s = vmax_f16(cst, s);
	cst = vmov_n_f16(255);
	s = vmin_f16(s, cst);
	
	// get the result
	float16_t res[4];
	vst1_f16(res, s);
	*yp = res[0];
	*cbp = res[1];
	*crp = res[2];
}
#else
void convert_to_ycbcr(jpegf r, jpegf g, jpegf b, jpegf *yp, jpegf *cbp, jpegf *crp)
{
	jpegf y = get_min(get_max(0, get_round(0.299 * r + 0.587 * g + 0.114 * b)), 255);
	jpegf cb = get_min(get_max(0, get_round((-0.299 * r - 0.587 * g + 0.886 * b) / 1.772 + 128)), 255);
	jpegf cr = get_min(get_max(0, get_round((0.701 * r - 0.587 * g - 0.114 * b) / 1.402 + 128)), 255);
	*yp = y;
	*cbp = cb;
	*crp = cr;
}
#endif

// defined in Annex A 3.3 of the standard
#ifdef USE_NEON
static jpegf fdct_helper(jpegf **mat, uint8_t ii, uint8_t jj)
{
	float16x8_t crow = vmovq_n_f16(c[ii]);	
	float16x8_t ccol = vmovq_n_f16(c[jj]);	
	float16x8_t sumcol = vmovq_n_f16(0);
	float16x8_t part;
	float16_t r[4];

	// compute the innermost sum per column
	for (uint8_t i = 0; i < 8; i++)
	{
		float16x8_t rowcos = vmovq_n_f16(fcostable[ii][i]);
		float16x8_t colcos = vld1q_f16(fcostable[jj]);
		float16x8_t matrow = vld1q_f16(mat[i]);
		
		part = vmulq_f16(crow, ccol);
		part = vmulq_f16(part, matrow);
		part = vmulq_f16(part, rowcos);	
		part = vmulq_f16(part, colcos);
		 
		sumcol = vaddq_f16(sumcol, part);
	}
	
	float16x4_t hsumcol = vget_high_f16(sumcol);
	float16x4_t lsumcol = vget_low_f16(sumcol); 
	float16x4_t sum = vadd_f16(lsumcol, hsumcol);
	vst1_f16(r, sum);
	
	float16_t ret = r[0] + r[1] + r[2] + r[3];
	return ret / 4;
}
#else
static jpegf fdct_helper(jpegf **mat, uint8_t ii, uint8_t jj)
{
	jpegf ret = 0;
	for (uint8_t i = 0; i < 8; i++)
	{
		for (uint8_t j = 0; j < 8; j++)
		{
			jpegf cosi = icostable[i][ii];	
			jpegf cosj = icostable[j][jj];	
			ret += c[ii] * c[jj] * mat[i][j] * cosi * cosj;
		}
	}	
	ret /= 4;
	return ret;
}
#endif

void fdct(jpegf **mat)
{
	jpegf res[8][8];
	for (uint8_t i = 0; i < 8; i++)
		for (uint8_t j = 0; j < 8; j++)
			res[i][j] = fdct_helper(mat, i, j);
	for (uint8_t i = 0; i < 8; i++)
		for (uint8_t j = 0; j < 8; j++)
			mat[i][j] = res[i][j];
}


void zigzag(jpegf **mat, int16_t *b)
{
	for (uint8_t i = 0; i < 8; i++)
		for (uint8_t j = 0; j < 8; j++)
			b[fzigzag[i][j]] = mat[i][j];	
}

#ifdef USE_NEON
void quantize (int16_t *b, int16_t *q)
{
	for (uint8_t i = 0; i < 8; i++)
	{
		int16_t *bp = b + i * 8;
		int16_t *qp = q + i * 8;

		// load the integers
		int16x8_t ivb = vld1q_s16(bp);
		int16x8_t ivq = vld1q_s16(qp);
		// convert to float
		float16x8_t fvb = vcvtq_f16_s16(ivb);
		float16x8_t fvq = vcvtq_f16_s16(ivq);
		// do the quantization
		fvb = vdivq_f16(fvb, fvq);
		fvb = vrndq_f16(fvb);
		// convert to int
		ivb = vcvtq_s16_f16(fvb);
		// store the result
		vst1q_s16(bp, ivb);
	}
}
#else
void quantize(int16_t *b, int16_t *q)
{
	for (uint8_t i = 0; i < 64; i++)
		b[i] = get_round((jpegf) b[i] / q[i]);
}
#endif
