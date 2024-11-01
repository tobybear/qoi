/*
ROI:
Copyright (c) 2024, Matthew Ling


-- About

ROI is a simple byte format for storing lossless images. There are a handful
of ways it does this, at its core each pixel is diffed from the previous
pixel and stored in up to a 4 byte encoding for RGB or up to a 6 byte encoding
for RGBA.

The format has 6 ops defined:
* ROI_OP_LUMA232, ROI_OP_LUMA464, ROI_OP_LUMA777, ROI_OP_RGB: RGB ops, encode in
  1/2/3/4 bytes respectively
* ROI_OP_RUN: 1 byte RLE repeating the previous pixel 1..30 times
* ROI_OP_RGBA: 2 byte encoding used whenever alpha changes, followed by an RGB
  op to encode the RGB elements

In detail:

vr, vg, vb are red green blue diffed from the previous pixel respectively

vg_r, vg_b are vr and vb respectively diffed from vg

LUMA op values are stored with a bias, for example a 3 bit value is in the range
-4..3 inclusive, which is stored as 0..7 by adding 4

ROI_OP_RUN: xxxxx111
	1 byte op defining a run of repeating pixels, x=0..29 indicates runs of 1..30
	respectively. x=30 and x=31 is reserved for use by ROI_OP_RGB and ROI_OP_RGBA

ROI_OP_LUMA232: bbrrggg0
  1 byte op that stores vg_r and vg_b in 2 bits, vg in 3 bits

ROI_OP_LUMA464: gggggg01 bbbbrrrr
  2 byte op that stores vg_r and vg_b in 4 bits, vg in 6 bits

ROI_OP_LUMA777: ggggg011 rrrrrrgg bbbbbbbr
  3 byte op that stores vg_r, vg_b and vg in 7 bits

ROI_OP_RGB: 11110111 gggggggg rrrrrrrr bbbbbbbb
  4 byte op that stores vg_r, vg_b and vg in 8 bits, without any bias

ROI_OP_RGBA: 11111111 aaaaaaaa
	2 byte op that stores the current alpha value. Always followed by an RGB op
	to fully define a pixel

The byte stream's end is marked with 7 0x00 bytes followed a single 0x01 byte.

Unlike most qoi-like formats roi stores values within ops in little endian.
This allows for optimisations on little-endian hardware, most hardware.

*/

/* -----------------------------------------------------------------------------
Header - Public functions */

#ifndef ROI_H
#define ROI_H

#ifndef ROIDEF
#define ROIDEF static
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* A pointer to a roi_desc struct has to be supplied to all of roi's functions.
It describes either the input format (for roi_write and roi_encode), or is
filled with the description read from the file header (for roi_read and
roi_decode).

The colorspace in this roi_desc is an enum where
	0 = sRGB, i.e. gamma scaled RGB channels and a linear alpha channel
	1 = all channels are linear
You may use the constants ROI_SRGB or ROI_LINEAR. The colorspace is purely
informative. It will be saved to the file header, but does not affect
how chunks are en-/decoded. */

#define ROI_SRGB   0
#define ROI_LINEAR 1

typedef struct {
	unsigned int width;
	unsigned int height;
	unsigned char channels;
	unsigned char colorspace;
} roi_desc;

enum codepath {scalar, sse};
typedef struct{
	enum codepath path;
	unsigned char mlut;
} roi_options;

#define ROI_HEADER_SIZE 14
#define UNUSED(x) { x = x; }

#ifndef ROI_NO_STDIO

#ifndef ROI_NO_ENCODER

/* Encode raw RGB or RGBA pixels into a ROI image and write it to the file
system. The roi_desc struct must be filled with the image width, height,
number of channels (3 = RGB, 4 = RGBA) and the colorspace.

The function returns 0 on failure (invalid parameters, or fopen or malloc
failed) or the number of bytes written on success. */
ROIDEF
int roi_write(const char *filename, const void *data, const roi_desc *desc, const options *opt);

/* Encode directly from PAM/PPM file to file

The function returns 0 on failure (invalid parameters, or fopen or malloc
failed) or 1 on success. */
ROIDEF
int roi_write_from_pam(const char *pam_f, const char *roi_f, const options *opt);
ROIDEF
int roi_write_from_ppm(const char *ppm_f, const char *roi_f, const options *opt);

/* Read and decode a ROI image from the file system. If channels is 0, the
number of channels from the file header is used. If channels is 3 or 4 the
output format will be forced into this number of channels.

The function either returns NULL on failure (invalid data, or malloc or fopen
failed) or a pointer to the decoded pixels. On success, the roi_desc struct
will be filled with the description from the file header.

The returned pixel data should be ROI_FREE()d after use. */

#endif

#ifndef ROI_NO_DECODER
ROIDEF
void *roi_read(const char *filename, roi_desc *desc, int channels, const options *opt);

/* Decode directly from file to PAM/PPM file

The function returns 0 on failure (invalid parameters, or fopen or malloc
failed) or 1 on success. */
ROIDEF
int roi_read_to_pam(const char *roi_f, const char *pam_f, const options *opt);
ROIDEF
int roi_read_to_ppm(const char *roi_f, const char *ppm_f, const options *opt);
#endif

#endif /* ROI_NO_STDIO */

/* Encode raw RGB or RGBA pixels into a ROI image in memory.

The function either returns NULL on failure (invalid parameters or malloc
failed) or a pointer to the encoded data on success. On success the out_len
is set to the size in bytes of the encoded data.

The returned roi data should be ROI_FREE()d after use. */
#ifndef ROI_NO_ENCODER
ROIDEF
void *roi_encode(const void *data, const roi_desc *desc, int *out_len, const roi_options *opt);
#endif

/* Decode a ROI image from memory.

The function either returns NULL on failure (invalid parameters or malloc
failed) or a pointer to the decoded pixels. On success, the roi_desc struct
is filled with the description from the file header.

The returned pixel data should be ROI_FREE()d after use. */
#ifndef ROI_NO_DECODER
ROIDEF
void *roi_decode(const void *data, int size, roi_desc *desc, int channels);
#endif
#ifdef __cplusplus
}
#endif
#endif /* ROI_H */

/* -----------------------------------------------------------------------------
Implementation */

#ifdef ROI_IMPLEMENTATION
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#ifdef ROI_SSE
#include <immintrin.h>
#endif

#ifndef ROI_MALLOC
	#ifdef ROI_SSE
		#define ROI_MALLOC(sz) _mm_malloc(sz, 64)
		#define ROI_FREE(p)    _mm_free(p)
	#else
		#define ROI_MALLOC(sz) malloc(sz)
		#define ROI_FREE(p)    free(p)
	#endif
#endif
#ifndef ROI_ZEROARR
	#define ROI_ZEROARR(a) memset((a),0,sizeof(a))
#endif

/* 2GB is the max file size that this implementation can safely handle. We guard
against anything larger than that, assuming the worst case with 5 bytes per
pixel, rounded down to a nice clean value. 400 million pixels ought to be
enough for anybody. */
#define ROI_PIXELS_MAX ((unsigned int)400000000)

//the number of pixels to process per chunk when chunk processing
//must be a multiple of 64 for simd alignment
#define CHUNK 131072

#define ENC_READ_RGBA do{ \
	memcpy(&(s.px), s.pixels+s.px_pos, 4); \
}while(0)

typedef union {
	struct { unsigned char r, g, b, a; } rgba;
	unsigned int v;
} roi_rgba_t;

static const unsigned char roi_padding[8] = {0,0,0,0,0,0,0,1};

static void roi_write_32(unsigned char *bytes, unsigned int *p, unsigned int v) {
	bytes[(*p)++] = (0xff000000 & v) >> 24;
	bytes[(*p)++] = (0x00ff0000 & v) >> 16;
	bytes[(*p)++] = (0x0000ff00 & v) >> 8;
	bytes[(*p)++] = (0x000000ff & v);
}

static unsigned int roi_read_32(const unsigned char *bytes, unsigned int *p) {
	unsigned int a = bytes[(*p)++];
	unsigned int b = bytes[(*p)++];
	unsigned int c = bytes[(*p)++];
	unsigned int d = bytes[(*p)++];
	return a << 24 | b << 16 | c << 8 | d;
}

#define ROI_OP_LUMA232 0x00 /* xxxxxxx0 */
#define ROI_OP_LUMA464 0x01 /* xxxxxx01 */
#define ROI_OP_LUMA777 0x03 /* xxxxx011 */
#define ROI_OP_RUN     0x07 /* xxxxx111 */
#define ROI_OP_RGB     0xf7 /* 11110111 */
#define ROI_OP_RGBA    0xff /* 11111111 */

#define ROI_OP_RUN_FULL 0xef /* 11101111 */
#define ROI_RUN_FULL_VAL (30)

#define ROI_MASK_1     0x01 /* 00000001 */
#define ROI_MASK_2     0x03 /* 00000011 */
#define ROI_MASK_3     0x07 /* 00000111 */

#define ROI_MAGIC \
	(((unsigned int)'r') << 24 | ((unsigned int)'o') << 16 | \
	 ((unsigned int)'i') <<  8 | ((unsigned int)'f'))

#define EXT_STR "roi"

#define ROI_PIXEL_WORST_CASE (desc->channels==4?6:4)

#define DUMP_RUN(rrr) do{ \
	for(;rrr>=ROI_RUN_FULL_VAL;rrr-=ROI_RUN_FULL_VAL) \
		s.bytes[s.b++] = ROI_OP_RUN_FULL; \
	if (rrr) { \
		s.bytes[s.b++] = ROI_OP_RUN | ((rrr - 1)<<3); \
		rrr = 0; \
	} \
}while(0)

//	px.rgba.r = pixels[px_pos + 0];
//	px.rgba.g = pixels[px_pos + 1];
//	px.rgba.b = pixels[px_pos + 2];
#define ENC_READ_RGB do{ \
	memcpy(&(s.px), s.pixels+s.px_pos, 4); \
	s.px.v&=0x00FFFFFF; \
}while(0)

//optimised encode functions////////////////////////////////////////////////////

#define ROI_RGB_ENC_SCALAR do{\
	signed char vr = s.px.rgba.r - px_prev.rgba.r;\
	signed char vg = s.px.rgba.g - px_prev.rgba.g;\
	signed char vb = s.px.rgba.b - px_prev.rgba.b;\
	signed char vg_r = vr - vg;\
	signed char vg_b = vb - vg;\
	unsigned char ar = (vg_r<0)?(-vg_r)-1:vg_r;\
	unsigned char ag = (vg<0)?(-vg)-1:vg;\
	unsigned char ab = (vg_b<0)?(-vg_b)-1:vg_b;\
	unsigned char arb = ar|ab;\
	if ( arb < 2 && ag  < 4 ) {\
		s.bytes[s.b++]=ROI_OP_LUMA232|((vg_b+2)<<6)|((vg_r+2)<<4)|((vg+4)<<1);\
	} else if ( arb <  8 && ag  < 32 ) {\
		*(unsigned int*)(s.bytes+s.b)=ROI_OP_LUMA464|((vg_b+8)<<12)|((vg_r+8)<<8)|((vg+32)<<2); \
		s.b+=2; \
	} else if ( (arb|ag) < 64 ) {\
		*(unsigned int*)(s.bytes+s.b)=ROI_OP_LUMA777|((vg_b+64)<<17)|((vg_r+64)<<10)|((vg+64)<<3); \
		s.b+=3; \
	} else {\
		s.bytes[s.b++]=ROI_OP_RGB; \
		s.bytes[s.b++]=vg; \
		s.bytes[s.b++]=vg_r; \
		s.bytes[s.b++]=vg_b; \
	}\
}while(0)

typedef struct{
	unsigned char *bytes, *pixels;
	roi_rgba_t px;
	unsigned int b, px_pos, run, pixel_cnt;
} enc_state;

static enc_state roi_encode_chunk3_mlut(enc_state s);
static enc_state roi_encode_chunk4_mlut(enc_state s);
static enc_state roi_encode_chunk3_scalar(enc_state s);
static enc_state roi_encode_chunk4_scalar(enc_state s);
static enc_state roi_encode_chunk3_sse(enc_state s);
static enc_state roi_encode_chunk4_sse(enc_state s);

//pointers to optimised functions
#define ENC_ARR_INDEX ((opt->path<<1)|(desc->channels-3))
static enc_state (*enc_arr[])(enc_state)={
	roi_encode_chunk3_scalar, roi_encode_chunk4_scalar, roi_encode_chunk3_sse, roi_encode_chunk4_sse
};

int gen_mlut(const char *path){
	signed char vg, vg_r, vg_b;
	unsigned char ar, ag, ab, arb, *mlut, t;
	unsigned int val;
	FILE *fo;
	roi_rgba_t d={0};
	mlut = (unsigned char *)malloc(256*256*256*5);
	memset(mlut, 0, 256*256*256*5);
	for(int rr=-128;rr<128;++rr){
	for(int gg=-128;gg<128;++gg){
	for(int bb=-128;bb<128;++bb){
		d.rgba.r=rr;
		d.rgba.g=gg;
		d.rgba.b=bb;
		vg=d.rgba.g;
		vg_r = d.rgba.r - vg;
		vg_b = d.rgba.b - vg;
		ar = (vg_r<0)?(-vg_r)-1:vg_r;
		ag = (vg<0)?(-vg)-1:vg;
		ab = (vg_b<0)?(-vg_b)-1:vg_b;
		arb = ar|ab;
		if ( arb < 2 && ag  < 4 ) {
			mlut[(5*d.v)]=1;
			mlut[(5*d.v)+1]=ROI_OP_LUMA232|((vg_b+2)<<6)|((vg_r+2)<<4)|((vg+4)<<1);
		} else if ( arb <  8 && ag  < 32 ) {
			mlut[(5*d.v)]=2;
			val=ROI_OP_LUMA464|((vg_b+8)<<12)|((vg_r+8)<<8)|((vg+32)<<2);
			mlut[(5*d.v)+1]=val&255;
			mlut[(5*d.v)+2]=(val>>8)&255;
		} else if ( (arb|ag) < 64 ) {
			mlut[(5*d.v)]=3;
			val=ROI_OP_LUMA777|((vg_b+64)<<17)|((vg_r+64)<<10)|((vg+64)<<3);
			mlut[(5*d.v)+1]=val&255;
			mlut[(5*d.v)+2]=(val>>8)&255;
			mlut[(5*d.v)+3]=(val>>16)&255;
		} else {
			mlut[(5*d.v)]=4;
			mlut[(5*d.v)+1]=ROI_OP_RGB;
			t=vg;
			mlut[(5*d.v)+2]=t;
			t=vg_r;
			mlut[(5*d.v)+3]=t;
			t=vg_b;
			mlut[(5*d.v)+4]=t;
		}
	}
	}
	}
	fo=fopen(path, "wb");
	fwrite(mlut, 1, 256*256*256*5, fo);
	fclose(fo);
	free(mlut);
	return 0;
}

#ifdef ROI_MLUT_EMBED
extern char _binary_roi_mlut_start[];
unsigned char *roi_mlut=(unsigned char*)_binary_roi_mlut_start;
#else
unsigned char *roi_mlut=NULL;
#endif

static enc_state roi_encode_chunk3_mlut(enc_state s){
	roi_rgba_t px_prev=s.px, diff={0};
	unsigned int px_end=(s.pixel_cnt-1)*3;
	for (; s.px_pos <= px_end; s.px_pos += 3) {
		ENC_READ_RGB;
		while(s.px.v == px_prev.v) {
			++s.run;
			if(s.px_pos == px_end){
				for(;s.run>=ROI_RUN_FULL_VAL;s.run-=ROI_RUN_FULL_VAL)
					s.bytes[s.b++] = ROI_OP_RUN_FULL;
				s.px_pos+=3;
				return s;
			}
			s.px_pos+=3;
			ENC_READ_RGB;
		}
		DUMP_RUN(s.run);
		diff.rgba.r=s.px.rgba.r-px_prev.rgba.r;
		diff.rgba.g=s.px.rgba.g-px_prev.rgba.g;
		diff.rgba.b=s.px.rgba.b-px_prev.rgba.b;
		*(unsigned int*)(s.bytes+s.b)=*(unsigned int*)(roi_mlut+(diff.v*5)+1);
		s.b+=roi_mlut[diff.v*5];
		px_prev = s.px;
	}
	return s;
}

static enc_state roi_encode_chunk4_mlut(enc_state s){
	roi_rgba_t px_prev=s.px, diff={0};
	unsigned int px_end=(s.pixel_cnt-1)*4;
	for (; s.px_pos <= px_end; s.px_pos += 4) {
		ENC_READ_RGBA;
		while(s.px.v == px_prev.v) {
			++s.run;
			if(s.px_pos == px_end) {
				for(;s.run>=ROI_RUN_FULL_VAL;s.run-=ROI_RUN_FULL_VAL)
					s.bytes[s.b++] = ROI_OP_RUN_FULL;
				s.px_pos+=4;
				return s;
			}
			s.px_pos+=4;
			ENC_READ_RGBA;
		}
		DUMP_RUN(s.run);
		if(s.px.rgba.a!=px_prev.rgba.a){
			s.bytes[s.b++] = ROI_OP_RGBA;
			s.bytes[s.b++] = s.px.rgba.a;
		}
		diff.rgba.r=s.px.rgba.r-px_prev.rgba.r;
		diff.rgba.g=s.px.rgba.g-px_prev.rgba.g;
		diff.rgba.b=s.px.rgba.b-px_prev.rgba.b;
		*(unsigned int*)(s.bytes+s.b)=*(unsigned int*)(roi_mlut+(diff.v*5)+1);
		s.b+=roi_mlut[diff.v*5];
		px_prev = s.px;
	}
	return s;
}

static enc_state roi_encode_chunk3_scalar(enc_state s){
	roi_rgba_t px_prev=s.px;
	unsigned int px_end=(s.pixel_cnt-1)*3;
	for (; s.px_pos <= px_end; s.px_pos += 3) {
		ENC_READ_RGB;
		while(s.px.v == px_prev.v) {
			++s.run;
			if(s.px_pos == px_end){
				for(;s.run>=ROI_RUN_FULL_VAL;s.run-=ROI_RUN_FULL_VAL)
					s.bytes[s.b++] = ROI_OP_RUN_FULL;
				s.px_pos+=3;
				return s;
			}
			s.px_pos+=3;
			ENC_READ_RGB;
		}
		DUMP_RUN(s.run);
		ROI_RGB_ENC_SCALAR;
		px_prev = s.px;
	}
	return s;
}

static enc_state roi_encode_chunk4_scalar(enc_state s){
	roi_rgba_t px_prev=s.px;
	unsigned int px_end=(s.pixel_cnt-1)*4;
	for (; s.px_pos <= px_end; s.px_pos += 4) {
		ENC_READ_RGBA;
		while(s.px.v == px_prev.v) {
			++s.run;
			if(s.px_pos == px_end) {
				for(;s.run>=ROI_RUN_FULL_VAL;s.run-=ROI_RUN_FULL_VAL)
					s.bytes[s.b++] = ROI_OP_RUN_FULL;
				s.px_pos+=4;
				return s;
			}
			s.px_pos+=4;
			ENC_READ_RGBA;
		}
		DUMP_RUN(s.run);
		if(s.px.rgba.a!=px_prev.rgba.a){
			s.bytes[s.b++] = ROI_OP_RGBA;
			s.bytes[s.b++] = s.px.rgba.a;
		}
		ROI_RGB_ENC_SCALAR;
		px_prev = s.px;
	}
	return s;
}

#ifdef ROI_SSE
//load the next 16 bytes, diff pixels
#define LOAD16(raw, diff, prev, offset, lshift, rshift) do{ \
	raw=_mm_loadu_si128((__m128i const*)(s.pixels+s.px_pos+offset)); \
	diff=_mm_slli_si128(raw, lshift); \
	prev=_mm_srli_si128(prev, rshift); \
	diff=_mm_or_si128(diff, prev); \
	diff=_mm_sub_epi8(raw, diff); \
}while(0)

//de-interleave one plane from 3 vectors containing RGB
#define PLANAR_SHUFFLE(plane, source1, source2, source3, shufflemask) do{ \
	plane=_mm_blendv_epi8(source1, source2, blend1); \
	plane=_mm_blendv_epi8(plane, source3, blend2); \
	plane=_mm_shuffle_epi8(plane, shufflemask); \
}while(0)

//do (x<0)?(-x)-1:x for a single plane
#define ABSOLUTER(plane, absolute) do{ \
	w2=_mm_cmpgt_epi8(_mm_setzero_si128(), plane); \
	w1=_mm_and_si128(w2, plane); \
	w1=_mm_add_epi8(w1, _mm_set1_epi8(1)); \
	w1=_mm_abs_epi8(w1); \
	absolute=_mm_blendv_epi8(plane, w1, w2); \
}while(0)

//the following 2 macros:
// normalise value depending on opcode
// shift value to where it is in the op
// combine into 4 result vectors
#define NORMALISE_SHIFT16_EMBIGGEN(plane, opmask, value, shift) do{ \
	w1=_mm_add_epi8(plane, value); \
	w1=_mm_and_si128(w1, opmask); \
	w2=_mm_unpacklo_epi8(w1, _mm_setzero_si128()); \
	w2=_mm_slli_epi16(w2, shift); \
	w3=_mm_unpacklo_epi16(w2, _mm_setzero_si128()); \
	res0=_mm_or_si128(w3, res0); \
	w3=_mm_unpackhi_epi16(w2, _mm_setzero_si128()); \
	res1=_mm_or_si128(w3, res1); \
	w2=_mm_unpackhi_epi8(w1, _mm_setzero_si128()); \
	w2=_mm_slli_epi16(w2, shift); \
	w3=_mm_unpacklo_epi16(w2, _mm_setzero_si128()); \
	res2=_mm_or_si128(w3, res2); \
	w3=_mm_unpackhi_epi16(w2, _mm_setzero_si128()); \
	res3=_mm_or_si128(w3, res3); \
}while(0)

#define NORMALISE_SHIFT32_EMBIGGEN(plane, opmask, value, shift) do{ \
	w1=_mm_add_epi8(plane, value); \
	w1=_mm_and_si128(w1, opmask); \
	w2=_mm_unpacklo_epi8(w1, _mm_setzero_si128()); \
	w3=_mm_unpacklo_epi16(w2, _mm_setzero_si128()); \
	w3=_mm_slli_epi32(w3, shift); \
	res0=_mm_or_si128(w3, res0); \
	w3=_mm_unpackhi_epi16(w2, _mm_setzero_si128()); \
	w3=_mm_slli_epi32(w3, shift); \
	res1=_mm_or_si128(w3, res1); \
	w2=_mm_unpackhi_epi8(w1, _mm_setzero_si128()); \
	w3=_mm_unpacklo_epi16(w2, _mm_setzero_si128()); \
	w3=_mm_slli_epi32(w3, shift); \
	res2=_mm_or_si128(w3, res2); \
	w3=_mm_unpackhi_epi16(w2, _mm_setzero_si128()); \
	w3=_mm_slli_epi32(w3, shift); \
	res3=_mm_or_si128(w3, res3); \
}while(0)

#define NORMALISE_INPLACE(plane, opmask, value) do{ \
	w1=_mm_add_epi8(plane, value); \
	plane=_mm_blendv_epi8(plane, w1, opmask); \
}while(0)

//sse lut
static const uint8_t writer_lut[4096] = {//shuffle used bytes in output vector to the left ready for writing
	0,4,8,12,0,0,0,0,0,0,0,0,0,0,0,0, 0,1,4,8,12,0,0,0,0,0,0,0,0,0,0,0, 0,1,2,4,8,12,0,0,0,0,0,0,0,0,0,0, 0,1,2,3,4,8,12,0,0,0,0,0,0,0,0,0,
	0,4,5,8,12,0,0,0,0,0,0,0,0,0,0,0, 0,1,4,5,8,12,0,0,0,0,0,0,0,0,0,0, 0,1,2,4,5,8,12,0,0,0,0,0,0,0,0,0, 0,1,2,3,4,5,8,12,0,0,0,0,0,0,0,0,
	0,4,5,6,8,12,0,0,0,0,0,0,0,0,0,0, 0,1,4,5,6,8,12,0,0,0,0,0,0,0,0,0, 0,1,2,4,5,6,8,12,0,0,0,0,0,0,0,0, 0,1,2,3,4,5,6,8,12,0,0,0,0,0,0,0,
	0,4,5,6,7,8,12,0,0,0,0,0,0,0,0,0, 0,1,4,5,6,7,8,12,0,0,0,0,0,0,0,0, 0,1,2,4,5,6,7,8,12,0,0,0,0,0,0,0, 0,1,2,3,4,5,6,7,8,12,0,0,0,0,0,0,
	0,4,8,9,12,0,0,0,0,0,0,0,0,0,0,0, 0,1,4,8,9,12,0,0,0,0,0,0,0,0,0,0, 0,1,2,4,8,9,12,0,0,0,0,0,0,0,0,0, 0,1,2,3,4,8,9,12,0,0,0,0,0,0,0,0,
	0,4,5,8,9,12,0,0,0,0,0,0,0,0,0,0, 0,1,4,5,8,9,12,0,0,0,0,0,0,0,0,0, 0,1,2,4,5,8,9,12,0,0,0,0,0,0,0,0, 0,1,2,3,4,5,8,9,12,0,0,0,0,0,0,0,
	0,4,5,6,8,9,12,0,0,0,0,0,0,0,0,0, 0,1,4,5,6,8,9,12,0,0,0,0,0,0,0,0, 0,1,2,4,5,6,8,9,12,0,0,0,0,0,0,0, 0,1,2,3,4,5,6,8,9,12,0,0,0,0,0,0,
	0,4,5,6,7,8,9,12,0,0,0,0,0,0,0,0, 0,1,4,5,6,7,8,9,12,0,0,0,0,0,0,0, 0,1,2,4,5,6,7,8,9,12,0,0,0,0,0,0, 0,1,2,3,4,5,6,7,8,9,12,0,0,0,0,0,
	0,4,8,9,10,12,0,0,0,0,0,0,0,0,0,0, 0,1,4,8,9,10,12,0,0,0,0,0,0,0,0,0, 0,1,2,4,8,9,10,12,0,0,0,0,0,0,0,0, 0,1,2,3,4,8,9,10,12,0,0,0,0,0,0,0,
	0,4,5,8,9,10,12,0,0,0,0,0,0,0,0,0, 0,1,4,5,8,9,10,12,0,0,0,0,0,0,0,0, 0,1,2,4,5,8,9,10,12,0,0,0,0,0,0,0, 0,1,2,3,4,5,8,9,10,12,0,0,0,0,0,0,
	0,4,5,6,8,9,10,12,0,0,0,0,0,0,0,0, 0,1,4,5,6,8,9,10,12,0,0,0,0,0,0,0, 0,1,2,4,5,6,8,9,10,12,0,0,0,0,0,0, 0,1,2,3,4,5,6,8,9,10,12,0,0,0,0,0,
	0,4,5,6,7,8,9,10,12,0,0,0,0,0,0,0, 0,1,4,5,6,7,8,9,10,12,0,0,0,0,0,0, 0,1,2,4,5,6,7,8,9,10,12,0,0,0,0,0, 0,1,2,3,4,5,6,7,8,9,10,12,0,0,0,0,
	0,4,8,9,10,11,12,0,0,0,0,0,0,0,0,0, 0,1,4,8,9,10,11,12,0,0,0,0,0,0,0,0, 0,1,2,4,8,9,10,11,12,0,0,0,0,0,0,0, 0,1,2,3,4,8,9,10,11,12,0,0,0,0,0,0,
	0,4,5,8,9,10,11,12,0,0,0,0,0,0,0,0, 0,1,4,5,8,9,10,11,12,0,0,0,0,0,0,0, 0,1,2,4,5,8,9,10,11,12,0,0,0,0,0,0, 0,1,2,3,4,5,8,9,10,11,12,0,0,0,0,0,
	0,4,5,6,8,9,10,11,12,0,0,0,0,0,0,0, 0,1,4,5,6,8,9,10,11,12,0,0,0,0,0,0, 0,1,2,4,5,6,8,9,10,11,12,0,0,0,0,0, 0,1,2,3,4,5,6,8,9,10,11,12,0,0,0,0,
	0,4,5,6,7,8,9,10,11,12,0,0,0,0,0,0, 0,1,4,5,6,7,8,9,10,11,12,0,0,0,0,0, 0,1,2,4,5,6,7,8,9,10,11,12,0,0,0,0, 0,1,2,3,4,5,6,7,8,9,10,11,12,0,0,0,
	0,4,8,12,13,0,0,0,0,0,0,0,0,0,0,0, 0,1,4,8,12,13,0,0,0,0,0,0,0,0,0,0, 0,1,2,4,8,12,13,0,0,0,0,0,0,0,0,0, 0,1,2,3,4,8,12,13,0,0,0,0,0,0,0,0,
	0,4,5,8,12,13,0,0,0,0,0,0,0,0,0,0, 0,1,4,5,8,12,13,0,0,0,0,0,0,0,0,0, 0,1,2,4,5,8,12,13,0,0,0,0,0,0,0,0, 0,1,2,3,4,5,8,12,13,0,0,0,0,0,0,0,
	0,4,5,6,8,12,13,0,0,0,0,0,0,0,0,0, 0,1,4,5,6,8,12,13,0,0,0,0,0,0,0,0, 0,1,2,4,5,6,8,12,13,0,0,0,0,0,0,0, 0,1,2,3,4,5,6,8,12,13,0,0,0,0,0,0,
	0,4,5,6,7,8,12,13,0,0,0,0,0,0,0,0, 0,1,4,5,6,7,8,12,13,0,0,0,0,0,0,0, 0,1,2,4,5,6,7,8,12,13,0,0,0,0,0,0, 0,1,2,3,4,5,6,7,8,12,13,0,0,0,0,0,
	0,4,8,9,12,13,0,0,0,0,0,0,0,0,0,0, 0,1,4,8,9,12,13,0,0,0,0,0,0,0,0,0, 0,1,2,4,8,9,12,13,0,0,0,0,0,0,0,0, 0,1,2,3,4,8,9,12,13,0,0,0,0,0,0,0,
	0,4,5,8,9,12,13,0,0,0,0,0,0,0,0,0, 0,1,4,5,8,9,12,13,0,0,0,0,0,0,0,0, 0,1,2,4,5,8,9,12,13,0,0,0,0,0,0,0, 0,1,2,3,4,5,8,9,12,13,0,0,0,0,0,0,
	0,4,5,6,8,9,12,13,0,0,0,0,0,0,0,0, 0,1,4,5,6,8,9,12,13,0,0,0,0,0,0,0, 0,1,2,4,5,6,8,9,12,13,0,0,0,0,0,0, 0,1,2,3,4,5,6,8,9,12,13,0,0,0,0,0,
	0,4,5,6,7,8,9,12,13,0,0,0,0,0,0,0, 0,1,4,5,6,7,8,9,12,13,0,0,0,0,0,0, 0,1,2,4,5,6,7,8,9,12,13,0,0,0,0,0, 0,1,2,3,4,5,6,7,8,9,12,13,0,0,0,0,
	0,4,8,9,10,12,13,0,0,0,0,0,0,0,0,0, 0,1,4,8,9,10,12,13,0,0,0,0,0,0,0,0, 0,1,2,4,8,9,10,12,13,0,0,0,0,0,0,0, 0,1,2,3,4,8,9,10,12,13,0,0,0,0,0,0,
	0,4,5,8,9,10,12,13,0,0,0,0,0,0,0,0, 0,1,4,5,8,9,10,12,13,0,0,0,0,0,0,0, 0,1,2,4,5,8,9,10,12,13,0,0,0,0,0,0, 0,1,2,3,4,5,8,9,10,12,13,0,0,0,0,0,
	0,4,5,6,8,9,10,12,13,0,0,0,0,0,0,0, 0,1,4,5,6,8,9,10,12,13,0,0,0,0,0,0, 0,1,2,4,5,6,8,9,10,12,13,0,0,0,0,0, 0,1,2,3,4,5,6,8,9,10,12,13,0,0,0,0,
	0,4,5,6,7,8,9,10,12,13,0,0,0,0,0,0, 0,1,4,5,6,7,8,9,10,12,13,0,0,0,0,0, 0,1,2,4,5,6,7,8,9,10,12,13,0,0,0,0, 0,1,2,3,4,5,6,7,8,9,10,12,13,0,0,0,
	0,4,8,9,10,11,12,13,0,0,0,0,0,0,0,0, 0,1,4,8,9,10,11,12,13,0,0,0,0,0,0,0, 0,1,2,4,8,9,10,11,12,13,0,0,0,0,0,0, 0,1,2,3,4,8,9,10,11,12,13,0,0,0,0,0,
	0,4,5,8,9,10,11,12,13,0,0,0,0,0,0,0, 0,1,4,5,8,9,10,11,12,13,0,0,0,0,0,0, 0,1,2,4,5,8,9,10,11,12,13,0,0,0,0,0, 0,1,2,3,4,5,8,9,10,11,12,13,0,0,0,0,
	0,4,5,6,8,9,10,11,12,13,0,0,0,0,0,0, 0,1,4,5,6,8,9,10,11,12,13,0,0,0,0,0, 0,1,2,4,5,6,8,9,10,11,12,13,0,0,0,0, 0,1,2,3,4,5,6,8,9,10,11,12,13,0,0,0,
	0,4,5,6,7,8,9,10,11,12,13,0,0,0,0,0, 0,1,4,5,6,7,8,9,10,11,12,13,0,0,0,0, 0,1,2,4,5,6,7,8,9,10,11,12,13,0,0,0, 0,1,2,3,4,5,6,7,8,9,10,11,12,13,0,0,
	0,4,8,12,13,14,0,0,0,0,0,0,0,0,0,0, 0,1,4,8,12,13,14,0,0,0,0,0,0,0,0,0, 0,1,2,4,8,12,13,14,0,0,0,0,0,0,0,0, 0,1,2,3,4,8,12,13,14,0,0,0,0,0,0,0,
	0,4,5,8,12,13,14,0,0,0,0,0,0,0,0,0, 0,1,4,5,8,12,13,14,0,0,0,0,0,0,0,0, 0,1,2,4,5,8,12,13,14,0,0,0,0,0,0,0, 0,1,2,3,4,5,8,12,13,14,0,0,0,0,0,0,
	0,4,5,6,8,12,13,14,0,0,0,0,0,0,0,0, 0,1,4,5,6,8,12,13,14,0,0,0,0,0,0,0, 0,1,2,4,5,6,8,12,13,14,0,0,0,0,0,0, 0,1,2,3,4,5,6,8,12,13,14,0,0,0,0,0,
	0,4,5,6,7,8,12,13,14,0,0,0,0,0,0,0, 0,1,4,5,6,7,8,12,13,14,0,0,0,0,0,0, 0,1,2,4,5,6,7,8,12,13,14,0,0,0,0,0, 0,1,2,3,4,5,6,7,8,12,13,14,0,0,0,0,
	0,4,8,9,12,13,14,0,0,0,0,0,0,0,0,0, 0,1,4,8,9,12,13,14,0,0,0,0,0,0,0,0, 0,1,2,4,8,9,12,13,14,0,0,0,0,0,0,0, 0,1,2,3,4,8,9,12,13,14,0,0,0,0,0,0,
	0,4,5,8,9,12,13,14,0,0,0,0,0,0,0,0, 0,1,4,5,8,9,12,13,14,0,0,0,0,0,0,0, 0,1,2,4,5,8,9,12,13,14,0,0,0,0,0,0, 0,1,2,3,4,5,8,9,12,13,14,0,0,0,0,0,
	0,4,5,6,8,9,12,13,14,0,0,0,0,0,0,0, 0,1,4,5,6,8,9,12,13,14,0,0,0,0,0,0, 0,1,2,4,5,6,8,9,12,13,14,0,0,0,0,0, 0,1,2,3,4,5,6,8,9,12,13,14,0,0,0,0,
	0,4,5,6,7,8,9,12,13,14,0,0,0,0,0,0, 0,1,4,5,6,7,8,9,12,13,14,0,0,0,0,0, 0,1,2,4,5,6,7,8,9,12,13,14,0,0,0,0, 0,1,2,3,4,5,6,7,8,9,12,13,14,0,0,0,
	0,4,8,9,10,12,13,14,0,0,0,0,0,0,0,0, 0,1,4,8,9,10,12,13,14,0,0,0,0,0,0,0, 0,1,2,4,8,9,10,12,13,14,0,0,0,0,0,0, 0,1,2,3,4,8,9,10,12,13,14,0,0,0,0,0,
	0,4,5,8,9,10,12,13,14,0,0,0,0,0,0,0, 0,1,4,5,8,9,10,12,13,14,0,0,0,0,0,0, 0,1,2,4,5,8,9,10,12,13,14,0,0,0,0,0, 0,1,2,3,4,5,8,9,10,12,13,14,0,0,0,0,
	0,4,5,6,8,9,10,12,13,14,0,0,0,0,0,0, 0,1,4,5,6,8,9,10,12,13,14,0,0,0,0,0, 0,1,2,4,5,6,8,9,10,12,13,14,0,0,0,0, 0,1,2,3,4,5,6,8,9,10,12,13,14,0,0,0,
	0,4,5,6,7,8,9,10,12,13,14,0,0,0,0,0, 0,1,4,5,6,7,8,9,10,12,13,14,0,0,0,0, 0,1,2,4,5,6,7,8,9,10,12,13,14,0,0,0, 0,1,2,3,4,5,6,7,8,9,10,12,13,14,0,0,
	0,4,8,9,10,11,12,13,14,0,0,0,0,0,0,0, 0,1,4,8,9,10,11,12,13,14,0,0,0,0,0,0, 0,1,2,4,8,9,10,11,12,13,14,0,0,0,0,0, 0,1,2,3,4,8,9,10,11,12,13,14,0,0,0,0,
	0,4,5,8,9,10,11,12,13,14,0,0,0,0,0,0, 0,1,4,5,8,9,10,11,12,13,14,0,0,0,0,0, 0,1,2,4,5,8,9,10,11,12,13,14,0,0,0,0, 0,1,2,3,4,5,8,9,10,11,12,13,14,0,0,0,
	0,4,5,6,8,9,10,11,12,13,14,0,0,0,0,0, 0,1,4,5,6,8,9,10,11,12,13,14,0,0,0,0, 0,1,2,4,5,6,8,9,10,11,12,13,14,0,0,0, 0,1,2,3,4,5,6,8,9,10,11,12,13,14,0,0,
	0,4,5,6,7,8,9,10,11,12,13,14,0,0,0,0, 0,1,4,5,6,7,8,9,10,11,12,13,14,0,0,0, 0,1,2,4,5,6,7,8,9,10,11,12,13,14,0,0, 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,0,
	0,4,8,12,13,14,15,0,0,0,0,0,0,0,0,0, 0,1,4,8,12,13,14,15,0,0,0,0,0,0,0,0, 0,1,2,4,8,12,13,14,15,0,0,0,0,0,0,0, 0,1,2,3,4,8,12,13,14,15,0,0,0,0,0,0,
	0,4,5,8,12,13,14,15,0,0,0,0,0,0,0,0, 0,1,4,5,8,12,13,14,15,0,0,0,0,0,0,0, 0,1,2,4,5,8,12,13,14,15,0,0,0,0,0,0, 0,1,2,3,4,5,8,12,13,14,15,0,0,0,0,0,
	0,4,5,6,8,12,13,14,15,0,0,0,0,0,0,0, 0,1,4,5,6,8,12,13,14,15,0,0,0,0,0,0, 0,1,2,4,5,6,8,12,13,14,15,0,0,0,0,0, 0,1,2,3,4,5,6,8,12,13,14,15,0,0,0,0,
	0,4,5,6,7,8,12,13,14,15,0,0,0,0,0,0, 0,1,4,5,6,7,8,12,13,14,15,0,0,0,0,0, 0,1,2,4,5,6,7,8,12,13,14,15,0,0,0,0, 0,1,2,3,4,5,6,7,8,12,13,14,15,0,0,0,
	0,4,8,9,12,13,14,15,0,0,0,0,0,0,0,0, 0,1,4,8,9,12,13,14,15,0,0,0,0,0,0,0, 0,1,2,4,8,9,12,13,14,15,0,0,0,0,0,0, 0,1,2,3,4,8,9,12,13,14,15,0,0,0,0,0,
	0,4,5,8,9,12,13,14,15,0,0,0,0,0,0,0, 0,1,4,5,8,9,12,13,14,15,0,0,0,0,0,0, 0,1,2,4,5,8,9,12,13,14,15,0,0,0,0,0, 0,1,2,3,4,5,8,9,12,13,14,15,0,0,0,0,
	0,4,5,6,8,9,12,13,14,15,0,0,0,0,0,0, 0,1,4,5,6,8,9,12,13,14,15,0,0,0,0,0, 0,1,2,4,5,6,8,9,12,13,14,15,0,0,0,0, 0,1,2,3,4,5,6,8,9,12,13,14,15,0,0,0,
	0,4,5,6,7,8,9,12,13,14,15,0,0,0,0,0, 0,1,4,5,6,7,8,9,12,13,14,15,0,0,0,0, 0,1,2,4,5,6,7,8,9,12,13,14,15,0,0,0, 0,1,2,3,4,5,6,7,8,9,12,13,14,15,0,0,
	0,4,8,9,10,12,13,14,15,0,0,0,0,0,0,0, 0,1,4,8,9,10,12,13,14,15,0,0,0,0,0,0, 0,1,2,4,8,9,10,12,13,14,15,0,0,0,0,0, 0,1,2,3,4,8,9,10,12,13,14,15,0,0,0,0,
	0,4,5,8,9,10,12,13,14,15,0,0,0,0,0,0, 0,1,4,5,8,9,10,12,13,14,15,0,0,0,0,0, 0,1,2,4,5,8,9,10,12,13,14,15,0,0,0,0, 0,1,2,3,4,5,8,9,10,12,13,14,15,0,0,0,
	0,4,5,6,8,9,10,12,13,14,15,0,0,0,0,0, 0,1,4,5,6,8,9,10,12,13,14,15,0,0,0,0, 0,1,2,4,5,6,8,9,10,12,13,14,15,0,0,0, 0,1,2,3,4,5,6,8,9,10,12,13,14,15,0,0,
	0,4,5,6,7,8,9,10,12,13,14,15,0,0,0,0, 0,1,4,5,6,7,8,9,10,12,13,14,15,0,0,0, 0,1,2,4,5,6,7,8,9,10,12,13,14,15,0,0, 0,1,2,3,4,5,6,7,8,9,10,12,13,14,15,0,
	0,4,8,9,10,11,12,13,14,15,0,0,0,0,0,0, 0,1,4,8,9,10,11,12,13,14,15,0,0,0,0,0, 0,1,2,4,8,9,10,11,12,13,14,15,0,0,0,0, 0,1,2,3,4,8,9,10,11,12,13,14,15,0,0,0,
	0,4,5,8,9,10,11,12,13,14,15,0,0,0,0,0, 0,1,4,5,8,9,10,11,12,13,14,15,0,0,0,0, 0,1,2,4,5,8,9,10,11,12,13,14,15,0,0,0, 0,1,2,3,4,5,8,9,10,11,12,13,14,15,0,0,
	0,4,5,6,8,9,10,11,12,13,14,15,0,0,0,0, 0,1,4,5,6,8,9,10,11,12,13,14,15,0,0,0, 0,1,2,4,5,6,8,9,10,11,12,13,14,15,0,0, 0,1,2,3,4,5,6,8,9,10,11,12,13,14,15,0,
	0,4,5,6,7,8,9,10,11,12,13,14,15,0,0,0, 0,1,4,5,6,7,8,9,10,11,12,13,14,15,0,0, 0,1,2,4,5,6,7,8,9,10,11,12,13,14,15,0, 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
};
static const uint8_t writer_len[256] = {
	4, 5, 6, 7, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 10, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 10, 8, 9, 10, 11,
	6, 7, 8, 9, 7, 8, 9, 10, 8, 9, 10, 11, 9, 10, 11, 12,	7, 8, 9, 10, 8, 9, 10, 11, 9, 10, 11, 12, 10, 11, 12, 13,
	5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 10, 8, 9, 10, 11, 6, 7, 8, 9, 7, 8, 9, 10, 8, 9, 10, 11, 9, 10, 11, 12,
	7, 8, 9, 10, 8, 9, 10, 11, 9, 10, 11, 12, 10, 11, 12, 13, 8, 9, 10, 11, 9, 10, 11, 12, 10, 11, 12, 13, 11, 12, 13, 14,
	6, 7, 8, 9, 7, 8, 9, 10, 8, 9, 10, 11, 9, 10, 11, 12, 7, 8, 9, 10, 8, 9, 10, 11, 9, 10, 11, 12, 10, 11, 12, 13,
	8, 9, 10, 11, 9, 10, 11, 12, 10, 11, 12, 13, 11, 12, 13, 14, 9, 10, 11, 12, 10, 11, 12, 13, 11, 12, 13, 14, 12, 13, 14, 15,
	7, 8, 9, 10, 8, 9, 10, 11, 9, 10, 11, 12, 10, 11, 12, 13, 8, 9, 10, 11, 9, 10, 11, 12, 10, 11, 12, 13, 11, 12, 13, 14,
	9, 10, 11, 12, 10, 11, 12, 13, 11, 12, 13, 14, 12, 13, 14, 15, 10, 11, 12, 13, 11, 12, 13, 14, 12, 13, 14, 15, 13, 14, 15, 16
};

//SSE to the point that we have 4 output vectors ready to write
#define SSE_COMMON do{ \
	/*convert vr, vb to vg_r, vg_b respectively*/ \
	r=_mm_sub_epi8(r, g); \
	b=_mm_sub_epi8(b, g); \
	/*generate absolute vectors for each of r, g, b, (vg<0)?(-vg)-1:vg;*/ \
	ABSOLUTER(r, ar); \
	ABSOLUTER(g, ag); \
	ABSOLUTER(b, ab); \
	/*determine how to store pixels*/ \
	/* 1 byte if arb<2, ag<4*/ \
	/* 2 byte if arb<8, ag<32*/ \
	/* 3 byte if argb<64*/ \
	/* 4 byte otherwise*/ \
	arb=_mm_or_si128(ar, ab); \
	op1=_mm_subs_epu8(ag, _mm_set1_epi8(2)); \
	op1=_mm_or_si128(op1, arb); \
	op1=_mm_cmpgt_epi8(_mm_set1_epi8(2), op1);/*op1*/ \
	op2=_mm_subs_epu8(ag, _mm_set1_epi8(24)); \
	op2=_mm_or_si128(op2, arb); \
	op2=_mm_cmpgt_epi8(_mm_set1_epi8(8), op2);/*op1|op2*/ \
	op3=_mm_cmpgt_epi8(_mm_set1_epi8(64), _mm_or_si128(arb, ag));/*op1|op2|op3*/ \
	op4=_mm_andnot_si128(op3, _mm_set1_epi8(-1));/*op4*/ \
	op3=_mm_sub_epi8(op3, op2);/*op3*/ \
	op2=_mm_sub_epi8(op2, op1);/*op2*/ \
	res0=_mm_setzero_si128(); \
	res1=_mm_setzero_si128(); \
	res2=_mm_setzero_si128(); \
	res3=_mm_setzero_si128(); \
	/*build opcode vector*/ \
	opuse=_mm_and_si128(op2, _mm_set1_epi8(1)); \
	opuse=_mm_or_si128(opuse, _mm_and_si128(op3, _mm_set1_epi8(3))); \
	opuse=_mm_or_si128(opuse, _mm_and_si128(op4, _mm_set1_epi8(-9))); \
	/*apply opcodes to output*/ \
	w1=_mm_unpacklo_epi8(opuse, _mm_setzero_si128()); \
	w2=_mm_unpacklo_epi16(w1, _mm_setzero_si128()); \
	res0=_mm_or_si128(w2, res0); \
	w2=_mm_unpackhi_epi16(w1, _mm_setzero_si128()); \
	res1=_mm_or_si128(w2, res1); \
	w1=_mm_unpackhi_epi8(opuse, _mm_setzero_si128()); \
	w2=_mm_unpacklo_epi16(w1, _mm_setzero_si128()); \
	res2=_mm_or_si128(w2, res2); \
	w2=_mm_unpackhi_epi16(w1, _mm_setzero_si128()); \
	res3=_mm_or_si128(w2, res3); \
	/*bbrrggg0*/ \
	NORMALISE_SHIFT16_EMBIGGEN(g, op1, _mm_set1_epi8(4), 1); \
	NORMALISE_SHIFT16_EMBIGGEN(r, op1, _mm_set1_epi8(2), 4); \
	NORMALISE_SHIFT16_EMBIGGEN(b, op1, _mm_set1_epi8(2), 6); \
	/*bbbbrrrr gggggg01*/ \
	NORMALISE_SHIFT16_EMBIGGEN(g, op2, _mm_set1_epi8(32), 2); \
	NORMALISE_SHIFT16_EMBIGGEN(r, op2, _mm_set1_epi8(8), 8); \
	NORMALISE_SHIFT16_EMBIGGEN(b, op2, _mm_set1_epi8(8), 12); \
	/*bbbbbbbr rrrrrrgg ggggg011*/ \
	NORMALISE_SHIFT16_EMBIGGEN(g, op3, _mm_set1_epi8(64), 3); \
	NORMALISE_SHIFT32_EMBIGGEN(r, op3, _mm_set1_epi8(64), 10); \
	NORMALISE_SHIFT32_EMBIGGEN(b, op3, _mm_set1_epi8(64), 17); \
	/*bbbbbbbb rrrrrrrr gggggggg 11110111*/ \
	/*shift op4 g*/ \
	w1=_mm_and_si128(g, op4); \
	w2=_mm_unpacklo_epi8(_mm_setzero_si128(), w1);/*switched to end up at 2nd byte position*/ \
	w3=_mm_unpacklo_epi16(w2, _mm_setzero_si128()); \
	res0=_mm_or_si128(w3, res0); \
	w3=_mm_unpackhi_epi16(w2, _mm_setzero_si128()); \
	res1=_mm_or_si128(w3, res1); \
	w2=_mm_unpackhi_epi8(_mm_setzero_si128(), w1);/*switched*/ \
	w3=_mm_unpacklo_epi16(w2, _mm_setzero_si128()); \
	res2=_mm_or_si128(w3, res2); \
	w3=_mm_unpackhi_epi16(w2, _mm_setzero_si128()); \
	res3=_mm_or_si128(w3, res3); \
	/*shift op4 r*/ \
	w1=_mm_and_si128(r, op4); \
	w2=_mm_unpacklo_epi8(w1, _mm_setzero_si128()); \
	w3=_mm_unpacklo_epi16(_mm_setzero_si128(), w2);/*switch*/ \
	res0=_mm_or_si128(w3, res0); \
	w3=_mm_unpackhi_epi16(_mm_setzero_si128(), w2);/*switch*/ \
	res1=_mm_or_si128(w3, res1); \
	w2=_mm_unpackhi_epi8(w1, _mm_setzero_si128()); \
	w3=_mm_unpacklo_epi16(_mm_setzero_si128(), w2);/*switch*/ \
	res2=_mm_or_si128(w3, res2); \
	w3=_mm_unpackhi_epi16(_mm_setzero_si128(), w2);/*switch*/ \
	res3=_mm_or_si128(w3, res3); \
	/*shift op4 b*/ \
	w1=_mm_and_si128(b, op4); \
	w2=_mm_unpacklo_epi8(_mm_setzero_si128(), w1);/*switch*/ \
	w3=_mm_unpacklo_epi16(_mm_setzero_si128(), w2);/*switch*/ \
	res0=_mm_or_si128(w3, res0); \
	w3=_mm_unpackhi_epi16(_mm_setzero_si128(), w2);/*switch*/ \
	res1=_mm_or_si128(w3, res1); \
	w2=_mm_unpackhi_epi8(_mm_setzero_si128(), w1);/*switch*/ \
	w3=_mm_unpacklo_epi16(_mm_setzero_si128(), w2);/*switch*/ \
	res2=_mm_or_si128(w3, res2); \
	w3=_mm_unpackhi_epi16(_mm_setzero_si128(), w2);/*switch*/ \
	res3=_mm_or_si128(w3, res3); \
}while(0)

//write output vector scalar, RLE bytes are involved
#define ROI_SSE_WRITE_SAFE(data, len) do{ \
	_mm_storeu_si128((__m128i*)dump, data); \
	if(dump[0]==0xA8) \
		++s.run; \
	else{ \
		DUMP_RUN(s.run); \
		memcpy(s.bytes+s.b, dump+0, len[0]); \
		s.b+=len[0]; \
	} \
	if(dump[4]==0xA8) \
		++s.run; \
	else{ \
		DUMP_RUN(s.run); \
		memcpy(s.bytes+s.b, dump+4, len[1]); \
		s.b+=len[1]; \
	} \
	if(dump[8]==0xA8) \
		++s.run; \
	else{ \
		DUMP_RUN(s.run); \
		memcpy(s.bytes+s.b, dump+8, len[2]); \
		s.b+=len[2]; \
	} \
	if(dump[12]==0xA8) \
		++s.run; \
	else{ \
		DUMP_RUN(s.run); \
		memcpy(s.bytes+s.b, dump+12, len[3]); \
		s.b+=len[3]; \
	} \
}while(0)

//write output vector with SSE, RLE bytes not involved
//DUMP_RUN(s.run) not here, handled elsewhere because multiple paths
#define ROI_SSE_WRITE_QUICK(data, lut_offset) do{ \
	/*write first vec*/ \
	w1=_mm_loadu_si128((__m128i const*)(writer_lut)+((lut_index>>lut_offset)&255)); \
	w1=_mm_shuffle_epi8(data, w1); \
	_mm_storeu_si128((__m128i*)(s.bytes+s.b), w1); \
	s.b+=writer_len[(lut_index>>lut_offset)&255]; \
}while(0)

#define ROI_SSE_GETLUT0 do{ \
	/*get lut for first 8 pixels*/ \
	w2=_mm_unpacklo_epi8(op2, _mm_setzero_si128()); \
	w1=_mm_unpacklo_epi8(_mm_setzero_si128(), op3); \
	w2=_mm_or_si128(w2, w1); \
	w1=_mm_unpacklo_epi8(op4, op4); \
	w2=_mm_or_si128(w2, w1); \
	lut_index=_mm_movemask_epi8(w2); \
}while(0)

#define ROI_SSE_GETLUT1 do{ \
	/*get lut for next 8 pixels*/ \
	w2=_mm_unpackhi_epi8(op2, _mm_setzero_si128()); \
	w1=_mm_unpackhi_epi8(_mm_setzero_si128(), op3); \
	w2=_mm_or_si128(w2, w1); \
	w1=_mm_unpackhi_epi8(op4, op4); \
	w2=_mm_or_si128(w2, w1); \
	lut_index=_mm_movemask_epi8(w2); \
}while(0)

//RLE not involved, write all 4 output vectors with SSE branchless
#define ROI_SSE_QUICK do{ \
	DUMP_RUN(s.run); \
	ROI_SSE_GETLUT0; \
	ROI_SSE_WRITE_QUICK(res0, 0); \
	ROI_SSE_WRITE_QUICK(res1, 8); \
	ROI_SSE_GETLUT1; \
	ROI_SSE_WRITE_QUICK(res2, 0); \
	ROI_SSE_WRITE_QUICK(res3, 8); \
}while(0)

//RLE involved, branch to do scalar/sse output as necessary per vector
#define SSE_CAREFUL do{ \
	/*get op lengths*/ \
	opuse=_mm_and_si128(op1, _mm_set1_epi8(1)); \
	opuse=_mm_or_si128(opuse, _mm_and_si128(op2, _mm_set1_epi8(2))); \
	opuse=_mm_or_si128(opuse, _mm_and_si128(op3, _mm_set1_epi8(3))); \
	opuse=_mm_or_si128(opuse, _mm_and_si128(op4, _mm_set1_epi8(4))); \
	_mm_storeu_si128((__m128i*)opuse_, opuse); \
	if(!rle_status[0]||!rle_status[1]){ \
		ROI_SSE_GETLUT0; \
	} \
	if(!rle_status[0]){ \
		DUMP_RUN(s.run); \
		ROI_SSE_WRITE_QUICK(res0, 0); \
	} \
	else{ \
		ROI_SSE_WRITE_SAFE(res0, (opuse_+0)); \
	} \
	if(!rle_status[1]){ \
		DUMP_RUN(s.run); \
		ROI_SSE_WRITE_QUICK(res1, 8); \
	} \
	else{ \
		ROI_SSE_WRITE_SAFE(res1, (opuse_+4)); \
	} \
	if(!rle_status[2]||!rle_status[3]){ \
		ROI_SSE_GETLUT1; \
	} \
	if(!rle_status[2]){ \
		DUMP_RUN(s.run); \
		ROI_SSE_WRITE_QUICK(res2, 0); \
	} \
	else{ \
		ROI_SSE_WRITE_SAFE(res2, (opuse_+8)); \
	} \
	if(!rle_status[3]){ \
		DUMP_RUN(s.run); \
		ROI_SSE_WRITE_QUICK(res3, 8); \
	} \
	else{ \
		ROI_SSE_WRITE_SAFE(res3, (opuse_+12)); \
	} \
}while(0)

static enc_state roi_encode_chunk4_sse(enc_state s){
	__m128i ia, ib, ic, id, da, db, dc, dd, r, g, b, a, ar, ag, ab, arb, w1, w2, w3, w4, w5, w6, previous;
	__m128i gshuf, shuf1, shuf2, blend;
	__m128i op1, op2, op3, op4, opuse, res0, res1, res2, res3;
	int lut_index=0;
	unsigned char dump[16], opuse_[16];
	unsigned int rle_status[4];

	//constants
	shuf1=_mm_setr_epi8(0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15);
	shuf2=_mm_setr_epi8(1,5,9,13,0,4,8,12,3,7,11,15,2,6,10,14);
	gshuf=_mm_setr_epi8(8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7);
	blend=_mm_setr_epi8(0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1);

	id=_mm_setr_epi8(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, s.px.rgba.r, s.px.rgba.g, s.px.rgba.b, s.px.rgba.a);//prev pixel
	for (; s.px_pos < s.pixel_cnt*4; s.px_pos += 64) {
		//load pixels
		previous=_mm_and_si128(id, id);
		LOAD16(ia, da, id,  0, 4, 12);
		LOAD16(ib, db, ia, 16, 4, 12);
		LOAD16(ic, dc, ib, 32, 4, 12);
		LOAD16(id, dd, ic, 48, 4, 12);

		if(_mm_test_all_zeros( _mm_or_si128(_mm_or_si128(da, db), _mm_or_si128(dc, dd)), _mm_set1_epi8(-1))){//all RLE
			s.run+=16;
			continue;
		}

		//unpack into rgba planes
		w1=_mm_shuffle_epi8(da, shuf1);//r4g4b4a4
		w2=_mm_shuffle_epi8(db, shuf1);//r4g4b4a4
		w3=_mm_shuffle_epi8(dc, shuf2);//g4r4a4b4
		w4=_mm_shuffle_epi8(dd, shuf2);//g4r4a4b4
		w5=_mm_unpackhi_epi32(w1, w2);//b8a8
		w6=_mm_unpackhi_epi32(w3, w4);//a8b8
		a=_mm_blendv_epi8(w6, w5, blend);//out of order, irrelevant
		if(!_mm_test_all_zeros(a, _mm_set1_epi8(-1))){//alpha present, scalar this iteration
			//TODO ditch this scalar code, make a parallel alpha op vector and zip together with rgb vector
			unsigned int pixel_cnt_store=s.pixel_cnt;
			_mm_storeu_si128((__m128i*)dump, previous);
			s.px.rgba.r=dump[12];
			s.px.rgba.g=dump[13];
			s.px.rgba.b=dump[14];
			s.px.rgba.a=dump[15];
			s.pixel_cnt=(s.px_pos/4)+16;
			s=enc_arr[1](s);
			s.px_pos-=64;
			s.pixel_cnt=pixel_cnt_store;
			continue;
		}
		//no alpha, finish extracting planes then re-use rgb sse implementation
		b=_mm_blendv_epi8(w5, w6, blend);
		w1=_mm_unpacklo_epi32(w1, w2);//r8g8
		w2=_mm_unpacklo_epi32(w3, w4);//g8r8
		r=_mm_blendv_epi8(w1, w2, blend);
		g=_mm_blendv_epi8(w2, w1, blend);//out of order
		g=_mm_shuffle_epi8(g, gshuf);//in order

		SSE_COMMON;
		w1=_mm_cmpeq_epi8(_mm_or_si128(r, _mm_or_si128(g, b)), _mm_set1_epi8(0));
		if(!_mm_testz_si128(w1, w1)){//RLE, tread lightly
			_mm_storeu_si128((__m128i*)rle_status, w1);
			SSE_CAREFUL;
		}
		else{
			ROI_SSE_QUICK;
		}
	}
	_mm_storeu_si128((__m128i*)dump, id);
	s.px.rgba.r=dump[12];
	s.px.rgba.g=dump[13];
	s.px.rgba.b=dump[14];
	s.px.rgba.a=dump[15];
	return s;
}

static enc_state roi_encode_chunk3_sse(enc_state s){
	__m128i aa, bb, cc, da, db, dc, r, g, b, ar, ag, ab, arb, w1, w2, w3;
	__m128i rshuf, gshuf, bshuf, blend1, blend2;
	__m128i op1, op2, op3, op4, opuse, res0, res1, res2, res3;
	int lut_index=0;
	unsigned char dump[16], opuse_[16];
	unsigned int rle_status[4];

	//constants
	rshuf=_mm_setr_epi8(0,3,6,9,12,15, 2,5,8,11,14, 1,4,7,10,13);
	gshuf=_mm_setr_epi8(1,4,7,10,13, 0,3,6,9,12,15, 2,5,8,11,14);
	bshuf=_mm_setr_epi8(2,5,8,11,14, 1,4,7,10,13, 0,3,6,9,12,15);
	blend1=_mm_setr_epi8(0,0,-1,0,0,-1,0,0,-1,0,0,-1,0,0,-1,0);
	blend2=_mm_setr_epi8(0,-1,0,0,-1,0,0,-1,0,0,-1,0,0,-1,0,0);

	cc=_mm_setr_epi8(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, s.px.rgba.r, s.px.rgba.g, s.px.rgba.b);//prev pixel
	for (; s.px_pos < s.pixel_cnt*3; s.px_pos += 48) {
		//load and diff next 16 pixels
		LOAD16(aa, da, cc, 0, 3, 13);
		LOAD16(bb, db, aa, 16, 3, 13);
		LOAD16(cc, dc, bb, 32, 3, 13);

		if(_mm_test_all_zeros( _mm_or_si128(da, _mm_or_si128(db, dc)), _mm_set1_epi8(-1))){//all RLE
			s.run+=16;
			continue;
		}

		/*convert to rgb vectors*/
		PLANAR_SHUFFLE(r, da, db, dc, rshuf);
		PLANAR_SHUFFLE(g, db, dc, da, gshuf);
		PLANAR_SHUFFLE(b, dc, da, db, bshuf);

		SSE_COMMON;
		w1=_mm_cmpeq_epi8(_mm_or_si128(r, _mm_or_si128(g, b)), _mm_set1_epi8(0));
		if(!_mm_testz_si128(w1, w1)){//RLE, tread lightly
			_mm_storeu_si128((__m128i*)rle_status, w1);
			SSE_CAREFUL;
		}
		else{
			ROI_SSE_QUICK;
		}
	}
	_mm_storeu_si128((__m128i*)dump, cc);
	s.px.rgba.r=dump[13];
	s.px.rgba.g=dump[14];
	s.px.rgba.b=dump[15];
	return s;
}
#else
//not compiled with ROI_SSE, replace sse functions with scalar placeholders
static enc_state roi_encode_chunk3_sse(enc_state s){
	return roi_encode_chunk3_scalar(s);
}
static enc_state roi_encode_chunk4_sse(enc_state s){
	return roi_encode_chunk4_scalar(s);
}
#endif

//Optimised decode functions////////////////////////////////////////////////////
typedef struct{
	unsigned char *bytes, *pixels;
	roi_rgba_t px;
	unsigned int b, b_limit, b_present, p, p_limit, px_pos, run, pixel_cnt, pixel_curr;
} dec_state;

#define ROI_DECODE_COMMON \
	;int b1 = s.bytes[s.b++]; \
	if ((b1 & ROI_MASK_1) == ROI_OP_LUMA232) { \
		int vg = ((b1>>1)&7) - 6; \
		s.px.rgba.r += vg + ((b1 >> 4) & 3); \
		s.px.rgba.g += vg + 2; \
		s.px.rgba.b += vg + ((b1 >> 6) & 3); \
	} \
	else if ((b1 & ROI_MASK_2) == ROI_OP_LUMA464) { \
		int b2=s.bytes[s.b++]; \
		int vg = ((b1>>2)&63) - 40; \
		s.px.rgba.r += vg + ((b2     ) & 0x0f); \
		s.px.rgba.g += vg + 8; \
		s.px.rgba.b += vg + ((b2 >>4) & 0x0f); \
	} \
	else if ((b1 & ROI_MASK_3) == ROI_OP_LUMA777) { \
		int b2=s.bytes[s.b++]; \
		int b3=s.bytes[s.b++]; \
		int vg = (((b2&3)<<5)|((b1>>3)&31))-128; \
		s.px.rgba.r += vg + (((b3&1)<<6)|((b2>>2)&63)); \
		s.px.rgba.g += vg + 64; \
		s.px.rgba.b += vg + ((b3>>1)&127); \
	} \
	else if (b1 == ROI_OP_RGB) { \
		signed char vg=s.bytes[s.b++]; \
		signed char b3=s.bytes[s.b++]; \
		signed char b4=s.bytes[s.b++]; \
		s.px.rgba.r += vg + b3; \
		s.px.rgba.g += vg; \
		s.px.rgba.b += vg + b4; \
	}

static dec_state dec_in4out4(dec_state s){
	while( ((s.b+6)<s.b_present) && ((s.px_pos+4)<=s.p_limit) && (s.pixel_cnt!=s.pixel_curr) ){
		if (s.run)
			s.run--;
		else{
			OP_RGBA_GOTO:
			ROI_DECODE_COMMON
			else if (b1 == ROI_OP_RGBA) {
				s.px.rgba.a = s.bytes[s.b++];
				goto OP_RGBA_GOTO;
			}
			else// if ((b1 & ROI_MASK_3) == ROI_OP_RUN)
				s.run = ((b1>>3) & 0x1f);
		}
		s.pixels[s.px_pos + 0] = s.px.rgba.r;
		s.pixels[s.px_pos + 1] = s.px.rgba.g;
		s.pixels[s.px_pos + 2] = s.px.rgba.b;
		s.pixels[s.px_pos + 3] = s.px.rgba.a;
		s.px_pos+=4;
		s.pixel_curr++;
	}
	return s;
}

static dec_state dec_in4out3(dec_state s){
	while( ((s.b+6)<s.b_present) && ((s.px_pos+3)<=s.p_limit) && (s.pixel_cnt!=s.pixel_curr) ){
		if (s.run)
			s.run--;
		else{
			OP_RGBA_GOTO:
			ROI_DECODE_COMMON
			else if (b1 == ROI_OP_RGBA) {
				s.px.rgba.a = s.bytes[s.b++];
				goto OP_RGBA_GOTO;
			}
			else// if ((b1 & ROI_MASK_3) == ROI_OP_RUN)
				s.run = ((b1>>3) & 0x1f);
		}
		s.pixels[s.px_pos + 0] = s.px.rgba.r;
		s.pixels[s.px_pos + 1] = s.px.rgba.g;
		s.pixels[s.px_pos + 2] = s.px.rgba.b;
		s.px_pos+=3;
		s.pixel_curr++;
	}
	return s;
}

static dec_state dec_in3out4(dec_state s){
	while( ((s.b+6)<s.b_present) && ((s.px_pos+4)<=s.p_limit) && (s.pixel_cnt!=s.pixel_curr) ){
		if (s.run)
			s.run--;
		else{
			ROI_DECODE_COMMON
			else// if ((b1 & ROI_MASK_3) == ROI_OP_RUN)
				s.run = ((b1>>3) & 0x1f);
		}
		s.pixels[s.px_pos + 0] = s.px.rgba.r;
		s.pixels[s.px_pos + 1] = s.px.rgba.g;
		s.pixels[s.px_pos + 2] = s.px.rgba.b;
		s.pixels[s.px_pos + 3] = s.px.rgba.a;
		s.px_pos+=4;
		s.pixel_curr++;
	}
	return s;
}

static dec_state dec_in3out3(dec_state s){
	while( ((s.b+6)<s.b_present) && ((s.px_pos+3)<=s.p_limit) && (s.pixel_cnt!=s.pixel_curr) ){
		if (s.run)
			s.run--;
		else{
			ROI_DECODE_COMMON
			else// if ((b1 & ROI_MASK_3) == ROI_OP_RUN)
				s.run = ((b1>>3) & 0x1f);
		}
		s.pixels[s.px_pos + 0] = s.px.rgba.r;
		s.pixels[s.px_pos + 1] = s.px.rgba.g;
		s.pixels[s.px_pos + 2] = s.px.rgba.b;
		s.px_pos+=3;
		s.pixel_curr++;
	}
	return s;
}

#define DEC_ARR_INDEX (((desc->channels-3)<<1)|(channels-3))
static dec_state (*dec_arr[])(dec_state)={dec_in3out3, dec_in3out4, dec_in4out3, dec_in4out4};

static void roi_encode_init(const roi_desc *desc, unsigned char *bytes, unsigned int *p, roi_rgba_t *px_prev) {
	roi_write_32(bytes, p, ROI_MAGIC);
	roi_write_32(bytes, p, desc->width);
	roi_write_32(bytes, p, desc->height);
	bytes[(*p)++] = desc->channels;
	bytes[(*p)++] = desc->colorspace;
	px_prev->rgba.r = 0;
	px_prev->rgba.g = 0;
	px_prev->rgba.b = 0;
#ifdef ROI
	px_prev->rgba.a = desc->channels==3?0:255;//simplify ENC_READ_RGB
#else
	px_prev->rgba.a = 255;
#endif
}

#ifndef ROI_NO_ENCODER
ROIDEF
void *roi_encode(const void *data, const roi_desc *desc, int *out_len, const roi_options *opt) {
	enc_state s={0};
	int i, max_size;

	if (
		data == NULL || out_len == NULL || desc == NULL ||
		desc->width == 0 || desc->height == 0 ||
		desc->channels < 3 || desc->channels > 4 ||
		desc->colorspace > 1 ||
		desc->height >= ROI_PIXELS_MAX / desc->width ||
		opt->path>1
	)
		return NULL;
#ifdef ROI
	if(opt->mlut){
		enc_arr[0]=roi_encode_chunk3_mlut;
		enc_arr[1]=roi_encode_chunk4_mlut;
	}
#endif

	max_size =
		desc->width * desc->height * ROI_PIXEL_WORST_CASE +
		ROI_HEADER_SIZE + sizeof(roi_padding);

	if(!(s.bytes = (unsigned char *) ROI_MALLOC(max_size)))
		return NULL;
	s.pixels=(unsigned char *)data;

	roi_encode_init(desc, s.bytes, &(s.b), &(s.px));
	if((desc->width * desc->height)/CHUNK){//encode most of the input as the largest multiple of chunk size for simd
		s.pixel_cnt=(desc->width * desc->height)-((desc->width * desc->height)%CHUNK);
		s=enc_arr[ENC_ARR_INDEX](s);
	}
	if((desc->width * desc->height)%CHUNK){//encode the trailing input scalar
		s.pixel_cnt=(desc->width * desc->height);
		s=enc_arr[desc->channels-3](s);
	}
	DUMP_RUN(s.run);
	for (i = 0; i < (int)sizeof(roi_padding); i++)
		s.bytes[s.b++] = roi_padding[i];
	*out_len = s.b;
	return s.bytes;
}
#endif

#ifndef ROI_NO_DECODER
ROIDEF
void *roi_decode(const void *data, int size, roi_desc *desc, int channels) {
	unsigned int header_magic;
	dec_state s={0};

	if (
		data == NULL || desc == NULL ||
		(channels != 0 && channels != 3 && channels != 4) ||
		size < ROI_HEADER_SIZE + (int)sizeof(roi_padding)
	)
		return NULL;

	s.bytes=(unsigned char*)data;

	header_magic = roi_read_32(s.bytes, &(s.b));
	desc->width = roi_read_32(s.bytes, &(s.b));
	desc->height = roi_read_32(s.bytes, &(s.b));
	desc->channels = s.bytes[s.b++];
	desc->colorspace = s.bytes[s.b++];

	if (
		desc->width == 0 || desc->height == 0 ||
		desc->channels < 3 || desc->channels > 4 ||
		desc->colorspace > 1 ||
		header_magic != ROI_MAGIC ||
		desc->height >= ROI_PIXELS_MAX / desc->width
	)
		return NULL;

	if (channels == 0)
		channels = desc->channels;

	s.pixel_cnt=desc->width * desc->height;
	s.p_limit=s.pixel_cnt*channels;
	if(!(s.pixels = (unsigned char *)ROI_MALLOC(s.p_limit)))
		return NULL;
	s.b_limit=size;
	s.b_present=size;
	s.px.rgba.a=255;

	s=dec_arr[DEC_ARR_INDEX](s);

	return s.pixels;
}
#endif

#ifndef ROI_NO_STDIO
#include <stdio.h>

static inline FILE* roi_fopen(const char *path, const char *mode){
	if(0==strcmp(path, "-"))
		return *mode=='r'?stdin:stdout;
	else
		return fopen(path, mode);
}

static inline void roi_fclose(const char *path, FILE *stream){
	if(0!=strcmp(path, "-"))
		fclose(stream);
}

//decode to a format that contains raw pixels in RGB/A
static int roi_read_to_file(FILE *fi, const char *out_f, char *head, size_t head_len, roi_desc *desc, int channels, const options *opt){
	dec_state s={0};
	FILE *fo;
	UNUSED(opt);

	if(
		desc->width==0 || desc->height==0 ||
		desc->channels<3 || desc->channels>4 ||
		desc->colorspace>1
	)
		goto BADEXIT0;

	if(!(fo=roi_fopen(out_f, "wb")))
		goto BADEXIT0;

	if(head_len){
		if(head_len!=fwrite(head, 1, head_len, fo))
			goto BADEXIT1;
	}

	s.b_limit=CHUNK*(desc->channels==3?2:3);
	if(!(s.bytes = (unsigned char *)ROI_MALLOC(s.b_limit)))
		goto BADEXIT1;
	s.p_limit=CHUNK*channels;
	if(!(s.pixels = (unsigned char *)ROI_MALLOC(s.p_limit)))
		goto BADEXIT2;
	s.px.rgba.a=255;
	s.pixel_cnt=desc->width*desc->height;
	while(s.pixel_curr!=s.pixel_cnt){
		s.b_present+=fread(s.bytes+s.b_present, 1, s.b_limit-s.b_present, fi);
		s=dec_arr[DEC_ARR_INDEX](s);
		if(!s.px_pos)//truncated input
			goto BADEXIT3;
		if(s.px_pos!=fwrite(s.pixels, 1, s.px_pos, fo))
			goto BADEXIT3;
		memmove(s.bytes, s.bytes+s.b, s.b_present-s.b);
		s.b_present-=s.b;
		s.b=0;
		s.px_pos=0;
	}

	ROI_FREE(s.pixels);
	ROI_FREE(s.bytes);
	roi_fclose(out_f, fo);
	return 0;
	BADEXIT3:
	ROI_FREE(s.pixels);
	BADEXIT2:
	ROI_FREE(s.bytes);
	BADEXIT1:
	roi_fclose(out_f, fo);
	BADEXIT0:
	return 1;
}

static int file_to_desc(FILE *fi, roi_desc *desc){
	unsigned char head[14];
	if(14!=fread(head, 1, 14, fi))
		return 1;
	if(ROI_MAGIC!=(head[0] << 24 | head[1] << 16 | head[2] << 8 | head[3]))
		return 1;
	desc->width = head[4] << 24 | head[5] << 16 | head[6] << 8 | head[7];
	desc->height = head[8] << 24 | head[9] << 16 | head[10] << 8 | head[11];
	desc->channels = head[12];
	desc->colorspace = head[13];
	return 0;
}

ROIDEF
int roi_read_to_pam(const char *roi_f, const char *pam_f, const options *opt) {
	char head[128];
	FILE *fi;
	roi_desc desc;
	if(!(fi=roi_fopen(roi_f, "rb")))
		goto BADEXIT0;
	if(file_to_desc(fi, &desc))
		goto BADEXIT1;

	sprintf(head, "P7\nWIDTH %u\nHEIGHT %u\nDEPTH %u\nMAXVAL 255\nTUPLTYPE RGB%s\nENDHDR\n", desc.width, desc.height, desc.channels, desc.channels==3?"":"_ALPHA");

	if(roi_read_to_file(fi, pam_f, head, strlen(head), &desc, desc.channels, opt))
		goto BADEXIT1;

	roi_fclose(roi_f, fi);
	return 0;
	BADEXIT1:
	roi_fclose(roi_f, fi);
	BADEXIT0:
	return 1;
}

ROIDEF
int roi_read_to_ppm(const char *roi_f, const char *ppm_f, const options *opt) {
	char head[128];
	FILE *fi;
	roi_desc desc;
	if(!(fi=roi_fopen(roi_f, "rb")))
		goto BADEXIT0;
	if(file_to_desc(fi, &desc))
		goto BADEXIT1;

	sprintf(head, "P6 %u %u 255\n", desc.width, desc.height);
	if(roi_read_to_file(fi, ppm_f, head, strlen(head), &desc, 3, opt))
		goto BADEXIT1;

	roi_fclose(roi_f, fi);
	return 0;
	BADEXIT1:
	roi_fclose(roi_f, fi);
	BADEXIT0:
	return 1;
}

//process from an opened raw file directly
static inline int roi_write_from_file(FILE *fi, const char *roi_f, roi_desc *desc, const options *opt){
	enc_state s={0};
	FILE *fo;
	unsigned int i, totpixels;

	if(opt->path>1)
		goto BADEXIT0;
	if(!(fo=roi_fopen(roi_f, "wb")))
		goto BADEXIT0;

	if(!(s.pixels = (unsigned char *)ROI_MALLOC((CHUNK*desc->channels)+1)))
		goto BADEXIT1;
	if(!(s.bytes = (unsigned char *)ROI_MALLOC(CHUNK*ROI_PIXEL_WORST_CASE)))
		goto BADEXIT2;

	roi_encode_init(desc, s.bytes, &(s.b), &(s.px));
	if(s.b!=fwrite(s.bytes, 1, s.b, fo))
		goto BADEXIT3;

#ifdef ROI
	if(opt->mlut){
		enc_arr[0]=roi_encode_chunk3_mlut;
		enc_arr[1]=roi_encode_chunk4_mlut;
	}
#endif

	totpixels=desc->width*desc->height;
	s.pixel_cnt=CHUNK;
	for(i=0;(i+CHUNK)<=totpixels;i+=CHUNK){
		if((CHUNK*desc->channels)!=fread(s.pixels, 1, CHUNK*desc->channels, fi))
			goto BADEXIT3;
		s.b=0;
		s.px_pos=0;
		s=enc_arr[ENC_ARR_INDEX](s);
		if(s.b!=fwrite(s.bytes, 1, s.b, fo))
			goto BADEXIT3;
	}
	if(i<totpixels){//finish scalar
		if(((totpixels-i)*desc->channels)!=fread(s.pixels, 1, (totpixels-i)*desc->channels, fi))
			goto BADEXIT3;
		s.b=0;
		s.px_pos=0;
		s.pixel_cnt=totpixels-i;
		s=enc_arr[desc->channels-3](s);
		if(s.b!=fwrite(s.bytes, 1, s.b, fo))
			goto BADEXIT3;
	}
	s.b=0;
	DUMP_RUN(s.run);
	if(s.b && s.b!=fwrite(s.bytes, 1, s.b, fo))
		goto BADEXIT3;
	if(sizeof(roi_padding)!=fwrite(roi_padding, 1, sizeof(roi_padding), fo))
		goto BADEXIT3;

	ROI_FREE(s.bytes);
	ROI_FREE(s.pixels);
	roi_fclose(roi_f, fo);
	return 0;
	BADEXIT3:
	ROI_FREE(s.bytes);
	BADEXIT2:
	ROI_FREE(s.pixels);
	BADEXIT1:
	roi_fclose(roi_f, fo);
	BADEXIT0:
	return 1;
}

//PAM and PPM reading macros
#define roi_isspace(num) (num==' '||((num>=0x09) && (num<=0x0d)))
#define roi_isdigit(num) ((num>='0') && (num<='9'))

#define PAM_READ1 do{ \
	if(1!=fread(&t, 1, 1, fi)) \
		goto BADEXIT1; \
}while(0)

#define PAM_SPACE_NUM(var) do{ \
	if(!roi_isspace(t)) \
		goto BADEXIT1; \
	do { \
		PAM_READ1; \
	} while(roi_isspace(t)); \
	if(!roi_isdigit(t)) \
		goto BADEXIT1; \
	while(roi_isdigit(t)){ \
		var*=10; \
		var+=(t-'0'); \
		PAM_READ1; \
	} \
}while(0);

#define PAM_EXPECT(val) do{ \
	PAM_READ1; \
	if(t!=val) \
		goto BADEXIT1; \
}while(0)

#define PAM_COMMENT do{ \
	while(t!='\n'){ \
		PAM_READ1; \
	} \
}while(0)

#ifndef ROI_NO_ENCODER
ROIDEF
int roi_write_from_pam(const char *pam_f, const char *roi_f, const options *opt) {
	roi_desc desc;
	char *token[]={"WIDTH", "HEIGHT", "DEPTH", "MAXVAL", "ENDHDR\n"};
	unsigned int hval[4]={0};
	unsigned char t;
	unsigned int i, j;
	FILE *fi;

	if(!(fi=roi_fopen(pam_f, "rb")))
		goto BADEXIT0;

	PAM_EXPECT('P');
	PAM_EXPECT('7');
	PAM_EXPECT('\n');

	while(1){//read header line by line
		PAM_READ1;
		if(t=='\n')//empty line
			continue;
		if(t=='#'){//comment
			PAM_COMMENT;
			continue;
		}
		for(i=0;i<5;++i){
			if(t==token[i][0])
				break;
		}
		if(i==5){//irrelevant token
			PAM_COMMENT;
			continue;
		}
		for(j=1;token[i][j];++j){
			PAM_READ1;
			if(t!=token[i][j])
				break;
		}
		if(token[i][j]){
			PAM_COMMENT;
			continue;
		}
		if(i==4)//ENDHDR
			break;
		//WIDTH HEIGHT DEPTH MAXVAL
		if(hval[i])//there can be only one
			goto BADEXIT1;
		PAM_READ1;
		PAM_SPACE_NUM(hval[i]);
	}
	if(hval[0]==0 || hval[1]==0 || hval[2]<3 || hval[2]>4 || hval[3]>255 )
		goto BADEXIT1;
	desc.width=hval[0];
	desc.height=hval[1];
	desc.channels=hval[2];
	desc.colorspace=0;

	if(roi_write_from_file(fi, roi_f, &desc, opt))
		goto BADEXIT1;

	roi_fclose(pam_f, fi);
	return 0;
	BADEXIT1:
	roi_fclose(pam_f, fi);
	BADEXIT0:
	return 1;
}

ROIDEF
int roi_write_from_ppm(const char *ppm_f, const char *roi_f, const options *opt) {
	roi_desc desc={0};
	unsigned char t;
	unsigned int maxval=0;
	FILE *fi;

	if(!(fi=roi_fopen(ppm_f, "rb")))
		goto BADEXIT0;

	PAM_EXPECT('P');
	PAM_EXPECT('6');
	PAM_READ1;
	PAM_SPACE_NUM(desc.width);
	PAM_SPACE_NUM(desc.height);
	PAM_SPACE_NUM(maxval);
	if(t=='#'){
		PAM_COMMENT;
	}
	if(!roi_isspace(t))
		goto BADEXIT1;
	if(maxval>255)
		goto BADEXIT1;
	desc.channels=3;
	desc.colorspace=0;
	if(roi_write_from_file(fi, roi_f, &desc, opt))
		goto BADEXIT1;

	roi_fclose(ppm_f, fi);
	return 0;
	BADEXIT1:
	roi_fclose(ppm_f, fi);
	BADEXIT0:
	return 1;
}

ROIDEF
int roi_write(const char *filename, const void *data, const roi_desc *desc, const options *opt) {
	FILE *f = fopen(filename, "wb");
	int size, err;
	void *encoded;

	if (!f)
		return 0;

	encoded = roi_encode(data, desc, &size, opt);
	if (!encoded) {
		fclose(f);
		return 0;
	}

	fwrite(encoded, 1, size, f);
	fflush(f);
	err = ferror(f);
	fclose(f);

	ROI_FREE(encoded);
	return err ? 0 : size;
}
#endif

#ifndef ROI_NO_DECODER
ROIDEF
void *roi_read(const char *filename, roi_desc *desc, int channels, const options *opt) {
	FILE *f = fopen(filename, "rb");
	int size, bytes_read;
	void *pixels, *data;
	UNUSED(opt);

	if (!f)
		return NULL;

	fseek(f, 0, SEEK_END);
	size = ftell(f);
	if (size <= 0 || fseek(f, 0, SEEK_SET) != 0) {
		fclose(f);
		return NULL;
	}

	if (!(data = ROI_MALLOC(size))) {
		fclose(f);
		return NULL;
	}

	bytes_read = fread(data, 1, size, f);
	fclose(f);
	pixels = (bytes_read != size) ? NULL : roi_decode(data, bytes_read, desc, channels);
	ROI_FREE(data);
	return pixels;
}
#endif

#endif /* ROI_NO_STDIO */
#endif /* ROI_IMPLEMENTATION */
