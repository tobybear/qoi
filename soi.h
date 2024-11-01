/*

Copyright (c) 2021, Dominic Szablewski - https://phoboslab.org
SPDX-License-Identifier: MIT


SOI - The "Special OK Image" format for fast, lossless image compression

-- About

SOI encodes and decodes images in a lossless format. Compared to stb_image and
stb_image_write SOI offers 20x-50x faster encoding, 3x-4x faster decoding and
20% better compression.


-- Synopsis

// Define `SOI_IMPLEMENTATION` in *one* C/C++ file before including this
// library to create the implementation.

#define SOI_IMPLEMENTATION
#include "soi.h"

// Encode and store an RGBA buffer to the file system. The soi_desc describes
// the input pixel data.
soi_write("image_new.soi", rgba_pixels, &(soi_desc){
	.width = 1920,
	.height = 1080,
	.channels = 4,
	.colorspace = SOI_SRGB
});

// Load and decode a SOI image from the file system into a 32bbp RGBA buffer.
// The soi_desc struct will be filled with the width, height, number of channels
// and colorspace read from the file header.
soi_desc desc;
void *rgba_pixels = soi_read("image.soi", &desc, 4);



-- Documentation

This library provides the following functions;
- soi_read    -- read and decode a SOI file
- soi_decode  -- decode the raw bytes of a SOI image from memory
- soi_write   -- encode and write a SOI file
- soi_encode  -- encode an rgba buffer into a SOI image in memory

See the function declaration below for the signature and more information.

If you don't want/need the soi_read and soi_write functions, you can define
SOI_NO_STDIO before including this library.

This library uses malloc() and free(). To supply your own malloc implementation
you can define SOI_MALLOC and SOI_FREE before including this library.


-- Data Format

A SOI file has a 14 byte header, followed by any number of data "chunks" and an
8-byte end marker.

struct soi_header_t {
	char     magic[4];   // magic bytes "soif"
	uint32_t width;      // image width in pixels (BE)
	uint32_t height;     // image height in pixels (BE)
	uint8_t  channels;   // 3 = RGB, 4 = RGBA
	uint8_t  colorspace; // 0 = sRGB with linear alpha, 1 = all channels linear
};

Images are encoded row by row, left to right, top to bottom. The decoder and
encoder start with {r: 0, g: 0, b: 0, a: 255} as the previous pixel value. An
image is complete when all pixels specified by width * height have been covered.

Pixels are encoded as
 - a run of the previous pixel
 - a difference to the previous pixel value in r,g,b in a 1/2/3 byte encoding
 - full r,g,b or r,g,b,a values

The color channels are assumed to not be premultiplied with the alpha channel
("un-premultiplied alpha").

Each chunk starts with a 1, 2, 3 or 8-bit tag, followed by a number of data bits.
The bit length of chunks is divisible by 8 - i.e. all chunks are byte aligned. All
values encoded in these data bits have the most significant bit on the left.

The 8-bit tags have precedence over the other tags. A decoder must check for the
presence of an 8-bit tag first.

The byte stream's end is marked with 7 0x00 bytes followed a single 0x01 byte.


The possible chunks are:


.- SOI_OP_LUMA232 --------.
|         Byte[0]         |
|  7  6  5  4  3  2  1  0 |
|---+------+-------+------|
|  0| dr-dg|  dg   | db-dg|
`-------------------------`
1-bit tag b0
2-bit   red channel difference minus green channel difference -2..1
3-bit green channel difference from the previous pixel between -4..3
2-bit  blue channel difference minus green channel difference -2..1

The green channel is used to indicate the general direction of change and is
encoded in 3 bits. The red and blue channels (dr and db) base their diffs off
of the green channel difference and are encoded in 2 bits. I.e.:
	dr_dg = (cur_px.r - prev_px.r) - (cur_px.g - prev_px.g)
	db_dg = (cur_px.b - prev_px.b) - (cur_px.g - prev_px.g)

The difference to the current channel values are using a wraparound operation,
so "10 - 13" will result in 253, while "250 + 7" will result in 1.

Values are stored as unsigned integers with a bias of 4 for the green channel
and a bias of 2 for the red and blue channel.

The alpha value remains unchanged from the previous pixel.


.- SOI_OP_LUMA464 ----------------------------------.
|         Byte[0]         |         Byte[1]         |
|  7  6  5  4  3  2  1  0 |  7  6  5  4  3  2  1  0 |
|-------+-----------------+-------------+-----------|
|  1  0 |  green diff     |   dr - dg   |  db - dg  |
`---------------------------------------------------`
2-bit tag b10
6-bit green channel difference from the previous pixel -32..31
4-bit   red channel difference minus green channel difference -8..7
4-bit  blue channel difference minus green channel difference -8..7

The green channel is used to indicate the general direction of change and is
encoded in 6 bits. The red and blue channels (dr and db) base their diffs off
of the green channel difference and are encoded in 4 bits. I.e.:
	dr_dg = (cur_px.r - prev_px.r) - (cur_px.g - prev_px.g)
	db_dg = (cur_px.b - prev_px.b) - (cur_px.g - prev_px.g)

The difference to the current channel values are using a wraparound operation,
so "10 - 13" will result in 253, while "250 + 7" will result in 1.

Values are stored as unsigned integers with a bias of 32 for the green channel
and a bias of 8 for the red and blue channel.

The alpha value remains unchanged from the previous pixel.


.- SOI_OP_LUMA777 ------------------------------------------------------------.
|         Byte[0]         |         Byte[1]         |         Byte[2]         |
|  7  6  5  4  3  2  1  0 |  7  6  5  4  3  2  1  0 |  7  6  5  4  3  2  1  0 |
|-------+-------------------------+-----+----------------+--------------------|
|  1  1  0 |     db - dg          |      dr - dg         |       dg           |
`-----------------------------------------------------------------------------`
3-bit tag b10
7-bit  blue channel difference minus green channel difference -64..63
7-bit   red channel difference minus green channel difference -64..63
7-bit green channel difference from the previous pixel -64..63

The green channel is used to indicate the general direction of change and is
encoded in 8 bits. The red and blue channels (dr and db) base their diffs off
of the green channel difference and are encoded in 4 bits. I.e.:
	dr_dg = (cur_px.r - prev_px.r) - (cur_px.g - prev_px.g)
	db_dg = (cur_px.b - prev_px.b) - (cur_px.g - prev_px.g)

The difference to the current channel values are using a wraparound operation,
so "10 - 13" will result in 253, while "250 + 7" will result in 1.

Values are stored as unsigned integers with a bias of 64 for each channel.

The alpha value remains unchanged from the previous pixel.


.- SOI_OP_RUN ------------.
|         Byte[0]         |
|  7  6  5  4  3  2  1  0 |
|----------+--------------|
|  1  1  0 |     run      |
`-------------------------`
2-bit tag b11
6-bit run-length repeating the previous pixel: 1..30

The run-length is stored with a bias of -1. Note that the run-lengths 31 and 32
(b111110 and b111111) are illegal as they are occupied by the SOI_OP_RGB and
SOI_OP_RGBA tags.


.- SOI_OP_RGB ------------------------------------------.
|         Byte[0]         | Byte[1] | Byte[2] | Byte[3] |
|  7  6  5  4  3  2  1  0 | 7 .. 0  | 7 .. 0  | 7 .. 0  |
|-------------------------+---------+---------+---------|
|  1  1  1  1  1  1  1  0 |   red   |  green  |  blue   |
`-------------------------------------------------------`
8-bit tag b11111110
8-bit   red channel value
8-bit green channel value
8-bit  blue channel value

The alpha value remains unchanged from the previous pixel.


.- SOI_OP_RGBA ---------------------------.
|         Byte[0]         | Byte[1] | ... |
|  7  6  5  4  3  2  1  0 | 7 .. 0  | ... |
|-------------------------+---------+-----|
|  1  1  1  1  1  1  1  1 |  alpha  | ... |
`-----------------------------------------`
8-bit tag b11111111
8-bit alpha channel value
RGB encoded as one of SOI_OP_LUMA222, SOI_OP_LUMA555, SOI_OP_LUMA777, SOI_OP_RGB

*/


/* -----------------------------------------------------------------------------
Header - Public functions */

#ifndef SOI_H
#define SOI_H

#ifndef SOIDEF
#define SOIDEF static
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* A pointer to a soi_desc struct has to be supplied to all of soi's functions.
It describes either the input format (for soi_write and soi_encode), or is
filled with the description read from the file header (for soi_read and
soi_decode).

The colorspace in this soi_desc is an enum where
	0 = sRGB, i.e. gamma scaled RGB channels and a linear alpha channel
	1 = all channels are linear
You may use the constants SOI_SRGB or SOI_LINEAR. The colorspace is purely
informative. It will be saved to the file header, but does not affect
how chunks are en-/decoded. */

#define SOI_SRGB   0
#define SOI_LINEAR 1

typedef struct {
	unsigned int width;
	unsigned int height;
	unsigned char channels;
	unsigned char colorspace;
} soi_desc;

#ifndef SOI_NO_STDIO

/* Encode raw RGB or RGBA pixels into a SOI image and write it to the file
system. The soi_desc struct must be filled with the image width, height,
number of channels (3 = RGB, 4 = RGBA) and the colorspace.

The function returns 0 on failure (invalid parameters, or fopen or malloc
failed) or the number of bytes written on success. */

SOIDEF
int soi_write(const char *filename, const void *data, const soi_desc *desc);


/* Read and decode a SOI image from the file system. If channels is 0, the
number of channels from the file header is used. If channels is 3 or 4 the
output format will be forced into this number of channels.

The function either returns NULL on failure (invalid data, or malloc or fopen
failed) or a pointer to the decoded pixels. On success, the soi_desc struct
will be filled with the description from the file header.

The returned pixel data should be free()d after use. */

SOIDEF
void *soi_read(const char *filename, soi_desc *desc, int channels);

#endif /* SOI_NO_STDIO */


/* Encode raw RGB or RGBA pixels into a SOI image in memory.

The function either returns NULL on failure (invalid parameters or malloc
failed) or a pointer to the encoded data on success. On success the out_len
is set to the size in bytes of the encoded data.

The returned soi data should be free()d after use. */

SOIDEF
void *soi_encode(const void *data, const soi_desc *desc, int *out_len);


/* Decode a SOI image from memory.

The function either returns NULL on failure (invalid parameters or malloc
failed) or a pointer to the decoded pixels. On success, the soi_desc struct
is filled with the description from the file header.

The returned pixel data should be free()d after use. */

SOIDEF
void *soi_decode(const void *data, int size, soi_desc *desc, int channels);


#ifdef __cplusplus
}
#endif
#endif /* SOI_H */


/* -----------------------------------------------------------------------------
Implementation */

#ifdef SOI_IMPLEMENTATION
#include <stdlib.h>
#include <string.h>

#ifndef SOI_MALLOC
	#define SOI_MALLOC(sz) malloc(sz)
	#define SOI_FREE(p)    free(p)
#endif
#ifndef SOI_ZEROARR
	#define SOI_ZEROARR(a) memset((a),0,sizeof(a))
#endif

#define SOI_OP_LUMA555 0x00 /* 0xxxxxxx */
#define SOI_OP_LUMA222 0x80 /* 10xxxxxx */
#define SOI_OP_LUMA777 0xc0 /* 110xxxxx */
#define SOI_OP_RUN     0xe0 /* 111xxxxx */
#define SOI_OP_RGB     0xfe /* 11111110 */
#define SOI_OP_RGBA    0xff /* 11111111 */

#define SOI_MASK_1     0x80 /* 10000000 */
#define SOI_MASK_2     0xc0 /* 11000000 */
#define SOI_MASK_3     0xe0 /* 11100000 */

#define SOI_MAGIC \
	(((unsigned int)'s') << 24 | ((unsigned int)'o') << 16 | \
	 ((unsigned int)'i') <<  8 | ((unsigned int)'f'))
#define SOI_HEADER_SIZE 14

/* 2GB is the max file size that this implementation can safely handle. We guard
against anything larger than that, assuming the worst case with 5 bytes per
pixel, rounded down to a nice clean value. 400 million pixels ought to be
enough for anybody. */
#define SOI_PIXELS_MAX ((unsigned int)400000000)

typedef union {
	struct { unsigned char r, g, b, a; } rgba;
	unsigned int v;
} soi_rgba_t;

static const unsigned char soi_padding[8] = {0,0,0,0,0,0,0,1};

static void soi_write_32(unsigned char *bytes, int *p, unsigned int v) {
	bytes[(*p)++] = (0xff000000 & v) >> 24;
	bytes[(*p)++] = (0x00ff0000 & v) >> 16;
	bytes[(*p)++] = (0x0000ff00 & v) >> 8;
	bytes[(*p)++] = (0x000000ff & v);
}

static unsigned int soi_read_32(const unsigned char *bytes, int *p) {
	unsigned int a = bytes[(*p)++];
	unsigned int b = bytes[(*p)++];
	unsigned int c = bytes[(*p)++];
	unsigned int d = bytes[(*p)++];
	return a << 24 | b << 16 | c << 8 | d;
}

const unsigned char optable[128]={0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};

#define SOI_RGB_ENC_SCALAR do{\
	signed char vr = px.rgba.r - px_prev.rgba.r;\
	signed char vg = px.rgba.g - px_prev.rgba.g;\
	signed char vb = px.rgba.b - px_prev.rgba.b;\
	signed char vg_r = vr - vg;\
	signed char vg_b = vb - vg;\
	unsigned char ar = (vg_r<0)?(-vg_r)-1:vg_r;\
	unsigned char ag = (vg<0)?(-vg)-1:vg;\
	unsigned char ab = (vg_b<0)?(-vg_b)-1:vg_b;\
	unsigned char argb = ar|ag|ab;\
	switch(optable[argb]){\
		case 0:\
			bytes[p++] = SOI_OP_LUMA222 | ((vg_r + 2) << 4) | ((vg_b + 2) << 2) | (vg + 2);\
			break;\
		case 1:\
			bytes[p++] = SOI_OP_LUMA555    | ((vg_b   + 16) << 2) | ((vg_r + 16)>>3);\
			bytes[p++] = (((vg_r + 16) & 7) << 5) | (vg +  16);\
			break;\
		case 2:\
			bytes[p++] = SOI_OP_LUMA777     | ((vg_b + 64)>>2);\
			bytes[p++] = (((vg_b+64)&3)<<6) | ((vg_r + 64)>>1);\
			bytes[p++] = (((vg_r+64)&1)<<7) | (vg+64);\
			break;\
		case 3:\
			bytes[p++] = SOI_OP_RGB;\
			bytes[p++] = px.rgba.r;\
			bytes[p++] = px.rgba.g;\
			bytes[p++] = px.rgba.b;\
			break;\
	}\
}while(0)

SOIDEF
void *soi_encode(const void *data, const soi_desc *desc, int *out_len) {
	int i, max_size, p, run;
	int px_len, px_end, px_pos, channels;
	unsigned char *bytes;
	const unsigned char *pixels;
	soi_rgba_t px, px_prev;

	if (
		data == NULL || out_len == NULL || desc == NULL ||
		desc->width == 0 || desc->height == 0 ||
		desc->channels < 3 || desc->channels > 4 ||
		desc->colorspace > 1 ||
		desc->height >= SOI_PIXELS_MAX / desc->width
	) {
		return NULL;
	}

	max_size =
		desc->width * desc->height * (desc->channels + 1) +
		SOI_HEADER_SIZE + sizeof(soi_padding);

	p = 0;
	bytes = (unsigned char *) SOI_MALLOC(max_size);
	if (!bytes) {
		return NULL;
	}

	soi_write_32(bytes, &p, SOI_MAGIC);
	soi_write_32(bytes, &p, desc->width);
	soi_write_32(bytes, &p, desc->height);
	bytes[p++] = desc->channels;
	bytes[p++] = desc->colorspace;

	pixels = (const unsigned char *)data;

	run = 0;
	px_prev.rgba.r = 0;
	px_prev.rgba.g = 0;
	px_prev.rgba.b = 0;
	px_prev.rgba.a = 255;
	px = px_prev;

	px_len = desc->width * desc->height * desc->channels;
	px_end = px_len - desc->channels;
	channels = desc->channels;

	if (channels == 4) {
		for (px_pos = 0; px_pos < px_len; px_pos += 4) {
			px.rgba.r = pixels[px_pos + 0];
			px.rgba.g = pixels[px_pos + 1];
			px.rgba.b = pixels[px_pos + 2];
			px.rgba.a = pixels[px_pos + 3];

			while(px.v == px_prev.v) {
				++run;
				if(px_pos == px_end) {
					bytes[p++] = SOI_OP_RUN | (run - 1);
					goto DONE;
				}
				else if (run == 30) {
					bytes[p++] = SOI_OP_RUN | (run - 1);
					run = 0;
				}
				px_pos+=4;
				px.rgba.r = pixels[px_pos + 0];
				px.rgba.g = pixels[px_pos + 1];
				px.rgba.b = pixels[px_pos + 2];
				px.rgba.a = pixels[px_pos + 3];
			}
			if (run) {
				bytes[p++] = SOI_OP_RUN | (run - 1);
				run = 0;
			}

			if(px.rgba.a!=px_prev.rgba.a){
				bytes[p++] = SOI_OP_RGBA;
				bytes[p++] = px.rgba.a;
			}

			SOI_RGB_ENC_SCALAR;
			px_prev = px;
		}
	}
	else {
		for (px_pos = 0; px_pos < px_len; px_pos += 3) {
			px.rgba.r = pixels[px_pos + 0];
			px.rgba.g = pixels[px_pos + 1];
			px.rgba.b = pixels[px_pos + 2];

			while(px.v == px_prev.v) {
				++run;
				if(px_pos == px_end) {
					bytes[p++] = SOI_OP_RUN | (run - 1);
					goto DONE;
				}
				else if (run == 30) {
					bytes[p++] = SOI_OP_RUN | (run - 1);
					run = 0;
				}
				px_pos+=3;
				px.rgba.r = pixels[px_pos + 0];
				px.rgba.g = pixels[px_pos + 1];
				px.rgba.b = pixels[px_pos + 2];
			}
			if (run) {
				bytes[p++] = SOI_OP_RUN | (run - 1);
				run = 0;
			}

			SOI_RGB_ENC_SCALAR;
			px_prev = px;
		}
	}
	DONE:

	for (i = 0; i < (int)sizeof(soi_padding); i++) {
		bytes[p++] = soi_padding[i];
	}

	*out_len = p;
	return bytes;
}

SOIDEF
void *soi_decode(const void *data, int size, soi_desc *desc, int channels) {
	const unsigned char *bytes;
	unsigned int header_magic;
	unsigned char *pixels;
	soi_rgba_t px;
	int px_len, chunks_len, px_pos;
	int p = 0, run = 0;

	if (
		data == NULL || desc == NULL ||
		(channels != 0 && channels != 3 && channels != 4) ||
		size < SOI_HEADER_SIZE + (int)sizeof(soi_padding)
	) {
		return NULL;
	}

	bytes = (const unsigned char *)data;

	header_magic = soi_read_32(bytes, &p);
	desc->width = soi_read_32(bytes, &p);
	desc->height = soi_read_32(bytes, &p);
	desc->channels = bytes[p++];
	desc->colorspace = bytes[p++];

	if (
		desc->width == 0 || desc->height == 0 ||
		desc->channels < 3 || desc->channels > 4 ||
		desc->colorspace > 1 ||
		header_magic != SOI_MAGIC ||
		desc->height >= SOI_PIXELS_MAX / desc->width
	) {
		return NULL;
	}

	if (channels == 0) {
		channels = desc->channels;
	}

	px_len = desc->width * desc->height * channels;
	pixels = (unsigned char *) SOI_MALLOC(px_len);
	if (!pixels) {
		return NULL;
	}

	px.rgba.r = 0;
	px.rgba.g = 0;
	px.rgba.b = 0;
	px.rgba.a = 255;

	chunks_len = size - (int)sizeof(soi_padding);
	for (px_pos = 0; px_pos < px_len; px_pos += channels) {
		if (run > 0) {
			run--;
		}
		else if (p < chunks_len) {
			OP_RGBA_GOTO:
			int b1 = bytes[p++];
			if (b1 == SOI_OP_RGB) {
				px.rgba.r = bytes[p++];
				px.rgba.g = bytes[p++];
				px.rgba.b = bytes[p++];
			}
			else if (b1 == SOI_OP_RGBA) {
				px.rgba.a = bytes[p++];
				goto OP_RGBA_GOTO;
			}
			else if ((b1 & SOI_MASK_2) == SOI_OP_LUMA222) {
				int vg = (b1 & 3) - 2;
				px.rgba.r += vg - 2 + ((b1 >> 4) & 3);
				px.rgba.g += vg;
				px.rgba.b += vg - 2 + ((b1 >> 2) & 3);
			}
			else if ((b1 & SOI_MASK_1) == SOI_OP_LUMA555) {
				int b2 = bytes[p++];
				int vg = (b2 & 31) - 16;
				px.rgba.r += vg - 16 + (((b1&3)<<3) | (b2>>5));
				px.rgba.g += vg;
				px.rgba.b += vg - 16 +  ((b1 >>2)&31);
			}
			else if ((b1 & SOI_MASK_3) == SOI_OP_LUMA777) {
				int b2 = bytes[p++];
				int b3 = bytes[p++];
				int vg = (b3 & 0x7f) - 64;
				px.rgba.r += vg - 64 + ((b2 & 0x3f)<<1) + (b3>>7);
				px.rgba.g += vg;
				px.rgba.b += vg - 64 + ((b1 & 0x1f)<<2) + (b2>>6);
			}
			else if ((b1 & SOI_MASK_3) == SOI_OP_RUN) {
				run = (b1 & 0x1f);
			}
		}
		pixels[px_pos + 0] = px.rgba.r;
		pixels[px_pos + 1] = px.rgba.g;
		pixels[px_pos + 2] = px.rgba.b;
		
		if (channels == 4) {
			pixels[px_pos + 3] = px.rgba.a;
		}
	}

	return pixels;
}

#ifndef SOI_NO_STDIO
#include <stdio.h>

SOIDEF
int soi_write(const char *filename, const void *data, const soi_desc *desc) {
	FILE *f = fopen(filename, "wb");
	int size, err;
	void *encoded;

	if (!f) {
		return 0;
	}

	encoded = soi_encode(data, desc, &size);
	if (!encoded) {
		fclose(f);
		return 0;
	}

	fwrite(encoded, 1, size, f);
	fflush(f);
	err = ferror(f);
	fclose(f);

	SOI_FREE(encoded);
	return err ? 0 : size;
}

SOIDEF
void *soi_read(const char *filename, soi_desc *desc, int channels) {
	FILE *f = fopen(filename, "rb");
	int size, bytes_read;
	void *pixels, *data;

	if (!f) {
		return NULL;
	}

	fseek(f, 0, SEEK_END);
	size = ftell(f);
	if (size <= 0 || fseek(f, 0, SEEK_SET) != 0) {
		fclose(f);
		return NULL;
	}

	data = SOI_MALLOC(size);
	if (!data) {
		fclose(f);
		return NULL;
	}

	bytes_read = fread(data, 1, size, f);
	fclose(f);
	pixels = (bytes_read != size) ? NULL : soi_decode(data, bytes_read, desc, channels);
	SOI_FREE(data);
	return pixels;
}

#endif /* SOI_NO_STDIO */
#endif /* SOI_IMPLEMENTATION */
