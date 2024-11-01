/*

Copyright (c) 2021, Dominic Szablewski - https://phoboslab.org
SPDX-License-Identifier: MIT


Command line tool to convert between png <> qoi format

Requires:
	-"stb_image.h" (https://github.com/nothings/stb/blob/master/stb_image.h)
	-"stb_image_write.h" (https://github.com/nothings/stb/blob/master/stb_image_write.h)
	-"qoi.h" (https://github.com/phoboslab/qoi/blob/master/qoi.h)

Compile with: 
	gcc qoiconv.c -std=c99 -O3 -o qoiconv

*/
#include <fstream>
#include <string>
using namespace std;

#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_PNG
#define STBI_NO_LINEAR
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define QOI_IMPLEMENTATION
#define QOI_NO_STDIO
#include "qoi.h"

#define SOI_IMPLEMENTATION
#define SOI_NO_STDIO
#include "soi.h"

#define ROI_IMPLEMENTATION
#define ROI_NO_STDIO
#include "roi.h"

#define MINIZ_IMPLEMENTATION
#define MINIZ_NO_STDIO
#include "miniz.h"

#define STR_ENDS_WITH(S, E) (strcmp(S + strlen(S) - (sizeof(E)-1), E) == 0)

int compressData(const char* input, int inputSize, char* output, int& outputSize) {
	int err;
	z_stream zs;
	zs.zalloc = Z_NULL;
	zs.zfree = Z_NULL;
	zs.opaque = Z_NULL;
	zs.next_in = (Bytef *)input;
	zs.avail_in = (uInt)inputSize;
	zs.avail_out = (uInt)outputSize;
	zs.next_out = (Bytef *)output;
	err = deflateInit2(&zs, Z_BEST_COMPRESSION, Z_DEFLATED, 15, 8, Z_DEFAULT_STRATEGY);
	if (err != Z_OK) return err;
	err = deflate(&zs, Z_FINISH);
	deflateEnd(&zs);
	outputSize = zs.total_out;
	return err;
}

int uncompressData(const char* input, int inputSize, char* output, int& outputSize) {
	int err;
	z_stream zs;
	zs.zalloc = Z_NULL;
	zs.zfree = Z_NULL;
	zs.opaque = Z_NULL;
	zs.avail_in = (uInt)inputSize;
	zs.next_in = (Bytef *)input;
	zs.avail_out = (uInt)outputSize;
	zs.next_out = (Bytef *)output;
	err = inflateInit2(&zs, 15);
	if (err != Z_OK) return err;
	err = inflate(&zs, Z_FINISH);
	if (err != Z_STREAM_END) {
		inflateEnd(&zs);
		return err == Z_OK ? Z_BUF_ERROR : err;
	}
	outputSize = zs.total_out;
	err = inflateEnd(&zs);
	return err;
}

static string loadStringFromFile(const string filename) {
	string res;
	if (filename == "") return res;
	ifstream infile;
	infile.open(filename.c_str(), ios::in | ios::binary | ios::ate);
	if (infile.is_open()) {
		int size = (int)infile.tellg();
		char* bufferP = new char[size];
		infile.seekg (0, ios::beg);
		infile.read(bufferP, size);
		if (!infile) size = (int)infile.gcount();
		infile.close();
		res = string(reinterpret_cast<const char*>(bufferP), size);
		delete[] bufferP;
	}
	return res;
}

static int saveStringToFile(string s, const string filename) {
	ofstream myfile(filename.c_str(), ios::out | ios::binary);
	if (myfile.is_open()) {
		myfile.write(s.c_str(), s.size());
		myfile.close();
		return 0;
	}
	return -1;
}

int main(int argc, char **argv) {
	if (argc < 3) {
		puts("Usage: qoiconv <infile> <outfile>");
		puts("Examples:");
		puts("  qoiconv input.png output.roi");
		puts("  qoiconv input.qoi output.png");
		puts("Supported formats: qoi, qoz, roi, roz, soi, soz, png");
		exit(1);
	}

	void *pixels = NULL;
	int w, h, channels;
	if (STR_ENDS_WITH(argv[1], ".png")) {
		if(!stbi_info(argv[1], &w, &h, &channels)) {
			printf("Couldn't read header %s\n", argv[1]);
			exit(1);
		}

		// Force all odd encodings to be RGBA
		if(channels != 3) {
			channels = 4;
		}

		pixels = (void *)stbi_load(argv[1], &w, &h, NULL, channels);
	} else if (STR_ENDS_WITH(argv[1], ".qoi")) {
		qoi_desc desc;
		string s = loadStringFromFile(argv[1]);
		pixels = qoi_decode(s.c_str(), s.size(), &desc, 0);
		// pixels = qoi_read(argv[1], &desc, 0);
		channels = desc.channels;
		w = desc.width;
		h = desc.height;
	} else if (STR_ENDS_WITH(argv[1], ".qoz")) {
		qoi_desc desc;
		string s = loadStringFromFile(argv[1]);
		int outputSize = 0;
		char* c = (char *)(&outputSize);
		*c++ = s[s.size() - 4];
		*c++ = s[s.size() - 3];
		*c++ = s[s.size() - 2];
		*c = s[s.size() - 1];
		char* output = new char[outputSize];
		uncompressData(s.c_str(), s.size(), output, outputSize);
		pixels = qoi_decode(output, outputSize, &desc, 0);
		// pixels = qoi_read(argv[1], &desc, 0);
		channels = desc.channels;
		w = desc.width;
		h = desc.height;
		delete[] output;
	} else if (STR_ENDS_WITH(argv[1], ".roi")) {
		roi_desc desc;
		string s = loadStringFromFile(argv[1]);
		pixels = roi_decode(s.c_str(), s.size(), &desc, 0);
		channels = desc.channels;
		w = desc.width;
		h = desc.height;
	} else if (STR_ENDS_WITH(argv[1], ".roz")) {
		roi_desc desc;
		string s = loadStringFromFile(argv[1]);
		int outputSize = 0;
		char* c = (char *)(&outputSize);
		*c++ = s[s.size() - 4];
		*c++ = s[s.size() - 3];
		*c++ = s[s.size() - 2];
		*c = s[s.size() - 1];
		char* output = new char[outputSize];
		uncompressData(s.c_str(), s.size(), output, outputSize);
		pixels = roi_decode(output, outputSize, &desc, 0);
		channels = desc.channels;
		w = desc.width;
		h = desc.height;
		delete[] output;
	} else if (STR_ENDS_WITH(argv[1], ".soi")) {
		soi_desc desc;
		string s = loadStringFromFile(argv[1]);
		pixels = soi_decode(s.c_str(), s.size(), &desc, 0);
		channels = desc.channels;
		w = desc.width;
		h = desc.height;
	} else if (STR_ENDS_WITH(argv[1], ".soz")) {
		soi_desc desc;
		string s = loadStringFromFile(argv[1]);
		int outputSize = 0;
		char* c = (char *)(&outputSize);
		*c++ = s[s.size() - 4];
		*c++ = s[s.size() - 3];
		*c++ = s[s.size() - 2];
		*c = s[s.size() - 1];
		char* output = new char[outputSize];
		uncompressData(s.c_str(), s.size(), output, outputSize);
		pixels = soi_decode(output, outputSize, &desc, 0);
		channels = desc.channels;
		w = desc.width;
		h = desc.height;
		delete[] output;
	}

	if (pixels == NULL) {
		printf("Couldn't load/decode %s\n", argv[1]);
		exit(1);
	}

	int encoded = 0;
	if (STR_ENDS_WITH(argv[2], ".png")) {
		encoded = stbi_write_png(argv[2], w, h, channels, pixels, 0);
	} else if (STR_ENDS_WITH(argv[2], ".qoi")) {
		qoi_desc d;
		d.width = w;
		d.height = h;
		d.channels = channels;
		d.colorspace = QOI_SRGB;
		int l = 0;
		char* data = (char *)qoi_encode(pixels, &d, &l);
		if (l > 0) {
			encoded = 1;
			string s = string(data, l);
			saveStringToFile(s, argv[2]);
		}
	} else if (STR_ENDS_WITH(argv[2], ".qoz")) {
		qoi_desc d;
		d.width = w;
		d.height = h;
		d.channels = channels;
		d.colorspace = QOI_SRGB;
		int l = 0;
		char* data = (char *)qoi_encode(pixels, &d, &l);
		if (l > 0) {
			encoded = 1;
			string s = string(data, l);
			int inputSize = s.size();
			int outputSize = compressBound(inputSize);
			char* output = new char[outputSize];
			compressData(s.c_str(), inputSize, output, outputSize);
			string s2 = string(output, outputSize) + "xxxx";
			char* c = (char *)(&inputSize);
			s2[s2.size() - 4] = *c++;
			s2[s2.size() - 3] = *c++;
			s2[s2.size() - 2] = *c++;
			s2[s2.size() - 1] = *c;
			saveStringToFile(s2, argv[2]);
			delete[] output;
		}
	} else if (STR_ENDS_WITH(argv[2], ".roi")) {
		roi_desc d;
		d.width = w;
		d.height = h;
		d.channels = channels;
		d.colorspace = ROI_SRGB;
		int l = 0;
		roi_options opt;
		opt.path = scalar;
		opt.mlut = 0;
		char* data = (char *)roi_encode(pixels, &d, &l, &opt);
		if (l > 0) {
			encoded = 1;
			string s = string(data, l);
			saveStringToFile(s, argv[2]);
		}
	} else if (STR_ENDS_WITH(argv[2], ".roz")) {
		roi_desc d;
		d.width = w;
		d.height = h;
		d.channels = channels;
		d.colorspace = ROI_SRGB;
		int l = 0;
		roi_options opt;
		opt.path = scalar;
		opt.mlut = 0;
		char* data = (char *)roi_encode(pixels, &d, &l, &opt);
		if (l > 0) {
			encoded = 1;
			string s = string(data, l);
			int inputSize = s.size();
			int outputSize = compressBound(inputSize);
			char* output = new char[outputSize];
			compressData(s.c_str(), inputSize, output, outputSize);
			string s2 = string(output, outputSize) + "xxxx";
			char* c = (char *)(&inputSize);
			s2[s2.size() - 4] = *c++;
			s2[s2.size() - 3] = *c++;
			s2[s2.size() - 2] = *c++;
			s2[s2.size() - 1] = *c;
			saveStringToFile(s2, argv[2]);
			delete[] output;
		}
	} else if (STR_ENDS_WITH(argv[2], ".soi")) {
		soi_desc d;
		d.width = w;
		d.height = h;
		d.channels = channels;
		d.colorspace = SOI_SRGB;
		// encoded = soi_write((string(argv[2]) + "x").c_str(), pixels, &d);
		int l = 0;
		char* data = (char *)soi_encode(pixels, &d, &l);
		if (l > 0) {
			encoded = 1;
			string s = string(data, l);
			saveStringToFile(s, argv[2]);
		}
	} else if (STR_ENDS_WITH(argv[2], ".soz")) {
		soi_desc d;
		d.width = w;
		d.height = h;
		d.channels = channels;
		d.colorspace = SOI_SRGB;
		int l = 0;
		char* data = (char *)soi_encode(pixels, &d, &l);
		if (l > 0) {
			encoded = 1;
			string s = string(data, l);
			int inputSize = s.size();
			int outputSize = compressBound(inputSize);
			char* output = new char[outputSize];
			compressData(s.c_str(), inputSize, output, outputSize);
			string s2 = string(output, outputSize) + "xxxx";
			char* c = (char *)(&inputSize);
			s2[s2.size() - 4] = *c++;
			s2[s2.size() - 3] = *c++;
			s2[s2.size() - 2] = *c++;
			s2[s2.size() - 1] = *c;
			saveStringToFile(s2, argv[2]);
			delete[] output;
		}
	}

	if (!encoded) {
		printf("Couldn't write/encode %s\n", argv[2]);
		exit(1);
	}

	free(pixels);
	return 0;
}
