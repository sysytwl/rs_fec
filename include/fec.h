/**
 * zfec -- fast forward error correction library with Python interface
 *
 * Copyright (C) 2007-2008 Allmydata, Inc.
 * Author: Zooko Wilcox-O'Hearn
 *
 * This file is part of zfec.
 *
 * See README.rst for licensing information.
 */

/*
 * Much of this work is derived from the "fec" software by Luigi Rizzo, et
 * al., the copyright notice and licence terms of which are included below
 * for reference.
 *
 * fec.h -- forward error correction based on Vandermonde matrices
 * 980614
 * (C) 1997-98 Luigi Rizzo (luigi@iet.unipi.it)
 *
 * Portions derived from code by Phil Karn (karn@ka9q.ampr.org),
 * Robert Morelos-Zaragoza (robert@spectra.eng.hawaii.edu) and Hari
 * Thirumoorthy (harit@spectra.eng.hawaii.edu), Aug 1995
 *
 * Modifications by Dan Rubenstein (see Modifications.txt for
 * their description.
 * Modifications (C) 1998 Dan Rubenstein (drubenst@cs.umass.edu)
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:

 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above
 *    copyright notice, this list of conditions and the following
 *    disclaimer in the documentation and/or other materials
 *    provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
 * TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 */

/* Modified by JackShenYt */
#pragma once
#include <stddef.h>
#include <cstdio>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

static constexpr unsigned BLOCK_NUMS[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};



class ZFE_FEC{
public:
  typedef unsigned char gf;

  struct fec_t {
    unsigned long magic;
    unsigned short k, n;                     /* parameters of the code */
    gf* enc_matrix;
  };
  //fec_t* fec_type;

  void init_fec (void);

  /**
   * @param src the "primary blocks" i.e. the chunks of the input data
   * @param fecs buffers into which the secondary blocks will be written
   * @param block_nums the numbers of the desired check blocks (the id >= k) which fec_encode() will produce and store into the buffers of the fecs parameter
   * @param num_block_nums the length of the block_nums array
   * @param sz size of a packet in bytes
   */
  void fec_encode(const fec_t* code, const gf* const* const src, gf* const* const fecs, const unsigned* const block_nums, size_t num_block_nums, size_t sz);

  void fec_encode_block(const fec_t* code, const uint8_t *const *const src, uint8_t* const fec, const unsigned* const block_nums, int fec_block_index, size_t sz);

  /**
   * @param inpkts an array of packets (size k); If a primary block, i, is present then it must be at index i. Secondary blocks can appear anywhere.
   * @param outpkts an array of buffers into which the reconstructed output packets will be written (only packets which are not present in the inpkts input will be reconstructed and written to outpkts)
   * @param index an array of the blocknums of the packets in inpkts
   * @param sz size of a packet in bytes
   */
  void fec_decode(const fec_t* code, const uint8_t*const*const inpkts, uint8_t*const*const outpkts, const unsigned*const index, size_t sz);

  void fec_free (fec_t** p);

  /**
   * @param k the number of blocks required to reconstruct
   * @param m the total number of blocks created
   */
  ZFE_FEC::fec_t* fec_new(unsigned short k, unsigned short n);


  int fec_initialized = 0;

  uint8_t gf_exp[510];  /* index->poly form conversion table    */
  gf gf_log[256];  /* poly->index form conversion table    */
  gf inverse[256]; /* inverse of field elem.               */ /* inv[\alpha**i]=\alpha**(GF_SIZE-i-1) */

  /*
  * gf_mul(x,y) multiplies two numbers.  It is much faster to use a
  * multiplication table.
  *
  * USE_GF_MULC, GF_MULC0(c) and GF_ADDMULC(x) can be used when multiplying
  * many numbers by the same constant. In this case the first call sets the
  * constant, and others perform the multiplications.  A value related to the
  * multiplication is held in a local variable declared with USE_GF_MULC . See
  * usage in _addmul1().
  */
  typedef uint8_t gf_mul_table_t[256][256];
  typedef gf_mul_table_t* gf_mul_table_p;

  gf_mul_table_p gf_mul_table;


  private:
  /*
  * Primitive polynomials - see Lin & Costello, Appendix A,
  static const char * const Pp;
  */
  const char * const Pp = "101110001";

  /*
  * modnn(x) computes x % GF_SIZE, where GF_SIZE is 2**GF_BITS - 1,
  * without a slow divide.
  */
  gf modnn(int x);

  /*
  * To speed up computations, we have tables for logarithm, exponent and
  * inverse of a number.  We use a table for multiplication as well (it takes
  * 64K, no big deal even on a PDA, especially because it can be
  * pre-initialized an put into a ROM!), otherwhise we use a table of
  * logarithms. In any case the macro gf_mul(x,y) takes care of
  * multiplications.
  */
  /*
  * initialize the data structures used for computations in GF.
  */
  void generate_gf (void);

  /*
  * Generate GF(2**m) from the irreducible polynomial p(X) in p[0]..p[m]
  * Lookup tables:
  *     index->polynomial form		gf_exp[] contains j= \alpha^i;
  *     polynomial form -> index form	gf_log[ j = \alpha^i ] = i
  * \alpha=x is the primitive element of GF(2^m)
  *
  * For efficiency, gf_exp[] has size 2*GF_SIZE, so that a simple
  * multiplication of two numbers can be resolved without calling modnn
  */
  void _init_mul_table(void);

  void _addmul1(gf* dst, const gf* src, gf c, size_t sz);



/*
 * computes C = AB where A is n*k, B is k*m, C is n*m
 */
void _matmul(gf * a, gf * b, gf * c, unsigned n, unsigned k, unsigned m);

/*
 * _invert_mat() takes a matrix and produces its inverse
 * k is the size of the matrix.
 * (Gauss-Jordan, adapted from Numerical Recipes in C)
 * Return non-zero if singular.
 */
void _invert_mat(gf* src, unsigned k);

/*
 * fast code for inverting a vandermonde matrix.
 *
 * NOTE: It assumes that the matrix is not singular and _IS_ a vandermonde
 * matrix. Only uses the second column of the matrix, containing the p_i's.
 *
 * Algorithm borrowed from "Numerical recipes in C" -- sec.2.8, but largely
 * revised for my purposes.
 * p = coefficients of the matrix (p_i)
 * q = values of the polynomial (known)
 */
void _invert_vdm (gf* src, unsigned k);

/**
 * Build decode matrix into some memory space.
 *
 * @param matrix a space allocated for a k by k matrix
 */
void build_decode_matrix_into_space(const fec_t* const code, const unsigned*const index, const unsigned k, gf* const matrix);

};

#if defined(_MSC_VER)
#define alloca _alloca
#else
#ifdef __GNUC__
#ifndef alloca
#define alloca(x) __builtin_alloca(x)
#endif
#else
#include <alloca.h>
#endif
#endif
