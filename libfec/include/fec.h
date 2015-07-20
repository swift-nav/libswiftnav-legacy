/* User include file for libfec
 * Copyright 2004, Phil Karn, KA9Q
 * May be used under the terms of the GNU Lesser General Public License (LGPL)
 */

#ifndef _FEC_H_
#define _FEC_H_

/* r=1/2 k=7 convolutional encoder polynomials
 * The NASA-DSN convention is to use V27POLYA inverted, then V27POLYB
 * The CCSDS/NASA-GSFC convention is to use V27POLYB, then V27POLYA inverted
 */
#define V27POLYA  0x4f
#define V27POLYB  0x6d
#define SBAS_BLK_LEN 250

typedef union {
  unsigned int w[64];
} metric_t;

typedef union {
  unsigned long w[2];
} decision_t;

union branchtab27 {
  unsigned char c[32];
};

/* State info for instance of Viterbi decoder
 */
struct v27 {
  metric_t metrics1; /* path metric buffer 1 */
  metric_t metrics2; /* path metric buffer 2 */
  decision_t *dp;          /* Pointer to current decision */
  metric_t *old_metrics,*new_metrics; /* Pointers to path metrics, swapped on every bit */
  decision_t decisions[SBAS_BLK_LEN * 6 + 6];   /* Beginning of decisions for block */
};

void *create_viterbi27(int len);
void set_viterbi27_polynomial(int polys[2]);
void init_viterbi27(struct v27 *vp, int starting_state);
int update_viterbi27_blk(struct v27 *p, unsigned char sym[], int npairs);
int chainback_viterbi27(struct v27 *p, unsigned char *data, unsigned int nbits,
                        unsigned int endstate);

static inline int parity(int x)
{
  x ^= x >> 16;
  x ^= x >> 8;
  x ^= x >> 4;
  x &= 0xf;
  return (0x6996 >> x) & 1;
}

#endif
