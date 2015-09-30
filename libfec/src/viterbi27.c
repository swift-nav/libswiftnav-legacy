/* K=7 r=1/2 Viterbi decoder in portable C
 * Copyright Feb 2004, Phil Karn, KA9Q
 * May be used under the terms of the GNU Lesser General Public License (LGPL)
 */

#include <stdlib.h>
#include "fec.h"

static union branchtab27 Branchtab27[2] __attribute__ ((aligned(16)));

void set_viterbi27_polynomial(int polys[2])
{
  int state;

  for(state = 0; state < 32; state++) {
    Branchtab27[0].c[state] = (polys[0] < 0) ^ par((2*state) & abs(polys[0])) ? 255 : 0;
    Branchtab27[1].c[state] = (polys[1] < 0) ^ par((2*state) & abs(polys[1])) ? 255 : 0;
  }
}

/* Create a new instance of a Viterbi decoder */
void init_viterbi27(struct v27 *vp, int starting_state)
{
  int i;
  int polys[2] = { V27POLYA, V27POLYB };

  set_viterbi27_polynomial(polys);

  for(i = 0; i < 64; i++)
    vp->metrics1.w[i] = 63;

  vp->old_metrics = &vp->metrics1;
  vp->new_metrics = &vp->metrics2;
  vp->dp = vp->decisions;
  vp->old_metrics->w[starting_state & 63] = 0; /* Bias known start state */
}

/* Viterbi chainback */
int chainback_viterbi27(
      struct v27 *vp,
      unsigned char *data, /* Decoded output data */
      unsigned int nbits, /* Number of data bits */
      unsigned int endstate) /* Terminal encoder state */
{
  decision_t *d;

  if(vp == NULL)
    return -1;

  d = vp->decisions;
  /* Make room beyond the end of the encoder register so we can
   * accumulate a full byte of decoded data
   */
  endstate %= 64;
  endstate <<= 2;

  /* The store into data[] only needs to be done every 8 bits.
   * But this avoids a conditional branch, and the writes will
   * combine in the cache anyway
   */
  d += 6; /* Look past tail */
  while(nbits-- != 0) {
    int k;

    k = (d[nbits].w[(endstate>>2)/32] >> ((endstate>>2)%32)) & 1;
    data[nbits>>3] = endstate = (endstate >> 1) | (k << 7);
  }

  return 0;
}

/* C-language butterfly */
#define BFLY(i) {\
unsigned int metric,m0,m1,decision;\
    metric = (Branchtab27[0].c[i] ^ sym0) + (Branchtab27[1].c[i] ^ sym1);\
    m0 = vp->old_metrics->w[i] + metric;\
    m1 = vp->old_metrics->w[i+32] + (510 - metric);\
    decision = (signed int)(m0-m1) > 0;\
    vp->new_metrics->w[2*i] = decision ? m1 : m0;\
    d->w[i/16] |= decision << ((2*i)&31);\
    m0 -= (metric+metric-510);\
    m1 += (metric+metric-510);\
    decision = (signed int)(m0-m1) > 0;\
    vp->new_metrics->w[2*i+1] = decision ? m1 : m0;\
    d->w[i/16] |= decision << ((2*i+1)&31);\
}

/* Update decoder with a block of demodulated symbols
 * Note that nbits is the number of decoded data bits, not the number
 * of symbols!
 */
int update_viterbi27_blk(struct v27 *vp, const unsigned char *syms, int nbits)
{
  void *tmp;
  decision_t *d;

  if(vp == NULL)
    return -1;

  d = (decision_t *)vp->dp;

  while(nbits--) {
    unsigned char sym0,sym1;

    d->w[0] = d->w[1] = 0;
    sym0 = *syms++;
    sym1 = *syms++;

    BFLY(0);
    BFLY(1);
    BFLY(2);
    BFLY(3);
    BFLY(4);
    BFLY(5);
    BFLY(6);
    BFLY(7);
    BFLY(8);
    BFLY(9);
    BFLY(10);
    BFLY(11);
    BFLY(12);
    BFLY(13);
    BFLY(14);
    BFLY(15);
    BFLY(16);
    BFLY(17);
    BFLY(18);
    BFLY(19);
    BFLY(20);
    BFLY(21);
    BFLY(22);
    BFLY(23);
    BFLY(24);
    BFLY(25);
    BFLY(26);
    BFLY(27);
    BFLY(28);
    BFLY(29);
    BFLY(30);
    BFLY(31);
    d++;

    /* Swap pointers to old and new metrics */
    tmp = vp->old_metrics;
    vp->old_metrics = vp->new_metrics;
    vp->new_metrics = tmp;
  }

  vp->dp = d;
  return 0;
}


void set_decisions_viterbi27(struct v27 *vp, decision_t *dec)
{
  vp->decisions = dec;
}
