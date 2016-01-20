#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>

#include <libfec/fec.h>

int compare_files(FILE *f1, FILE *f2)
{
  int c1, c2;
  int count = 0;

  while (1) {
    c1 = fgetc(f1);
    c2 = fgetc(f2);

    if (c1 != c2)
      return count;

    if (c1 == EOF)
      return count;

    count++;
  }

  return count;
}

START_TEST(test_viterbi27)
{
  #define HISTORY_LENGTH_BITS 64
  #define DECODE_LENGTH_BITS 32

  FILE *waas_symbols = fopen("v27_sym_waas.bin", "r");
  fail_if(NULL == waas_symbols, "Could not load WAAS data file");
  FILE *waas_bits = fopen("v27_bits_waas.bin", "r");
  fail_if(NULL == waas_bits, "Could not load WAAS data file");
  const char *tmp_file = "tmp.bin";

  /* Read in viterbi symbols from file */
  fseek(waas_symbols, 0, SEEK_END);
  long fsize = ftell(waas_symbols);
  fseek(waas_symbols, 0, SEEK_SET);

  unsigned char *symbol_buffer = (unsigned char *)malloc((fsize + 1) *
                                                         sizeof(unsigned char));
  fread(symbol_buffer, fsize, 1, waas_symbols);
  fclose(waas_symbols);

  /* One symbol per byte in input file */
  int nbits = fsize/2;

  /* Initialize viterbi decoder */
  v27_t v;
  v27_decision_t decisions[HISTORY_LENGTH_BITS];
  v27_poly_t v27_poly;
  signed char poly_bytes[] = {V27POLYA, V27POLYB};
  int v27_bit_length = -32; /* hack to process 32 bits before outputting data */

  v27_poly_init(&v27_poly, poly_bytes);
  v27_init(&v, decisions, HISTORY_LENGTH_BITS, &v27_poly, 0);

  /* Perform decoding, writing output to file */
  int output = open(tmp_file, O_WRONLY | O_CREAT, 0644);
  int bit_offset = 0;
  int output_length_bits = 0;
  while (bit_offset < nbits) {
    int bits_to_process = HISTORY_LENGTH_BITS - v27_bit_length;
    if (bit_offset + bits_to_process > nbits)
      break;

    v27_update(&v, &symbol_buffer[2*bit_offset], bits_to_process);
    v27_bit_length += bits_to_process;
    bit_offset += bits_to_process;

    unsigned char data[HISTORY_LENGTH_BITS / 8];
    memset(data, 0x43, HISTORY_LENGTH_BITS / 8);
    v27_chainback_likely(&v, data, HISTORY_LENGTH_BITS);

    /* write out the first DECODE_LENGTH bytes */
    write(output, data, DECODE_LENGTH_BITS / 8);
    v27_bit_length -= DECODE_LENGTH_BITS;
    output_length_bits += DECODE_LENGTH_BITS;
  }
  close(output);

  /* compare output with expected result */
  FILE *tmp = fopen(tmp_file, "r");
  int match_length_bytes = compare_files(waas_bits, tmp);

  fail_unless(match_length_bytes == output_length_bits / 8,
      "Viterbi decoder produced incorrect file!\n"
      "%d %d", match_length_bytes, output_length_bits / 8);

  free(symbol_buffer);
  fclose(waas_bits);
  fclose(tmp);
}
END_TEST

Suite* viterbi_suite(void)
{
  Suite *s = suite_create("Viterbi decoder 2/7");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_viterbi27);
  suite_add_tcase(s, tc_core);

  return s;
}
