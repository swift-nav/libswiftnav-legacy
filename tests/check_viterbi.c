
#include <check.h>
#include <fec.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>

int compare_files(FILE *f1, FILE *f2)
{
  char ch1, ch2;
  int flag = 0;

  ch1 = fgetc(f1);
  ch2 = fgetc(f2);

  while (ch1 != EOF && ch2 != EOF) {
    if (ch1 != ch2) {
      flag = 1;
      return flag;
    }
  ch1 = fgetc(f1);
  ch2 = fgetc(f2);
  }
  if (ch1 == EOF && ch2 != ch1) {
    flag = 1;
    return flag;
  }

  if (ch2 == EOF && ch1 != ch2) {
    flag = 1;
    return flag;
  }

  return flag;
}

START_TEST(test_viterbi27)
{
  int i;
  struct v27 vp;
  FILE *waas_data = NULL;
  FILE *tmp = NULL;
  FILE *check = NULL;

  int output = open ("tmp.bin", O_WRONLY | O_CREAT, 0644);
  waas_data = fopen("waas_data.bin", "r");
  check = fopen("waas_check.bin", "r");

  fseek(waas_data, 0, SEEK_END);
  long fsize = ftell(waas_data);
  fseek(waas_data, 0, SEEK_SET);

  int nsymbols = fsize;
  int nbits = nsymbols/2 - 12;

  unsigned char data[nbits];
  memset(data, 0x42, nbits);

  unsigned char *buffer = (unsigned char *)malloc((fsize + 1) * sizeof(unsigned char));

  fread(buffer, fsize, 1, waas_data);
  fclose(waas_data);

  for(i = 0; i < fsize; i++) {
    char c = buffer[i];
    int digit = c - '0';
    if (digit == 1)
      buffer[i] = 0xff;
    else
      buffer[i] = 0x00;
  }

  init_viterbi27(&vp, 0);
  update_viterbi27_blk(&vp, buffer, 250 * 6);
  chainback_viterbi27(&vp, data, 250 * 6, 0);

  write(output, &data, 250 * 6);
  close(output);
  tmp = fopen("tmp.bin", "r");

  fail_unless(compare_files(check, tmp) == 0, "Viterbi decoder produced incorect file!");

  free(buffer);
  fclose(check);
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

