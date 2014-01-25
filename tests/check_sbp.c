
#include <stdio.h>
#include <check.h>

#include <sbp.h>

extern sbp_msg_callbacks_node_t *sbp_msg_callbacks_head;

u32 dummy_wr = 0;
u32 dummy_rd = 0;
u8 dummy_buff[1024];

u32 dummy_write(u8 *buff, u32 n)
{
 u32 real_n = n;//(dummy_n > n) ? n : dummy_n;
 memcpy(dummy_buff + dummy_wr, buff, real_n);
 dummy_wr += real_n;
 return real_n;
}
u32 dummy_read(u8 *buff, u32 n)
{
 u32 real_n = n;//(dummy_n > n) ? n : dummy_n;
 memcpy(buff, dummy_buff + dummy_rd, real_n);
 dummy_rd += real_n;
 return real_n;
}

typedef struct {
 u8 x, y;
} awesome_struct;

void printy_callback(u16 sender_id, u8 len, u8 msg[])
{
  printf("MSG id: 0x%04X, len: %d\n", sender_id, len);
  for (u8 i=0; i<len; i++)
    printf("0x%02X ", msg[i]);
  printf("\n");
}

void test_callback(u16 sender_id __attribute__((unused)),
                   u8 len __attribute__((unused)),
                   u8 msg[] __attribute__((unused)))
{
  /* Do nothing. */
}
void test_callback2(u16 sender_id __attribute__((unused)),
                    u8 len __attribute__((unused)),
                    u8 msg[] __attribute__((unused)))
{
  /* Do nothing. */
}

START_TEST(test_sbp_send_message)
{
  awesome_struct s = { 0x22, 0x33 };

  sbp_send_message(0x4455, 0x6677, sizeof(s), (u8*)&s, &dummy_write);

  static sbp_msg_callbacks_node_t test_node;
  sbp_register_callback(0x4456, &printy_callback, &test_node);

  sbp_state_t state;
  sbp_state_init(&state);
  s8 ret;
  while (dummy_wr > dummy_rd) {
    if ((ret = sbp_process(&state, &dummy_read)) != SBP_OK)
      printf("Uh oh! %d\n", ret);
  }

  /*
  ret = getbitu(test_data, 10, 13);
  fail_unless(ret == 0x11A2,
      "test case 6 expected 0x11A2, got 0x%04X", ret);
  */
}
END_TEST

START_TEST(test_callbacks)
{
  /* Start with no callbacks registered.  */
  sbp_clear_callbacks();

  fail_unless(sbp_find_callback(0x1234) == 0,
      "sbp_find_callback should return NULL if no callbacks registered");

  fail_unless(sbp_register_callback(0x2233, &test_callback, 0) == SBP_NULL_ERROR,
      "sbp_register_callback should return an error if node is NULL");

  /* Add a first callback. */

  static sbp_msg_callbacks_node_t n;

  fail_unless(sbp_register_callback(0x2233, 0, &n) == SBP_NULL_ERROR,
      "sbp_register_callback should return an error if cb is NULL");

  fail_unless(sbp_register_callback(0x2233, &test_callback, &n) == SBP_OK,
      "sbp_register_callback should return success if everything is groovy");

  fail_unless(sbp_register_callback(0x2233, &test_callback, &n)
        == SBP_CALLBACK_ERROR,
      "sbp_register_callback should return SBP_CALLBACK_ERROR if a callback "
      "of the same type is already registered");

  fail_unless(sbp_find_callback(0x1234) == 0,
      "sbp_find_callback should return NULL if callback not registered");

  fail_unless(sbp_find_callback(0x2233) == &test_callback,
      "sbp_find_callback didn't return the correct callback function pointer");

  /* Add a second callback. */

  static sbp_msg_callbacks_node_t m;

  fail_unless(sbp_register_callback(0x1234, &test_callback2, &m) == SBP_OK,
      "sbp_register_callback should return success if everything is groovy (2)");

  fail_unless(sbp_find_callback(0x2233) == &test_callback,
      "sbp_find_callback didn't return the correct callback function pointer (2)");

  fail_unless(sbp_find_callback(0x1234) == &test_callback2,
      "sbp_find_callback didn't return the correct callback function pointer (3)");

  fail_unless(sbp_register_callback(0x1234, &test_callback, &n)
        == SBP_CALLBACK_ERROR,
      "sbp_register_callback should return SBP_CALLBACK_ERROR if a callback "
      "of the same type is already registered (2)");

  fail_unless(sbp_find_callback(0x7788) == 0,
      "sbp_find_callback should return NULL if callback not registered (2)");

  /* Clear all the registered callbacks and check they can no longer be found. */
  sbp_clear_callbacks();

  fail_unless(sbp_find_callback(0x1234) == 0,
      "sbp_find_callback should return NULL if no callbacks registered (2)");

  fail_unless(sbp_find_callback(0x2233) == 0,
      "sbp_find_callback should return NULL if no callbacks registered (3)");
}
END_TEST

Suite* sbp_suite(void)
{
  (void)test_sbp_send_message;

  Suite *s = suite_create("SBP");

  TCase *tc_core = tcase_create("Core");

  /* Clear callbacks before and after every test. */
  tcase_add_checked_fixture(tc_core, sbp_clear_callbacks, sbp_clear_callbacks);

  tcase_add_test(tc_core, test_callbacks);
  /*tcase_add_test(tc_core, test_sbp_send_message);*/

  suite_add_tcase(s, tc_core);

  return s;
}

