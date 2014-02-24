
#include <stdio.h>
#include <check.h>

#include <sbp.h>


int DUMMY_MEMORY_FOR_CALLBACKS = 0xdeadbeef;
int DUMMY_MEMORY_FOR_IO = 0xdead0000;

u32 dummy_wr = 0;
u32 dummy_rd = 0;
u8 dummy_buff[1024];
void* last_io_context;

void dummy_reset()
{
  dummy_rd = dummy_wr = 0;
  memset(dummy_buff, 0, sizeof(dummy_buff));
}

u32 dummy_write(u8 *buff, u32 n, void* context)
{
 last_io_context = context;
 u32 real_n = n;//(dummy_n > n) ? n : dummy_n;
 memcpy(dummy_buff + dummy_wr, buff, real_n);
 dummy_wr += real_n;
 return real_n;
}

u32 dummy_read(u8 *buff, u32 n, void* context)
{
 last_io_context = context;
 u32 real_n = n;//(dummy_n > n) ? n : dummy_n;
 memcpy(buff, dummy_buff + dummy_rd, real_n);
 dummy_rd += real_n;
 return real_n;
}

u32 dummy_read_single_byte(u8 *buff, u32 n, void* context)
{
  (void)n;
  last_io_context = context;
  memcpy(buff, dummy_buff + dummy_rd, 1);
  dummy_rd += 1;
  return 1;
}

u32 dummy_write_single_byte(u8 *buff, u32 n, void* context)
{
  (void)n;
  last_io_context = context;
  memcpy(dummy_buff + dummy_wr, buff, 1);
  dummy_wr += 1;
  return 1;
}

void printy_callback(u16 sender_id, u8 len, u8 msg[], void* context)
{
  (void)context;
  printf("MSG id: 0x%04X, len: %d\n", sender_id, len);
  if (len > 0) {
    for (u8 i=0; i<len; i++)
      printf("0x%02X ", msg[i]);
  }
  printf("\n");
}

u32 n_callbacks_logged;
u16 last_sender_id;
u8 last_len;
u8 last_msg[256];
void* last_context;


void logging_reset()
{
  n_callbacks_logged = 0;
  last_context = 0;
  memset(last_msg, 0, sizeof(last_msg));
}

void logging_callback(u16 sender_id, u8 len, u8 msg[], void* context)
{
  n_callbacks_logged++;
  last_sender_id = sender_id;
  last_len = len;
  last_context = context;
  memcpy(last_msg, msg, len);
  
  /*printy_callback(sender_id, len, msg);*/
}

void test_callback(u16 sender_id, u8 len, u8 msg[], void* context)
{
  /* Do nothing. */
  (void)sender_id;
  (void)len;
  (void)msg;
  (void)context;
}
void test_callback2(u16 sender_id, u8 len, u8 msg[], void* context)
{
  /* Do nothing. */
  (void)sender_id;
  (void)len;
  (void)msg;
  (void)context;
}

START_TEST(test_sbp_process)
{
  /* TODO: Tests with different read function behaviour. */

  sbp_state_t s;
  sbp_state_init(&s);
  sbp_state_set_io_context(&s, &DUMMY_MEMORY_FOR_IO);

  static sbp_msg_callbacks_node_t n;
  static sbp_msg_callbacks_node_t n2;

  sbp_register_callback(&s, 0x2269, &logging_callback, &DUMMY_MEMORY_FOR_CALLBACKS, &n);

  u8 test_data[] = { 0x01, 0x02, 0x03, 0x04 };

  dummy_reset();
  sbp_send_message(&s, 0x2269, 0x42, sizeof(test_data), test_data, &dummy_write);

  while (dummy_rd < dummy_wr) {
    fail_unless(sbp_process(&s, &dummy_read) >= SBP_OK,
        "sbp_process threw an error!");
  }

  fail_unless(n_callbacks_logged == 1,
      "one callback should have been logged");
  fail_unless(last_sender_id == 0x42,
      "sender_id decoded incorrectly");
  fail_unless(last_len == sizeof(test_data),
      "len decoded incorrectly");
  fail_unless(memcmp(last_msg, test_data, sizeof(test_data))
        == 0,
      "test data decoded incorrectly");
  fail_unless(last_context == &DUMMY_MEMORY_FOR_CALLBACKS, 
      "context pointer incorrectly passed");

  sbp_register_callback(&s, 0x2270, &logging_callback, &DUMMY_MEMORY_FOR_CALLBACKS, &n2);
  fail_unless(sbp_find_callback(&s, 0x2270) != 0,
    "second callback not found");

  logging_reset();
  sbp_send_message(&s, 0x2269, 0x4243, 0, 0, &dummy_write);

  fail_unless(last_io_context == &DUMMY_MEMORY_FOR_IO,
      "io context pointer incorrectly passed");

  last_io_context = 0;

  while (dummy_rd < dummy_wr) {
    fail_unless(sbp_process(&s, &dummy_read) >= SBP_OK,
        "sbp_process threw an error! (2)");
  }

  fail_unless(last_io_context == &DUMMY_MEMORY_FOR_IO,
      "io context pointer incorrectly passed");


  fail_unless(n_callbacks_logged == 1,
      "one callback should have been logged (2)");
  fail_unless(last_sender_id == 0x4243,
      "sender_id decoded incorrectly (2)");
  fail_unless(last_len == 0,
      "len decoded incorrectly (2)");

  logging_reset();
  sbp_send_message(&s, 0x22, 0x4243, 0, 0, &dummy_write);

  s8 ret = 0;
  while (dummy_rd < dummy_wr) {
    ret |= sbp_process(&s, &dummy_read);
  }

  fail_unless(ret == SBP_OK_CALLBACK_UNDEFINED,
      "sbp_process should have returned SBP_OK_CALLBACK_UNDEFINED "
      "if no cb was registered for that message type");

  fail_unless(n_callbacks_logged == 0,
      "no callbacks should have been logged");

  u8 awesome_message[] = {0x55, 0x33, 0x22, 0x77, 0x66,
                          0x02, 0x22, 0x33, 0x8A, 0x33};
  logging_reset();
  dummy_reset();
  dummy_rd = 0;
  dummy_wr = sizeof(awesome_message);
  memcpy(dummy_buff, awesome_message, sizeof(awesome_message));

  static sbp_msg_callbacks_node_t m;
  sbp_register_callback(&s, 0x2233, &logging_callback, 0, &m);

  while (dummy_rd < dummy_wr) {
    fail_unless(sbp_process(&s, &dummy_read) >= SBP_OK,
        "sbp_process threw an error! (3)");
  }

  fail_unless(n_callbacks_logged == 1,
      "one callback should have been logged (3)");
  fail_unless(last_sender_id == 0x6677,
      "sender_id decoded incorrectly (3)");
  fail_unless(last_len == 2,
      "len decoded incorrectly (3)");
  fail_unless(memcmp(last_msg, &awesome_message[6], 2)
        == 0,
      "test data decoded incorrectly (3)");

  awesome_message[4] = 0xAA;
  logging_reset();
  dummy_reset();
  dummy_rd = 0;
  dummy_wr = sizeof(awesome_message);
  memcpy(dummy_buff, awesome_message, sizeof(awesome_message));

  ret = 0;
  while (dummy_rd < dummy_wr) {
    ret |= sbp_process(&s, &dummy_read);
  }

  fail_unless(ret == SBP_CRC_ERROR,
      "sbp_process should have returned SBP_CRC_ERROR "
      "for malformed message");

  fail_unless(n_callbacks_logged == 0,
      "no callbacks should have been logged (2)");

  /* Test sbp_process with a one-byte-at-a-time read process */

  u8 awesome_message2[] = {0x55, 0x33, 0x22, 0x77, 0x66,
                          0x02, 0x22, 0x33, 0x8A, 0x33};
  logging_reset();
  dummy_reset();
  sbp_clear_callbacks(&s);
  dummy_rd = 0;
  dummy_wr = sizeof(awesome_message2);
  memcpy(dummy_buff, awesome_message2, sizeof(awesome_message2));

  static sbp_msg_callbacks_node_t p;
  sbp_register_callback(&s, 0x2233, &logging_callback, 0, &p);

  while (dummy_rd < dummy_wr) {
    fail_unless(sbp_process(&s, &dummy_read_single_byte) >= SBP_OK,
        "sbp_process threw an error! (3)");
  }

  fail_unless(n_callbacks_logged == 1,
      "one callback should have been logged (3)");
  fail_unless(last_sender_id == 0x6677,
      "sender_id decoded incorrectly (3)");
  fail_unless(last_len == 2,
      "len decoded incorrectly (3)");
  fail_unless(memcmp(last_msg, &awesome_message2[6], 2)
        == 0,
      "test data decoded incorrectly (3)");


  /* Test sbp_process with a one-byte-at-a-time read that starts with garbage */

  u8 crappy_then_awesome_message[] = {0x99, 0x88, 0x77, 0x66, 0x55, 0x33, 0x22, 0x77, 0x66,
                          0x02, 0x22, 0x33, 0x8A, 0x33};
  logging_reset();
  dummy_reset();
  sbp_clear_callbacks(&s);
  dummy_rd = 0;
  dummy_wr = sizeof(crappy_then_awesome_message);
  memcpy(dummy_buff, crappy_then_awesome_message, sizeof(crappy_then_awesome_message));

  static sbp_msg_callbacks_node_t q;
  sbp_register_callback(&s, 0x2233, &logging_callback, 0, &q);

  while (dummy_rd < dummy_wr) {
    fail_unless(sbp_process(&s, &dummy_read_single_byte) >= SBP_OK,
        "sbp_process threw an error! (3)");
  }

  fail_unless(n_callbacks_logged == 1,
      "one callback should have been logged (3)");
  fail_unless(last_sender_id == 0x6677,
      "sender_id decoded incorrectly (3)");
  fail_unless(last_len == 2,
      "len decoded incorrectly (3)");
  fail_unless(memcmp(last_msg, &crappy_then_awesome_message[10], 2)
        == 0,
      "test data decoded incorrectly (3)");

}
END_TEST

START_TEST(test_sbp_send_message)
{
  /* TODO: Tests with different write function behaviour. */

  sbp_state_t s;
  sbp_state_init(&s);


  u8 smsg[] = { 0x22, 0x33 };

  fail_unless(sbp_send_message(&s, 0x2233, 0x4455, 0, smsg, 0) == SBP_NULL_ERROR,
      "sbp_send_message should return an error if write is NULL");

  dummy_reset();
  fail_unless(sbp_send_message(&s, 0x2233, 0x4455, 1, 0, &dummy_write)
        == SBP_NULL_ERROR,
      "sbp_send_message should return an error if payload is NULL and len != 0");

  dummy_reset();
  fail_unless(sbp_send_message(&s, 0x2233, 0x4455, 0, 0, &dummy_write)
        == SBP_OK,
      "sbp_send_message should return OK if payload is NULL and len == 0");

  u8 zero_len_message[] = {0x55, 0x33, 0x22, 0x55, 0x44, 0x00, 0x2C, 0x4C};

  fail_unless(memcmp(dummy_buff, zero_len_message, sizeof(zero_len_message))
        == 0,
      "sbp_send_message encode error for len = 0");

  dummy_reset();
  sbp_send_message(&s, 0x2233, 0x6677, sizeof(smsg), smsg, &dummy_write);

  u8 awesome_message[] = {0x55, 0x33, 0x22, 0x77, 0x66,
                          0x02, 0x22, 0x33, 0x8A, 0x33};

  fail_unless(memcmp(dummy_buff, awesome_message, sizeof(awesome_message))
        == 0,
      "sbp_send_message encode error for test message");
}
END_TEST

START_TEST(test_callbacks)
{

  sbp_state_t s;
  sbp_state_init(&s);

  /* Start with no callbacks registered.  */
  sbp_clear_callbacks(&s);

  fail_unless(sbp_find_callback(&s, 0x1234) == 0,
      "sbp_find_callback should return NULL if no callbacks registered");

  fail_unless(sbp_register_callback(&s, 0x2233, &test_callback, 0, 0) == SBP_NULL_ERROR,
      "sbp_register_callback should return an error if node is NULL");

  /* Add a first callback. */

  static sbp_msg_callbacks_node_t n;

  int NUMBER = 42;

  fail_unless(sbp_register_callback(&s, 0x2233, 0, 0, &n) == SBP_NULL_ERROR,
      "sbp_register_callback should return an error if cb is NULL");

  fail_unless(sbp_register_callback(&s, 0x2233, &test_callback, &NUMBER, &n) == SBP_OK,
      "sbp_register_callback should return success if everything is groovy");

  fail_unless(sbp_register_callback(&s, 0x2233, &test_callback, 0, &n)
        == SBP_CALLBACK_ERROR,
      "sbp_register_callback should return SBP_CALLBACK_ERROR if a callback "
      "of the same type is already registered");

  fail_unless(sbp_find_callback(&s, 0x1234) == 0,
      "sbp_find_callback should return NULL if callback not registered");

  fail_unless(sbp_find_callback(&s, 0x2233) == &n,
      "sbp_find_callback didn't return the correct callback node pointer");

  fail_unless(sbp_find_callback(&s, 0x2233)->context == &NUMBER,
      "sbp_find_callback didn't return the correct context pointer");

  /* Add a second callback. */

  static sbp_msg_callbacks_node_t m;

  int NUMBER2 = 84;

  fail_unless(sbp_register_callback(&s, 0x1234, &test_callback2, &NUMBER2, &m) == SBP_OK,
      "sbp_register_callback should return success if everything is groovy (2)");

  fail_unless(sbp_find_callback(&s, 0x2233) == &n,
      "sbp_find_callback didn't return the correct callback function pointer (2)");

  fail_unless(sbp_find_callback(&s, 0x2233)->context == &NUMBER,
      "sbp_find_callback didn't return the correct context pointer");

  fail_unless(sbp_find_callback(&s, 0x1234) == &m,
      "sbp_find_callback didn't return the correct callback function pointer (3)");

  fail_unless(sbp_find_callback(&s, 0x1234)->context == &NUMBER2,
      "sbp_find_callback didn't return the correct context pointer");

  fail_unless(sbp_register_callback(&s, 0x1234, &test_callback, 0, &n)
        == SBP_CALLBACK_ERROR,
      "sbp_register_callback should return SBP_CALLBACK_ERROR if a callback "
      "of the same type is already registered (2)");

  fail_unless(sbp_find_callback(&s, 0x7788) == 0,
      "sbp_find_callback should return NULL if callback not registered (2)");

  /* Clear all the registered callbacks and check they can no longer be found. */
  sbp_clear_callbacks(&s);

  fail_unless(sbp_find_callback(&s, 0x1234) == 0,
      "sbp_find_callback should return NULL if no callbacks registered (2)");

  fail_unless(sbp_find_callback(&s, 0x2233) == 0,
      "sbp_find_callback should return NULL if no callbacks registered (3)");

}
END_TEST

Suite* sbp_suite(void)
{
  Suite *s = suite_create("SBP");

  TCase *tc_core = tcase_create("Core");


  tcase_add_test(tc_core, test_callbacks);
  tcase_add_test(tc_core, test_sbp_send_message);
  tcase_add_test(tc_core, test_sbp_process);

  suite_add_tcase(s, tc_core);

  return s;
}

