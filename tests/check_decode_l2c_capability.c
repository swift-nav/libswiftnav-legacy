#include <check.h>
#include <stdlib.h>
#include <stdio.h>

#include  <libswiftnav/l2c_capability.h>

#include "check_utils.h"

/* This calculates position SV configuration (4 bits) in Word 3,
 * where s is number (0..3) of SV inside the Word */
#define P3(s) (30-(12+s*4))
/* same as above, but for other Words, s in range (0..5)*/
#define P(s) (30-(4+s*4))

/* This gathers Word 3,
 *  c1..4 configuration for SV1-4, should be in range 0..15 */
#define W3(c1,c2,c3,c4) (c1<<P3(0)|c2<<P3(1)|c3<<P3(2)|c4<<P3(3))
/* same as above, but for other Words,
 * c1..6 config for SV1-6 _inside the Word_, range 0..15 */
#define W(c1,c2,c3,c4,c5,c6) \
        (c1<<P(0)|c2<<P(1)|c3<<P(2)|c4<<P(3)|c5<<P(4)|c6<<P(5))


START_TEST(l2c_update_test)
{
  struct {
    u32 frame4_words[8];
    u32 l2c_capability_expected;
  } test_case[] = {
      { /* received configuration -- all SV able to broadcast L2C stream */
        {
          W3(2,2,2,2),
          W(2,2,2,2,2,2),
          W(2,2,2,2,2,2),
          W(2,2,2,2,2,2),
          W(2,2,2,2,2,2),
          W(2,2,2,2,0,0),
          0,0
        },
        0xffffffff,
      },
      { /* received configuration -- even SV broadcasts, odd does not */
        {
          W3(5,2,7,3),
          W(0,4,5,4,7,2),
          W(0,3,5,4,7,2),
          W(0,2,5,3,7,2),
          W(0,3,5,2,7,4),
          W(0,2,5,4,7,2),
          0,0
        },
        0xaaaaaaaa,
      },
      { /* received configuration -- odd SV broadcasts, even does not */
        {
          W3(11,13,10,14),
          W(11,15,12,13,10,14),
          W(10,14,12,13,11,15),
          W(12,13,11,14,10,15),
          W(10,15,11,13,12,14),
          W(11,15,12,14,10,14),
          0,0
        },
        0x55555555,
      },
    };

  u32 l2c_cp = 0xdeadbeef;

  /* go through rest of the test cases */
  for(u8 i = 0; i < sizeof(test_case)/sizeof(test_case[0]); i++ ) {
    decode_l2c_capability(test_case[i].frame4_words, &l2c_cp);
    fail_unless(
      test_case[i].l2c_capability_expected == l2c_cp,
      "Test %u failed. Expected 0x%x, got: value 0x%x\n",
       i,
       test_case[i].l2c_capability_expected,
       l2c_cp
    );
  }


}
END_TEST

Suite* l2c_capability_test_suite(void)
{
  Suite *s = suite_create("L2C capability API check");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, l2c_update_test);
  suite_add_tcase(s, tc_core);

  return s;
}
