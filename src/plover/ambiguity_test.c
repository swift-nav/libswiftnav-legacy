#include <stdio.h>
#include <single_diff.h>
#include <math.h>
#include <plover/ambiguity_test.h>
int hello_world(void)
{
  printf("hello world\n");
  return 0;
}
int assign_de_mtxP(int num_sats, sdiff_t sds[num_sats], double ref_ecef[3], double DE[(num_sats + -(1))][3])
{
  double e_0[3];
  for (int var1 = 0; var1 < 3; var1++) {
    double var5;
    var5 = sds[0].sat_pos[var1];
    double var9;
    var9 = -(ref_ecef[var1]);
    e_0[var1] = (var5 + var9);
  }
  double e_0_[3];
  for (int var2 = 0; var2 < 3; var2++) {
    double var7;
    double var12;
    var12 = 0.0;
    for (int var15 = 0; var15 < 3; var15++) {
      var12 = (var12 + (e_0[var15] * e_0[var15]));
    }
    var7 = sqrt(var12);
    e_0_[var2] = (e_0[var2] / var7);
  }
  for (int i = 0; i < (num_sats + -(1)); i++) {
    for (int var0 = 0; var0 < 3; var0++) {
      double var3;
      double var10;
      double var16;
      var16 = sds[i].sat_pos[var0];
      double var18;
      var18 = -(ref_ecef[var0]);
      var10 = (var16 + var18);
      double var14;
      double var17;
      var17 = 0.0;
      for (int var19 = 0; var19 < 3; var19++) {
        double var21;
        double var22;
        double var24;
        var24 = sds[i].sat_pos[var19];
        double var26;
        var26 = -(ref_ecef[var19]);
        var22 = (var24 + var26);
        double var23;
        double var25;
        var25 = sds[i].sat_pos[var19];
        double var27;
        var27 = -(ref_ecef[var19]);
        var23 = (var25 + var27);
        var21 = (var22 * var23);
        var17 = (var17 + var21);
      }
      var14 = sqrt(var17);
      var3 = (var10 / var14);
      double var6;
      var6 = -(e_0[var0]);
      DE[i][var0] = (var3 + var6);
    }
  }
  return 0;
}
