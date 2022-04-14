
#include "etomica.h"

double foo();

JNIEXPORT jdouble JNICALL Java_etomica_potential_amoeba_PotentialTinker_tinkerEnergy
  (JNIEnv *, jobject) {
    return foo();
}

JNIEXPORT jdouble JNICALL Java_etomica_potential_amoeba_PotentialTinker_tinkerEnergyN
  (JNIEnv *, jobject, jlong n) {
    double rv = 0;
    for (long i = 0; i < n; i++) {
      rv = foo();
    }
    return rv;
}
