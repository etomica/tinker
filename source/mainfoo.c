#include <sys/time.h>
#include <stdlib.h>

#include <stdio.h>

long long getTime() {
  struct timeval t;
  gettimeofday(&t, NULL);
  long long tmus = (long long)(t.tv_sec) * 1000000LL + (long long)(t.tv_usec);
  return 1000LL*tmus;
}


double foo();

int main(int argc, char** argv) {
  long t1 = getTime();
  long n = 10000;
  if (argc > 1) {
    n = atol(argv[1]);
  }
  for (int i=0; i<n; i++) {
    foo();
  }
  long t2 = getTime();
  printf("%10.5f\n", (t2-t1)/(double)n/1e9);
}
