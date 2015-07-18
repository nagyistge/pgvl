#include <stdio.h>
#include "config.h"
#include <Image.h>
#include <ImageProcessing.h>
#include <time.h>

int main() {
   Image<uint8_t> lena(TEST_IMAGE_DIR "lena_gray.pgm");
   Image<uint8_t> lpf(lena.rows(), lena.cols(), lena.channels());

   volatile int numLoops = 5;
   clock_t beg = clock();
   for(int i = 0; i < numLoops; ++i)
      lowpassFilter(lpf, lena, 4);
   clock_t end = clock();

   printf("Time per call: %.1fms\n", static_cast<float>(end-beg)*1000.f/CLOCKS_PER_SEC/numLoops);
   return 0;
}