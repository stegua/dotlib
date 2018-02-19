#include "OT_Utils.h"

/**
* @brief Read images from text file and return corresponfing matrix
*/
histogram_t read_image(size_t n, std::string filename) {
   auto m = [&n](size_t x, size_t y) {
      return x * n + y;
   };

   histogram_t image(n * n, 0);
   image.shrink_to_fit();

   std::ifstream in_file(filename);

   if (!in_file) {
      fprintf(stdout, "Cannot open file %s.\n", filename);
      return image;
   }

   for (size_t i = 0; i < n; ++i) {
      size_t j = 0;
      std::string         line;
      std::getline(in_file, line);
      std::stringstream   lineStream(line);
      std::string         cell;

      while (std::getline(lineStream, cell, ',')) {
         image[m(i, j)] = stod(cell);
         ++j;
      }
   }

   in_file.close();

   return image;
}
