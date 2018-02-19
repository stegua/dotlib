/**
* @fileoverview Copyright (c) 2017-2018, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#pragma once

#include "OT_BasicTypes.h"

/**
* @brief Read images from text file and return corresponfing matrix
*/
histogram_t read_image(size_t n, std::string filename);