/*
 * ReadSequence.cpp
 *
 *  Created on: 05/Feb/2013
 *      Author: itai Sharon, itai.sharon@gmail.com
 */

#include "ReadSequence.h"

namespace Bio {

ReadSequence::Qscore_scheme 	ReadSequence::quality_offset = phred64_quals;
int				ReadSequence::max_score = 40;

}
