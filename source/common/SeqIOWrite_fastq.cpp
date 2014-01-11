/*
 * SeqIOWrite_fastq.cpp
 *
 *  Created on: 06/Feb/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 */

#include <ctype.h>
#include "SeqIOWrite_fastq.h"
#include "common.h"
#include "bio_exceptions.h"

namespace Bio {

/****************************************************************************************************************/
void SeqIOWrite_fastq::write_seq(const ReadSequence& seq)
{
        m_fp << '@' <<  seq.display_id() << '\n' <<  string(seq.seq()) << "\n+" << seq.desc() << '\n' << seq.quality() << endl;
        m_line_num += 4;
        m_seq_num++;
}

}
