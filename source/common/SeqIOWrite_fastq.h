/*
 * SeqIOWrite_fastq.h
 *
 *  Created on: 06/Feb/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 *
 */

#ifndef SEQIOWRITE_FASTQ_H_
#define SEQIOWRITE_FASTQ_H_

#include "SeqIOWrite.h"
#include "ReadSequence.h"

using namespace std;

namespace Bio {

/****************************************************************************************************************
 * SeqIOWrite_fastq
 ****************************************************************************************************************/
class SeqIOWrite_fastq : public SeqIOWrite<ReadSequence> {
public:
	SeqIOWrite_fastq(string fastq_file) : SeqIOWrite(fastq_file) 					{}
	void								write_seq(const ReadSequence& seq);
};

}
#endif /* SEQIOWRITE_FASTQ_H_ */
