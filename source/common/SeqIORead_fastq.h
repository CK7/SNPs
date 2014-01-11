/*
 * SeqIO_fastq.h
 *
 *  Created on: 06/Feb/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 *
 */

#ifndef SEQIOREAD_FASTQ_H_
#define SEQIOREAD_FASTQ_H_

#include "SeqIORead.h"
#include "ReadSequence.h"

using namespace std;

namespace Bio {

/****************************************************************************************************************
 * SeqIORead_fastq
 ****************************************************************************************************************/
class SeqIORead_fastq : public SeqIORead<ReadSequence> {
public:
	SeqIORead_fastq(string fastq_file) : SeqIORead(fastq_file) 					{}
	ReadSequence*				 			next_seq();
};

}
#endif /* SEQIOREAD_FASTQ_H_ */
