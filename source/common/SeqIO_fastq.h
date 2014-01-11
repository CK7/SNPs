/*
 * SeqIO_fastq.h
 *
 *  Created on: 31/Jan/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 *
 *  Classes:
 *  - SeqIORead_fastq (no write version yet) is for reading sequences from a fastq file.
 */

#ifndef SEQIO_FASTQ_H_
#define SEQIO_FASTQ_H_

#include "SeqIO.h"
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

/****************************************************************************************************************
 * SeqIOWrite_fastq
 ****************************************************************************************************************/
class SeqIOWrite_fastq : public SeqIOWrite<ReadSequence> {
public:
	SeqIOWrite_fastq(string fastq_file) : SeqIOWrite(fastq_file) 					{}
	void								write_seq(const ReadSequence& seq);
};

}
#endif /* SEQIO_FASTQ_H_ */
