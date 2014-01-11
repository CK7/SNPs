/*
 * ReadSequence.h
 *
 *  Created on: 05/Feb/2013
 *      Author: itai Sharon, itai.sharon@gmail.com
 */

#ifndef READSEQUENCE_H_
#define READSEQUENCE_H_

#include "Sequence.h"
#include "String.h"

namespace Bio {

/****************************************************************************************************************
 * ReadSequence
 * Default quality score scheme is phred64_quals
 ****************************************************************************************************************/
class ReadSequence : public DNASequence {
public:
	typedef enum {phred33_quals = 33, phred64_quals = 64, solexa1_3_quals = 64} Qscore_scheme;
	static Qscore_scheme 	quality_offset;
	static int		max_score;
public:
	ReadSequence(string display_id, string desc, DNAString sequence, string quality_string) : 
		DNASequence(display_id, desc, sequence), m_quality(quality_string)
		{}
/*	ReadSequence(string display_id, DNAString sequence, string quality_string) : 
		DNASequence(display_id, "", sequence), m_quality(quality_string)
		{}*/
	ReadSequence(string display_id, string desc, DNAString sequence) : 
		DNASequence(display_id, desc, sequence)
		{m_quality.resize(sequence.size()); m_quality.assign(m_quality.size(), max_score+quality_offset);}
/*	ReadSequence(string display_id, DNAString sequence) : 
		DNASequence(display_id, "", sequence)
		{m_quality.resize(sequence.size()); m_quality.assign(m_quality.size(), max_score+quality_offset);}*/
	const string&	quality() const						{return m_quality;}
protected:
	string		m_quality;
};

}

#endif /* READSEQUENCE_H_ */
