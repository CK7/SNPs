/*
 * SeqIO_fastq.cpp
 *
 *  Created on: 31/Jan/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 */

#include <ctype.h>
#include "SeqIO_fastq.h"
#include "common.h"
#include "bio_exceptions.h"

namespace Bio {

/****************************************************************************************************************/
ReadSequence* SeqIORead_fastq::next_seq()
{
	if(!m_fp.good())
		return NULL;

//The following didn't work - check it out using a tiny file
	string display_id, quality, sequence, desc;
	while(m_fp.good() && (display_id.size()==0))
		getline(m_fp, display_id);
	if(!m_fp.good())
		return NULL;

	chomp(display_id);
	if(display_id.at(0) != '@')
		throw Bio::Bad_file(this->file_name(), "Unexpected first line for sequence in fastq file (should begin with @)", __PRETTY_FUNCTION__, this->line_num());
	display_id.erase(0, 1);
	m_line_num++;

	if(!m_fp.good())
		throw Bio::Bad_file(this->file_name(), "Unexpected fastq file end-of-file after sequence name line (@)", __PRETTY_FUNCTION__, this->line_num());

	getline(m_fp, sequence);
	chomp(sequence);
	m_line_num++;

	if(!m_fp.good())
		throw Bio::Bad_file(this->file_name(), "Unexpected fastq file end-of-file after sequence line", __PRETTY_FUNCTION__, this->line_num());

	getline(m_fp, desc);
	if(desc.at(0) != '+')
		throw Bio::Bad_file(this->file_name(), "3rd line in sequence description should begin with a + sign but begins with something else", __PRETTY_FUNCTION__, this->line_num());
	desc.erase(0, 1);

	m_line_num++;

	if(!m_fp.good())
		throw Bio::Bad_file(this->file_name(), "Unexpected fastq file end-of-file after + line", __PRETTY_FUNCTION__, this->line_num());

	getline(m_fp, quality);
	chomp(quality);
	m_line_num++;

	m_seq_num++;
	return new ReadSequence(display_id, desc, sequence, quality);
}


/****************************************************************************************************************/
void SeqIOWrite_fastq::write_seq(const ReadSequence& seq)
{
        m_fp << '@' <<  seq.display_id() << '\n' <<  string(seq.seq()) << "\n+" << seq.desc() << '\n' << seq.quality() << endl;
        m_line_num += 4;
        m_seq_num++;
}

}
