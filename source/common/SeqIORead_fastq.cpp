/*
 * SeqIORead_fastq.cpp
 *
 *  Created on: 06/Feb/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 */

#include <ctype.h>
#include "SeqIORead_fastq.h"
#include "common.h"
#include "bio_exceptions.h"

namespace Bio {

/****************************************************************************************************************/
ReadSequence* SeqIORead_fastq::next_seq()
{
	string display_id, quality, sequence, desc;

	// Get header line
	this->getline(true);

	// If reached EOF - return
	if(this->m_next_line[0] == 0)
		return NULL;

	display_id = this->m_next_line;
	chomp(display_id);
	if(display_id.at(0) != '@')
		throw Bio::Bad_file(this->file_name(), "Unexpected first line for sequence in fastq file (should begin with @)", __PRETTY_FUNCTION__, this->line_num());
	display_id.erase(0, 1);
	m_line_num++;

	// Get sequence
	this->getline(true);
	// If reached EOF in the middle of an entry - error
	if(this->m_next_line[0] == 0)
		throw Bio::Bad_file(this->file_name(), "Unexpected fastq file end-of-file after sequence name line (@)", __PRETTY_FUNCTION__, this->line_num());
	sequence = this->m_next_line;
	chomp(sequence);
	m_line_num++;

	// Get description line
	this->getline(true);
	// If reached EOF in the middle of an entry - error
	if((this->m_next_line[0] == 0) || (this->m_next_line[0] != '+'))
		throw Bio::Bad_file(this->file_name(), "3rd line in sequence description should begin with a + sign but begins with something else", __PRETTY_FUNCTION__, this->line_num());
	desc = this->m_next_line;
	desc.erase(0, 1);
	m_line_num++;


	// Get quality scores
	this->getline(true);
	// If reached EOF in the middle of an entry - error
	if(this->m_next_line[0] == 0)
		throw Bio::Bad_file(this->file_name(), "Unexpected fastq file end-of-file after + line", __PRETTY_FUNCTION__, this->line_num());
	quality = this->m_next_line;
	chomp(quality);
	m_line_num++;

	return new ReadSequence(display_id, desc, sequence, quality);
}

}
