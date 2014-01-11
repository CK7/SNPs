/*
 * ReadMappingReader.h
 *
 *  Created on: 05/Feb/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 *
 * This version is implemented using low-level C features (FILE*, pointers instead of ifstream and string) in order to 
 * increase speed. As a result this version takes ~40% less time than the first one.
 */

#ifndef READMAPPINGREADER_H_
#define READMAPPINGREADER_H_

#include "common.h"
#include "ReadMapping.h"
#include <stdio.h>

namespace Bio {

/***************************************************************************************************************************
 * ReadMappingReader
 * Base class for classes reading read mapping files of all formats
 ***************************************************************************************************************************/
class ReadMappingReader {
public:
	ReadMappingReader(string file);
	virtual 			~ReadMappingReader()		{if(m_fp!=NULL) {fclose(m_fp);}}
	virtual	ReadMapping*		next_mapping() = 0;
	void				restart();
	const string&			file_name() const		{return m_file_name;}
	bool				good() const			{return !feof(m_fp);} 
	double				percent_entries_read() const	{if(m_total_to_read==0) return 100.0; return (100000*m_total_read/m_total_to_read)/1000.0;}
protected:
	FILE*		m_fp;
	string		m_file_name;
	uint64_t	m_total_to_read;
	uint64_t	m_total_read;
};

/***************************************************************************************************************************
 * SAMReader
 * Implementation for the SAM format
 ***************************************************************************************************************************/
class SAMReader : public ReadMappingReader {
public:
	static const ReadMapping&	end()				{return _end;}
public:
	SAMReader(string sam_file) : ReadMappingReader(sam_file)
		{this->read_header();}
	// Will return the next entry. ReadMappingReader::end will be returned if reached the end of the file 
	ReadMapping*			next_mapping();
	// Currently does nothing, have to decide what to do
	void				read_header();
protected:
	static ReadMapping	_end;
};

}
#endif /* READMAPPINGREADER_H_ */
