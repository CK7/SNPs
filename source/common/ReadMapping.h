/*
 * ReadMapping.h
 *
 *  Created on: 06/Feb/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 */

#ifndef READMAPPING_H_
#define READMAPPING_H_

#include "common.h"
#include "String.h"
#include <tr1/memory>
#include <iostream>
namespace Bio {

/***************************************************************************************************************************
 * ReadMapping
 * Stores information about the mapping of a read to a reference sequence.
 * Uses 0 offset (like all other classes).
 ***************************************************************************************************************************/
class ReadMapping {
public:
	class Illegal_mapping;
public:
	struct SNP;
	typedef enum {forward, reverse} Direction;
public:
	// ref_pos is expected to be a 0-offset position.
	ReadMapping(const char* SAM_line);
//	ReadMapping(string read_name, string ref_name, size_t ref_pos, Direction dir, string snp_CIGAR_string, const DNAString& dna_str, const string& quality_str) :
//		m_ref_name(ref_name), m_read_name(read_name), m_read(dna_str), m_ref_pos(ref_pos), m_dir(dir)
//		{this->determine_snps(dna_str, quality_str, snp_CIGAR_string, m_ref_pos);}
	virtual 		~ReadMapping()					{/*std::cout << "Destroying " << m_ref_name << std::endl;*/}
	const string&		ref_name() const				{return m_ref_name;}
	const string&		read_name() const				{return m_read_name;}
	const size_t		read_length() const				{return m_read.size();}
	const Bio::DNAString	mapped_read() const				{return m_read;}
	size_t			ref_pos() const					{return m_ref_pos;}
	Direction		direction() const				{return m_dir;}
	size_t			num_snps() const				{return m_snps.size();}
	// Index may be 0 .. (num_snps()-1)
	const SNP&		snp(size_t index) const				{return m_snps[index];}
	bool			unmapped() const				{return (m_flag & 4);}
	bool			multiple_hits() const				{return (m_flag & 0x100); }
protected:
	void 			determine_snps(const DNAString& dna_str, const string& quality_str, string CIGAR_indels, string CIGAR_mismatches, size_t ref_pos);
protected:
	string			m_ref_name;
	size_t			m_flag;
	string			m_read_name;
	Bio::DNAString		m_read;
	string			m_quality;
	size_t			m_ref_pos;
	vector<SNP>		m_snps;
	Direction		m_dir;
};

typedef std::tr1::shared_ptr<ReadMapping> ReadMappingPtr;

/****************************************************************************************************************/
struct ReadMapping::SNP {
	typedef enum {INSERTION, DELETION, MISMATCH} SNP_type;
	SNP_type		snp_type;
	size_t			ref_pos, read_pos;
	char			ref_char, read_char;
	unsigned long		read_quality;
	SNP(SNP_type _snp_type, size_t _ref_pos, size_t _read_pos, char _ref_char, char _read_char, unsigned long _read_quality) :
		snp_type(_snp_type), ref_pos(_ref_pos), read_pos(_read_pos), ref_char(_ref_char), read_char(_read_char), read_quality(_read_quality)
		{}
	SNP() : ref_pos(0), read_pos(0), ref_char('-'), read_char('-'), read_quality(0)
		{}
};

/****************************************************************************************************************/
class ReadMapping::Illegal_mapping : public exception {
	string msg;
public:
	Illegal_mapping(string mapping_str, string err_msg, const char* func)
		{
			msg = string("Fatal error:\t") + err_msg + string("\nMapping string:\t") + mapping_str + string("\n");
			msg += string("Function:  \t") + string(func) + string("\n");

		}
	~Illegal_mapping() throw()
		{}
	virtual const char* what() const throw()
		{return msg.c_str();}
};

}
#endif /* READMAPPING_H_ */
