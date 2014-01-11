/*
 * ScafSNPs.h
 *
 *  Created on: 06/Feb/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 */

#ifndef SCAFSNPS_H_
#define SCAFSNPS_H_

#include "common/common.h"
#include "common/bio_exceptions.h"
#include "common/String.h"
#include "common/ReadMapping.h"
#include <vector>

/*****************************************************************************************************************
 * ScafSNPs
 * This class keeps SNP information for a scaffold
 ****************************************************************************************************************/
class ScafSNPs {
public:
	class SNPInfo;
public:
	ScafSNPs(Bio::DNASequence& scaf);
	virtual ~ScafSNPs()												{}
	const ScafSNPs&			operator += (const Bio::ReadMapping& mapping);
	double  			scaf_coverage() const;
	size_t				size() const						{return m_snps.size();}
	const SNPInfo&			operator [] (unsigned int i) const;
	const Bio::DNASequence&		scaf() const						{return m_scaf;}
protected:
	void				update_m_snps();
protected:
	Bio::DNASequence&		m_scaf;
	vector<size_t>			m_read_starts;
	vector<size_t>			m_read_ends;
	// The use of mutable is ugly, I am aware of that. The reason I am using it is because m_snps may be updated
	// within print_snps, when I add the sum of non-SNPs to each position. This is done in order to save the time
	// it will take to update all relevant positions for each read and not just the SNPs. So I think it is better
	// to declare m_snps as a mutable than to make print_snps() non-const.
	mutable vector<SNPInfo>		m_snps;
	bool				m_snps_updated;
};

/*****************************************************************************************************************
 * ScafSNPs::SNPInfo
 * This class keeps SNP information for a single position in a scaffold
 ****************************************************************************************************************/
class ScafSNPs::SNPInfo {
public:
	class Illegal_char : public exception_base {
	public:
		Illegal_char(char c, const char* func) : exception_base("", func) {
			string s = string("Illegal char: ");
			s.append(1, c);
			msg = s + msg;
		}
	};
public:
	SNPInfo() : _A(0), _C(0), _G(0), _T(0), _N(0), _insertion(0), _deletion(0)			{}
	// The following functions refer to mismatches
	size_t	A() const		{return _A;}
	size_t	C() const		{return _C;}
	size_t	G() const		{return _G;}
	size_t	T() const		{return _T;}
	size_t	N() const		{return _N;}
	size_t	total() const	{return (_A+_C+_G+_T+_N);}
	
	void	add(char c, size_t amount=1)
	{
		switch(c) {
			case 'A': _A += amount; break;
			case 'C': _C += amount; break;
			case 'G': _G += amount; break;
			case 'T': _T += amount; break;
			case 'N':
			default:
				_N += amount; break;
		}
	}
	void	set(char c, size_t amount=0)
	{
		switch(c) {
			case 'A': _A = amount; break;
			case 'C': _C = amount; break;
			case 'G': _G = amount; break;
			case 'T': _T = amount; break;
			case 'N':
			default:
				_N = amount;
		}
	}
	// These refer to indels
	size_t	insertion() const		{return _insertion;}
	size_t	deletion() const		{return _deletion;}
	void	add_insertion(size_t i=1)	{_insertion += i;}
	void	add_deletion(size_t i=1)	{_deletion += i;}
protected:
	size_t	_A, _C, _G, _T, _N;
	size_t	_insertion, _deletion;
};

#endif /* SCAFSNPS_H_ */
