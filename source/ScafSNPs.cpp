/*
 * ScafSNPs.cpp
 *
 *  Created on: 06/Feb/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 */

#include "ScafSNPs.h"
#include "common/common.h"
#include <ostream>
#include <iostream>

/*****************************************************************************************************************/
ScafSNPs::ScafSNPs(Bio::DNASequence& scaf) : m_scaf(scaf), m_read_starts(scaf.seq().size()), m_read_ends(scaf.seq().size()),
	m_snps(scaf.seq().size()), m_snps_updated(false)
{
	Bio::DNAString::iterator end = scaf.seq().end();
	for(Bio::DNAString::iterator it = scaf.seq().begin(); it != end; it++) {
		if((*it != 'A') && (*it != 'C') && (*it != 'G') && (*it != 'T') && (*it != 'N')) {
			cerr << time_since_epoch() << "+++ Wrning: " << scaf.display_id() << " conatains non-A/C/G/T/N letters +++" << endl;
			break;
		}
	}
}
#include <string.h>
/*****************************************************************************************************************
 * ScafSNPs::operator +=
 * add information for SNPs on scaf. Note that the offset of coordinates stored in ReadMapping is 1 and here we 
 * translate them to 0.
 *****************************************************************************************************************/
const ScafSNPs& ScafSNPs::operator += (const Bio::ReadMapping& mapping)
{
#ifdef _DEBUG
	if(mapping.ref_pos() >= m_read_starts.size()) {
		char info[512];
		sprintf(info, "\n\tRead:\t%s (%d bps)\n\tSNP #:\t%d\n\tPosition:\t%d\n\tm_read_starts:\t%d\n\tref_pos exceeds m_read_starts size", 
				mapping.read_name().c_str(), int(mapping.read_length()), int(i), int(snp.ref_pos), int(m_read_starts.size()));
		throw Software_bug(info, "", __PRETTY_FUNCTION__);
	}
#endif
	m_snps_updated = false;
	m_read_starts[mapping.ref_pos()]++;
	size_t end = mapping.ref_pos()+mapping.read_length();
	if(end<m_read_ends.size())
		m_read_ends[end]++;
	// Else - we don't need to mark the end
	size_t	ndeleted = 0;
	for(size_t i=0; i<mapping.num_snps(); i++) {
		try {
			const Bio::ReadMapping::SNP& snp = mapping.snp(i);

			// Sanity check
			if((snp.ref_pos < mapping.ref_pos()) || (snp.ref_pos >= (mapping.ref_pos()+mapping.read_length()+ndeleted))) {
				char info[512];
				sprintf(info, "\n\tRead:\t%s (%d bps)\n\tSNP #:\t%d\n\tPosition:\t%d\n\tRange:\t(%d, %d)\n\treported SNP is outside of read mapping", 
						mapping.read_name().c_str(), int(mapping.read_length()), int(i), int(snp.ref_pos), int(mapping.ref_pos()), 
						int(mapping.ref_pos()+mapping.read_length()+ndeleted));
				throw Software_bug(info, "", __PRETTY_FUNCTION__);
			}

			// Sanity check
			if(snp.ref_pos >= m_snps.size()) {
				char info[512];
				sprintf(info, "\n\tReference position (%d) is out of range (%d)", int(snp.ref_pos), int(m_snps.size()-1));
				throw Software_bug(info, "", __PRETTY_FUNCTION__);
			}
			if(snp.snp_type == Bio::ReadMapping::SNP::INSERTION) {
				m_snps[snp.ref_pos].add_insertion();
			}
			else if(snp.snp_type == Bio::ReadMapping::SNP::DELETION) {
				m_snps[snp.ref_pos].add_deletion();
				ndeleted++;
			}
			else {
				m_snps[snp.ref_pos].add(snp.read_char);
			}
		}
		catch(SNPInfo::Illegal_char c) {
			char info[512];
			sprintf(info, "\n\tRead:\t%s (%d bps)\n\tSNP #:\t%d\n\tChar:\t%s", mapping.read_name().c_str(), int(mapping.read_length()),
					int(i), c.what());
			throw Software_bug(info, "", __PRETTY_FUNCTION__);
		}
	}

	return *this;
}

/*****************************************************************************************************************/
void	ScafSNPs::update_m_snps()
{
	int cvg = 0;
	for(size_t i=0; i<m_snps.size(); i++) {
		SNPInfo&	snp_info = m_snps[i];

		snp_info.set(m_scaf.seq().get(i), 0);

		cvg += m_read_starts[i];
		cvg -= m_read_ends[i];

		// Sanity check
		if(cvg < 0) {
			char info[512];
			sprintf(info, "\n\tScaf:\t%s (%d bps)\n\tSNP #:\t%d\n\tError:\tCoverage has gone below 0: %d", m_scaf.display_id().c_str(), int(m_scaf.seq().size()),
					int(i), cvg);
			throw Software_bug(info, "", __PRETTY_FUNCTION__);			
		}

		// Note that the following can be done several times and along the process of reding mappings
		try {
			snp_info.add(m_scaf.seq().get(i), cvg-snp_info.total());
		}
		catch(SNPInfo::Illegal_char c) {
			char info[512];
			sprintf(info, "\n\tScaf:\t%s (%d bps)\n\tSNP #:\t%d\n\tChar:\t%s", m_scaf.display_id().c_str(), int(m_scaf.seq().size()),
					int(i), c.what());
			throw Software_bug(info, "", __PRETTY_FUNCTION__);
		}
	}
	m_snps_updated = true;
}

/*****************************************************************************************************************/
const ScafSNPs::SNPInfo& ScafSNPs::operator [] (unsigned int i) const
{
	if(!m_snps_updated)
		const_cast<ScafSNPs*>(this)->update_m_snps();

	return m_snps[i];
}

/*****************************************************************************************************************/
double	ScafSNPs::scaf_coverage() const
{
	if(!m_snps_updated)
		const_cast<ScafSNPs*>(this)->update_m_snps();

	if(m_snps.size() <= 1)
		return 0;

	// Compute coverage, perform sanity checks
	double total_cvg = 0;
	for(size_t i=1; i<m_snps.size(); i++) {
		total_cvg += m_snps[i].total();
	}	
	total_cvg /= (m_snps.size()-1);

	return total_cvg;
}
