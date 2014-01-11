//================================================================================================================
// Name        : SNPs.cpp
// Author      : Itai Sharon
// Version     : 1.02, 31/Dec/2013
// Description : Computes the number of calls for of A, C, G, T, N, insertions and deletion in every position. on every sequence and generates a report file in
//             : the format
//             : ><scaf-name> /size=<size> /coverage=<coverage> /gc=<%G+C> /Ns=<# of Ns>
//             : <pos>	<overage>	<scaf-char>	<# A's>	<# C's>	<# G's>	<# T's>	<# N's>	<# Insertions>	<# Deletions>	<additional info>
//             : where <additional info> can be
//             : X -> Y, if X is the reference char but Y appears more times
//             : SNP if the number of calls of the consensus char is less than (1-<SNP threshold>) of the A, C, G and T (but not N) calls
//             : INSERTION if the number of insertion calls is more than <SNP threshold> of A, C, G, T AND N calls for the position.
//             : DELETION if the number of deletion calls if more than <SNP threshold> of A, C, G, T, N AND deletion call for the position  
//             : Insertions refer to insertions coming right after the position in which they are recognized. Insertions can be of 1 or more
//             : bases, in this version only the presence of an insertion is reported.
//================================================================================================================

#define VERSION 	1.02
#define DATE		"12/31/13"

#include <iostream>
#include <set>
#include <stdlib.h>
#include <typeinfo>
#include <sys/stat.h>
#include "common/common.h"
#include "common/SeqIORead_fasta.h"
#include "common/SequenceDB.h"
#include "common/ReadSequence.h"
#include "common/ReadMappingReader.h"
#include "ScafSNPs.h"

using namespace std;

/****************************************************************************************************************/
void 		usage(const char* prog_name);
bool 		read_arguments(int argc, const char* argv[], set<string>& sam_files, string& assembly_file, string& output_file, bool& verbose, double& snp_threshold);
void 		do_mapping(map<string, string>& pe_read2sam, string assembly_file, ostream& os, string root_dir);
void 		print_SNPs(map<string, ScafSNPs*>& scaf_snps, ostream& ossnps, double snp_threshold);

/****************************************************************************************************************/
int main(int argc, const char* argv[])
{
	try {
		// These will hold user-defined parameters
		set<string>		sam_files;		// sam files
		string			assembly_file;		// Path of assembly file
		string			output_profile;		// Path for output directory
		map<string, ScafSNPs*>	scaf_snps;
		bool			verbose = true;
		double			snp_threshold = 0.2;

		// Path for output files
		// Read arguments. If the user just wanted to know what the command line looks like - finish here
		if(!read_arguments(argc, argv, sam_files, assembly_file, output_profile, verbose, snp_threshold))
			return -1;

		// Open the output file
		ofstream	ossnps(output_profile.c_str());
		if(ossnps.fail()) {
			throw invalid_argument(string("Failed to write to ") + output_profile);
		}

		// Load the assembly file
		if(verbose)
			cerr << time_since_epoch() << "Loading assembly file (" << assembly_file << ")" << endl;

		Bio::SeqIORead_fasta<Bio::DNASequence>			reader_fasta(assembly_file);
		Bio::SequenceDB<Bio::DNASequence>			assembly_db(reader_fasta);
		for(Bio::SequenceDB<Bio::DNASequence>::iterator it=assembly_db.begin(); it!=assembly_db.end(); it++)
			scaf_snps[it->display_id()] = new ScafSNPs(*it);

		if(verbose)
			cerr << time_since_epoch() << "ok, " << assembly_db.size() << " sequences read" << endl;

		// Start constructing libraries info
		for(set<string>::const_iterator it=sam_files.begin(); it!=sam_files.end(); it++) 
		{
			if(verbose)
				cerr << time_since_epoch() << "*** Reading " << *it << " ***" << endl;
			Bio::SAMReader		reader(*it);
			size_t			num_read = 0;
			size_t			percent_read = 0;

			if(verbose) {
				cerr << time_since_epoch();
				print_percent(0, 0);
			}

			// Go over all mappings, collect the first mapping (if multiple exist) for each read
			while(reader.good()) {
				Bio::ReadMapping* mapping = reader.next_mapping();
				if(mapping->read_name() == Bio::SAMReader::end().read_name())
					break;

				if(((++num_read)%(1000) == 0) && (size_t(reader.percent_entries_read()) > percent_read)) {
					if(verbose)
						print_percent(percent_read, reader.percent_entries_read());
					percent_read = reader.percent_entries_read();
				}
				if(mapping->unmapped()) {
					continue;
				}
				// If multiple mappings exist we'll take the first one only
				if(mapping->multiple_hits()) {
					continue;
				}
				if(scaf_snps.find(mapping->ref_name()) != scaf_snps.end()) 
				{
					*(scaf_snps[mapping->ref_name()]) += *mapping;
				}

				delete(mapping);
			}
			if(verbose) {
				print_percent(percent_read, 100);
				cerr << endl;
				cerr << time_since_epoch() << "Finished, read " << num_read << " legal reads" << endl;
			}
		}
		// Now write output
		print_SNPs(scaf_snps, ossnps, snp_threshold);

		ossnps.close();
	}
	catch(const invalid_argument& e) {
		cerr << endl << "Invalid argument: " << e.what() << endl << endl;
		return -1;
	}
	catch(exception& e) {
		cerr 	<< endl << "\n*** Exception occurred ***" << endl << "Exception:\t" << typeid(e).name() << endl << "Error:    \t" << e.what()
				<< endl << endl;
		return -1;
	}
	catch(...) {
		cerr 	<< endl << "\n*** Unrecognized Exception occurred ***" << endl << endl;
		return -1;
	}
	return 0;
}

/****************************************************************************************************************/
void print_scaf_SNPs(ostream& os, const ScafSNPs& scaf_snps, ostream& ossnps, double snp_threshold);
void print_SNPs(map<string, ScafSNPs*>& scaf_snps, ostream& ossnps, double snp_threshold)
{
	for(map<string, ScafSNPs*>::const_iterator it=scaf_snps.begin(); it!=scaf_snps.end(); it++) {
		print_scaf_SNPs(ossnps, *(it->second), ossnps, snp_threshold);
		ossnps << endl;
		delete(it->second);
	}
}

/*****************************************************************************************************************/
void print_scaf_SNPs(ostream& os, const ScafSNPs& scaf_snps, ostream& ossnps, double snp_threshold)
{
	double 			total_cvg = (100.0*scaf_snps.scaf_coverage())/100;
	unsigned int		coverage_threshold = 10;
	const Bio::DNASequence&	scaf = scaf_snps.scaf();

	ossnps 	<< ">" << scaf.display_id() << " /size=" << scaf.seq().size() << " /coverage=" << total_cvg  << " /gc=" << Bio::gc(scaf.seq()) << " /Ns=" << Bio::Ns(scaf.seq()) << endl;
	for(size_t i=0; i<scaf_snps.size(); i++) {
		const ScafSNPs::SNPInfo&	snp_info = scaf_snps[i];
		char				abundant = (snp_info.A() > snp_info.C())? 'A' : 'C';
		size_t				n_abundant = (snp_info.A() > snp_info.C())? snp_info.A() : snp_info.C();

		if(snp_info.G() > n_abundant) {
			abundant = 'G';
			n_abundant = snp_info.G();
		}
		if(snp_info.T() > n_abundant) {
			abundant = 'T';
			n_abundant = snp_info.T();
		}

		ossnps 	<< (i+1) << '\t' << scaf.seq().get(i) << '\t' << snp_info.total() << '\t'<< snp_info.A() << '\t'<< snp_info.C()
			<< '\t' << snp_info.G() << '\t'<< snp_info.T() << '\t' << snp_info.N() << '\t' << snp_info.insertion() << '\t' << snp_info.deletion();
		if((snp_info.total()-snp_info.N()) > coverage_threshold) {
			ossnps << '\t';
			if(abundant != scaf.seq().get(i)) {
				ossnps << scaf.seq().get(i) << " -> " << abundant;
			}
			ossnps << '\t';
			if((1.0-double(n_abundant)/(snp_info.total()-snp_info.N())) >= snp_threshold) {
				ossnps << "SNP";	//\t" << n_abundant << "\t" << snp_info.total() << "\t" << snp_info.N() << endl;
			}
		}
		else 
			ossnps << "\t\t";

		if(snp_info.total() > coverage_threshold) {
			ossnps << '\t';
			if(double(snp_info.insertion())/snp_info.total() >= snp_threshold) { 
				ossnps << "INSERTION";
			}
		}
		else 
			ossnps << "\t";

		if(snp_info.total()+snp_info.deletion() > coverage_threshold) {
			ossnps << '\t';
			if(double(snp_info.deletion())/(snp_info.total()+snp_info.deletion()) >= snp_threshold) { 
				ossnps << "DELETION";
			}
		}
		ossnps << endl;
	}
}

/****************************************************************************************************************/
void usage(const char* prog_name)
{
	const char* p = strrchr(prog_name, '/');
	p = p? (p+1) : prog_name;
	cerr << endl;
	cerr << p << " version " << VERSION << ", " << DATE << endl << endl;
	cerr << "Computes base calls, insertions and deletions for all positions on every sequence in the input files." << endl << endl;
	cerr << "Usage: " << prog_name << " -assembly <sequence-file> -out <output-file> [--phred33-quals|--phred64-quals] [--silent]" << endl;
	cerr << "           [-snp <snp-threshold>] [-sam  <read-mapping-file1> <read-mapping-file2> ... <read-mapping-filen>]" << endl;
	cerr << "Where" << endl;
	cerr << "  -assembly         :  path to assembly file" << endl;
	cerr << "  -out              :  path for the output file" << endl;
	cerr << "  -sam              :  mapping files are in SAM format (only option in current version)" << endl;
	cerr << "  -snp              :  threshold for determining SNPs (default: 0.2)" << endl;
	cerr << "  --phred33-quals   :  Quality scores of reads are ASCII characters equal to the Phred quality plus 33 (same as bowtie)" << endl; 
	cerr << "  --phred64-quals   :  Quality scores of reads are ASCII characters equal to the Phred quality plus 64 (default, same as bowtie)" << endl; 
	cerr << "  --silent          :  progrsm will do its work without status messages to stderr (default: verbose)" << endl << endl;
	cerr << "OUTPUT" << endl;
	cerr << "  output consists of two type of lines:" << endl;
	cerr << "  - Header lines    :  mark the beginning of a new sequence. Format is" << endl;
	cerr << "                    :  ><reference-seq-name> /size=<size> /coverage=<coverage> /gc=<%G+C> /Ns=<# of Ns>" << endl;
	cerr << "  - Base call lines :  present information for a certain position of the sequence whose header was last seen." << endl;    
	cerr << "                    :  Format is" << endl;
	cerr << "                    :  <pos> <reference-char> <coverage> <# A's> <# C's> <# G's> <# T's> <# N's> <# Insertions> <# Deletions> <additional-columns>" << endl;
	cerr << "                    :  where <additional-columns> are either empty or contains the following:" << endl;
	cerr << "                    :  + X -> Y, if X is the reference char but Y is the consensus" << endl;
	cerr << "                    :  + SNP if the number of calls of the non-consensus char is more than <SNP-threshold> of the A, C, G and T (but not N) calls" << endl;
	cerr << "                    :  + INSERTION if the number of insertion calls is more than <SNP-threshold> of A, C, G, T AND N calls for the position" << endl;
	cerr << "                    :  + DELETION if the number of deletion calls if more than <SNP-threshold> of A, C, G, T, N AND deletion call for the position" << endl;
	cerr << "                    :  All fields, including additional columns, are separated by tabs." << endl;
	cerr << "                    :  Insertions refer to insertions coming right after the position in which they are recognized. Insertions can be of 1 or more" << endl;
	cerr << "                    :  bases, in this version only the presence of an insertion is reported." << endl;
	cerr << endl << "Please report bugs to itai.sharon@gmail.com" << endl << endl;
}

/****************************************************************************************************************/
bool read_arguments(int argc, const char* argv[], set<string>& sam_files, string& assembly_file, string& output_file, bool& verbose, double& snp_threshold)
{
	if((argc == 1) || ((argc ==  2) && (!strcmp(argv[0], "--help"))))
	{
		usage(argv[0]);
		return false;
	}

	assembly_file = "";

	enum {undef, sam, assembly, out_file, snp_thresh}	curr_flag = undef;
	for(int i=1; i<argc; i++) {
		if(argv[i][0] == '-') {
			if(!strcmp(argv[i], "-sam")) {
				curr_flag = sam;
			}
			else if(!strcmp(argv[i], "-assembly")) {
				curr_flag = assembly;
			}
			else if(!strcmp(argv[i], "-out")) {
				curr_flag = out_file;
			}
			else if(!strcmp(argv[i], "-snp")) {
				curr_flag = snp_thresh;
			}
			else if(!strcmp(argv[i], "--phred33-quals")) {
				Bio::ReadSequence::Qscore_scheme quality_offset = Bio::ReadSequence::phred33_quals;
			}
			else if(!strcmp(argv[i], "--phred64-quals")) {
				Bio::ReadSequence::Qscore_scheme quality_offset = Bio::ReadSequence::phred64_quals;
			}
			else if(!strcmp(argv[i], "--silent")) {
				verbose = false;
			}
			else {
				throw invalid_argument(argv[i]);
			}
		}
		else {
			string arg = argv[i];
			if(curr_flag == sam) {
				if(!file_exists(arg.c_str()))
					throw invalid_argument(string("File does not exist: ") + arg);
				sam_files.insert(arg);
			}
			else if(curr_flag == assembly) {
				if(assembly_file.size() > 0)
					throw invalid_argument("assembly file is defined twice");
				if(!file_exists(arg.c_str()))
					throw invalid_argument(string("File does not exist: ") + arg);
				assembly_file = arg;
			}
			else if(curr_flag == out_file) {
				if(output_file.size() != 0)
					throw invalid_argument(string("output file specified twice: ") + arg + string(" and ") + output_file);
				output_file = arg;
			}
			else if(curr_flag == snp_thresh) {
				snp_threshold = atof(arg.c_str());
				if((snp_threshold <= 0) || (snp_threshold >= 1))
					throw invalid_argument(string("SNP threshold (-snp) must be higher than 0 and lower than 1 but is set to ") + arg);
				curr_flag = undef;
			}
			else if(curr_flag == undef) {
				throw invalid_argument(string("Did not specify any flag before ") + arg);
			}
		}
	}
	if(assembly_file.size() == 0) {
		throw invalid_argument("assembly file was not defined");
	}
	if(sam_files.size() == 0) {
		throw invalid_argument("No mapping files were supplied");
	}
	return true;
}
