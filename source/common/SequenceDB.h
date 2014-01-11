/*
 * SequenceDB.h
 *
 *  Created on: 06/Feb/2013
 *      Author: Itai Sharon, itai.sharon@gmail.com
 */

#ifndef SEQUENCEDB_H_
#define SEQUENCEDB_H_

#include <map>
#include "common.h"
#include "SeqIORead_fasta.h"
#include "String.h"
#include "bio_exceptions.h"


using namespace std;

namespace Bio {

/*****************************************************************************************************************
 * SequenceDB
 * A simple db for sequences that provides basic functionality. Will be developed further as required
 ****************************************************************************************************************/
template <class S>
class SequenceDB {
public:
	class iterator;
	class const_iterator;
public:
	SequenceDB(SeqIORead_fasta<S>& reader);
	SequenceDB()				{}
	virtual 				~SequenceDB();
	void					clear();
	const SequenceDB&			operator += (S& seq)			{m_sequences.insert(pair<string, S*>(seq.display_id(), &seq));}
	void					erase(const string& id);
	void					erase(iterator it);
	size_t					size() const				{return m_sequences.size();}
	string					db_file() const				{return m_db_file;}	
	iterator				find(string display_id);
	iterator				begin();
	iterator				end();
	const_iterator				find(string display_id) const;
	const_iterator				begin() const;
	const_iterator				end() const;
protected:
	map<string, S*>				m_sequences;
	string 					m_db_file;
};

/*****************************************************************************************************************
 * SequenceDB::iterator
 * This is not exactly STL's const_iterator but it contains the two important operations - operators * and ++.
 ****************************************************************************************************************/
template <class S>
class SequenceDB<S>::iterator {
protected:
	map<string, S*>*			sequences;
	typename map<string, S*>::iterator	it;
public:
	S&					operator *()				{return *(it->second);}
	S* 					operator->()				{return it->second;}
	// Not a mistake - returned value for the following two is const iterator
	const iterator				operator ++(int dummy)			{iterator _it = *this; it++; return _it;}
	const iterator				operator ++()				{it++; return *this;}
	bool operator == (const iterator it2) const					{return((sequences==it2.sequences) && (it==it2.it));}
	bool operator != (const iterator it2) const					{return !(*this == it2);}
protected:
	friend class SequenceDB;
	iterator(map<string, S*>& _sequences, typename map<string, S*>::iterator _it) :
		sequences(&_sequences), it(_it) 
		{}
	iterator(map<string, S*>& _sequences) :
		sequences(&_sequences), it(_sequences.begin()) 
		{}
};

/*****************************************************************************************************************
 * SequenceDB::const_iterator
 * This is not exactly STL's const_iterator but it contains the two important operations - operators * and ++.
 ****************************************************************************************************************/
template <class S>
class SequenceDB<S>::const_iterator {
protected:
	const map<string, S*>*				sequences;
	typename map<string, S*>::const_iterator	it;
public:
	const_iterator(const SequenceDB<S>::iterator& ite) : sequences(ite.sequences), it(ite.it)	{}
	const S&				operator *() const			{return *(it->second);}
	const S* 				operator->() const			{return it->second;}
	// Not a mistake - returned value for the following two is const const_iterator
	const const_iterator 			operator ++(int dummy)			{const_iterator _it = *this; it++; return _it;}
	const const_iterator 			operator ++()				{it++; return *this;}
	bool operator == (const const_iterator it2) const				{return((sequences==it2.sequences) && (it==it2.it));}
	bool operator == (const iterator it2) const					{return((sequences==it2.sequences) && (it==it2.it));}
	bool operator != (const const_iterator it2) const				{return !(*this == it2);}
	bool operator != (const iterator it2) const					{return !(*this == it2);}
protected:
	friend class SequenceDB;
	const_iterator(const map<string, S*>& _sequences, typename map<string, S*>::const_iterator _it) :
		sequences(&_sequences), it(_it) 
		{}
	const_iterator(const map<string, S*>& _sequences) :
		sequences(&_sequences), it(_sequences.begin()) 
		{}
};

/****************************************************************************************************************/
template <class S>
SequenceDB<S>::SequenceDB(SeqIORead_fasta<S>& reader) : m_db_file(reader.file_name())
{
	reader.restart();
	S* p_seq = NULL;

	while((p_seq = reader.next_seq()) != NULL) 
	{
		if(m_sequences.find(p_seq->display_id()) != m_sequences.end())
			throw Bio::Bad_file(reader.file_name(), string("Same sequence id appears more than obce: ") + p_seq->display_id(), __PRETTY_FUNCTION__);
		m_sequences.insert(pair<string, S*>(p_seq->display_id(), p_seq));
	}
}

/****************************************************************************************************************/
template <class S>
SequenceDB<S>::~SequenceDB()
{
	this->clear();
}

/****************************************************************************************************************/
template <class S>
void SequenceDB<S>::clear()
{
	for(typename map<string, S*>::iterator it=m_sequences.begin(); it != m_sequences.end(); it++) {
		delete(it->second);
		it->second = NULL;
	}
	m_sequences.clear();
}

/****************************************************************************************************************/
template <class S>
void SequenceDB<S>::erase(const string& id)
{
	iterator it = this->find(id);
	this->erase(it);
}

/****************************************************************************************************************/
template <class S>
void SequenceDB<S>::erase(iterator it)
{
	if(it.sequences != &m_sequences)
		return;

	if(it != this->end()) {
		m_sequences.erase(it.it);
	}
}

/****************************************************************************************************************/
template <class S>
typename SequenceDB<S>::const_iterator SequenceDB<S>::find(string display_id) const
{
	return SequenceDB<S>::const_iterator(m_sequences, m_sequences.find(display_id));
}

/****************************************************************************************************************/
template <class S>
typename SequenceDB<S>::iterator SequenceDB<S>::begin()
{
	return iterator(m_sequences);
}

/****************************************************************************************************************/
template <class S>
typename SequenceDB<S>::iterator SequenceDB<S>::end()
{
	return iterator(m_sequences, m_sequences.end());
}

/****************************************************************************************************************/
template <class S>
typename SequenceDB<S>::iterator SequenceDB<S>::find(string display_id)
{
	return SequenceDB<S>::iterator(m_sequences, m_sequences.find(display_id));
}

/****************************************************************************************************************/
template <class S>
typename SequenceDB<S>::const_iterator SequenceDB<S>::begin() const
{
	return const_iterator(m_sequences);
}

/****************************************************************************************************************/
template <class S>
typename SequenceDB<S>::const_iterator SequenceDB<S>::end() const
{
	return const_iterator(m_sequences, m_sequences.end());
}

}

#endif /* SEQUENCEDB_H_ */
