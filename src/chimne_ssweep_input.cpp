//	Copyright (C) 2008 Mark T. Holder
//
//	chimne_ssweep is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	chimne_ssweep is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with NCL; if not, write to the Free Software Foundation, Inc., 
//	59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//

/*******************************************************************************
 *	This file contains code for an executable that takes path to a NEXUS file
 *		as a command line argument and writes a "normalized" version of the 
 *		blocks.	 This is useful for testing. 
 *
 *	Ideally 2 equivalent files will produce the same normalized output. This 
 *		version of tthe program is less ambitious. The goal is to be able to run 
 *		(for any valid NEXUS in.nex file):
 *			$ normalizer in.nex > outOrig.nex
 *			$ normalizer outOrig.nex > outSecond.nex
 *			$ diff outOrig.nex outSecond.nex
 *		and find no differences.
 */
#include <iostream>
#include "ncl/ncl.h"
#include "ncl/nxsblock.h"
#include "ncl/nxspublicblocks.h"
#include "ncl/nxscharactersblock.h"
#include "chimne_ssweep_input.hpp"
#include "chimne_ssweep_kernel.hpp"

using namespace chim;

bool parseTaxonLabelIntoSample(const std::string & label, Sample * sample);

void MCMCKernel::addTaxLabels(const std::vector<std::string> & labels)
{
	const std::vector<std::string> & oldLabels = this->getTaxLabels();
	unsigned numOld = oldLabels.size();
	const unsigned nLabels = labels.size();
	std::vector<Sample> newSamples(nLabels);
	for (unsigned i = 0; i < nLabels; ++i) {
		Sample & sample = newSamples[i];
		sample.setIndex(i);
		sample.setLabel(labels[i]);
		if (numOld > i && this->samples[i] != sample) {
			cerr << sample  << '\n' << this->samples[i] << '\n';
			throw NxsException("Differing taxa blocks encountered");
		}
	}
	if (numOld == 0)
		addSamples(newSamples);
}


void chim::ChimneSweepKernel::addTaxLabels(const std::vector<std::string> & labels)
{
	const std::vector<std::string> & oldLabels = this->getTaxLabels();
	unsigned numOld = oldLabels.size();
	const unsigned nLabels = labels.size();
	std::vector<Sample> newSamples(nLabels);
	for (unsigned i = 0; i < nLabels; ++i) {
		Sample & sample = newSamples[i];
		sample.setIndex(i);
		try {
			sample.parseCodedName(labels[i]);
		}
		catch (const std::exception &e) {
			NxsString eMsg;
			eMsg << "Could not parse the name of taxon " << i + 1 << ":\n";
			eMsg << labels[i];
			eMsg << "\ninto a Sample.  Expecting patient.time.cID.Accesion;\n";
			eMsg << e.what();
			throw NxsException(eMsg);
		}
		if (numOld > i && this->samples[i] != sample) {
			cerr << sample  << '\n' << this->samples[i] << '\n';
			throw NxsException("Differing taxa blocks encountered");
		}
	}
	if (numOld == 0)
		addSamples(newSamples);
}

void Sample::parseCodedName(const std::string & label)
{
	const std::string p = this->sampleLabel;
	const std::string a = this->accession;
	const int ini = this->intTime;
	const std::string instrTime = this->strTime;
	const std::string inid = this->id;
	try {
		const size_t last_ind = label.length() - 1;
		size_t beg_pos = 0;
		
		/* read patients name (before first .) */
		size_t end_pos = label.find_first_of('.');
		if (end_pos == string::npos)
			throw ChimneSweepException("no \'.\' found \n");
		this->sampleLabel.assign(label, 0, end_pos);
		
		/* read time (between first and second .) */
		beg_pos = end_pos + 1;
		end_pos = label.find_first_of('.', beg_pos);
		if (end_pos == string::npos)
			throw ChimneSweepException("second \'.\' not found \n");
	
		this->accession.assign(label, beg_pos, end_pos - beg_pos);
		long tmp;
		if (! NxsString::to_long(this->accession.c_str(), &tmp))
			ChimneSweepException("time not a long \n");
		if (tmp > INT_MAX)
			throw ChimneSweepException(" time field exceeds capacity of int");
		this->strTime = this->accession;
		this->intTime = (int) tmp;
	
		/* read cID (between second and third .) */
		beg_pos = end_pos + 1;
		/*
		if (beg_pos > last_ind || label[beg_pos] != 'c')
			throw ChimneSweepException("c not found before ID\n");
		++beg_pos;
		*/
		end_pos = label.find_first_of('.', beg_pos);
		if (end_pos == string::npos)
			throw ChimneSweepException("third \'.\' not found \n");
		this->id.assign(label, beg_pos, end_pos - beg_pos);
		/* read accession (after third .) */
		this->accession.assign(label, end_pos + 1, last_ind);
	}
	catch (...) {
		this->sampleLabel = p;
		this->accession = a;
		this->strTime = instrTime;
		this->intTime = ini;
		this->id = inid;
		throw;
	}
}

void MCMCKernel::replaceData(const NxsCharactersBlock * charactersBPtr)
{
	cerr << "Replacing data and clearing trees constructed in memory\n";
	this->clearTrees();
	this->charMatrix.Initialize(charactersBPtr, true);
	this->charMatrixInitialized = true;
}

bool chim::doChimneSweepNexusCommand(const NxsSimpleCommandStrings &s, InferenceKernelBundle & bundle)
{
	if (doPrivateInferenceCommand(s, bundle))
		return true;
		
	ChimneSweepKernel * kernel = (ChimneSweepKernel *)(&(bundle.kernel));
	std::vector<std::string> opts;
	std::vector<double> dv;
	std::vector<double> secondDV;
	/******** EdgLenFromDepth *************************************************/
	if (s.cmdName == "edgelenfromdepth") {
		std::cerr << "Setting branch lengths of the last tree in memory from the node heights...\n";
		kernel->setEdgeLenFromDepth(-1);
		std::cerr << "Finished setEdgeLenFromDepth routine.\n";
		return true;
	}
	/******** Ultrametricize ********************************************************/
	if (s.cmdName == "ultrametricize") {
		opts = s.GetOptValue("minedgelen").second;
		double minE = 0.0;
		if (opts.size() > 0) {
			if ((!NxsString::to_double(opts.rbegin()->c_str(), &minE)|| (minE < 0.0)))
				throw ChimneSweepException("MinEdgeLen must be a positive number.");
		}
		std::cerr << "Making the last tree in memory ultrametric...\n";
		kernel->ultrametricize(-1, minE);
		std::cerr << "Finished ultrametricize routine.\n";
		return true;
	}
	NxsString e;
	e << "Unknown command \"" << s.cmdName << "\" encountered.\n";
	throw NxsException(e, s.cmdPos);
	return false;
}


	
