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
#include "nexus_shell.hpp"


static long gStrictLevel = 2;
static long gTreeNum = -1;

RNG InferenceKernel::rng;



static void updateKernelTaxa(InferenceKernelBundle & bundle, const NxsTaxaBlockAPI * taxaBPtr);
static void updateKernelCharacters(InferenceKernelBundle & bundle, const NxsCharactersBlock * charactersBPtr);
static void updateKernelTrees(InferenceKernelBundle & bundle, const NxsTreesBlock * treesBPtr);

#if ! defined(NDEBUG)
	void nexus_shell_assert(const char *e, const char * func, const char * fileN, long line)
	{
	std::cerr << "Assertion failed.\nAssertion = " << e << "\nFunction = " << func;
	std::cerr << "\nFile = " << fileN << "\nLine = " << line << '\n';
	exit(1);
	}
#endif

void updateKernelTaxa(InferenceKernelBundle & bundle, const NxsTaxaBlockAPI * taxaBPtr)
{
	InferenceKernel * kernel = &(bundle.kernel);
	bundle.errStream << bundle.outputPrefix << ": adding taxa" << endl;
	const std::vector<std::string> labels = taxaBPtr->GetAllLabels();
	const std::vector<std::string> * oldSamples = NULL;
	unsigned numOld = 0;
	if (kernel) {
		oldSamples = &(kernel->getTaxLabels());
		numOld = oldSamples->size();
		if (numOld && (labels.size() != numOld)) {
			throw NxsException("Differing number of taxa labels encountered");
		}
	}
	if (kernel)
		kernel->addTaxLabels(labels);
}

void updateKernelCharacters(InferenceKernelBundle & bundle, const NxsCharactersBlock * charactersBPtr)
{
	InferenceKernel * kernel = &(bundle.kernel);
	if (kernel && charactersBPtr) {
		NxsCharactersBlock::DataTypesEnum dt = charactersBPtr->GetDataType();
		if (dt == NxsCharactersBlock::nucleotide || dt == NxsCharactersBlock::dna || dt == NxsCharactersBlock::rna) {
			bundle.errStream << bundle.outputPrefix << ": adding character data" << endl;
			kernel->replaceData(charactersBPtr);
		}
		else {
			bundle.errStream << bundle.outputPrefix << ": non-nucleotide character data" << endl;
		}
	}
}	

void updateKernelTrees(InferenceKernelBundle & bundle, const NxsTreesBlock * treesBPtr)
{
	InferenceKernel * kernel = &(bundle.kernel);
	if (!kernel || ! treesBPtr)
		return;
	unsigned n = treesBPtr->GetNumTrees();
	bundle.errStream << bundle.outputPrefix << ": adding " << n << " tree(s)..." << endl;
	treesBPtr->ProcessAllTrees();
	for (unsigned i = 0; i < n; ++i)
		kernel->addTree(treesBPtr->GetFullTreeDescription(i));
	bundle.errStream << bundle.outputPrefix << ": " << kernel->getNumTreeDescriptions() << " tree description(s) now in memory" << endl;
}



void readFilepath(const char * filename, InferenceKernelBundle & bundle)
{
	NEXUS_SHELL_ASSERT(filename);
	BlockReaderList blocks;
	try {
		ExceptionRaisingNxsReader nexusReader(NxsReader::WARNINGS_TO_STDERR);
		if (gStrictLevel != 2)
			nexusReader.SetWarningToErrorThreshold((int)NxsReader::FATAL_WARNING + 1 - (int) gStrictLevel);
		ifstream inf(filename, ios::binary);
		if (!inf.good()) {
			NxsString err;
			err << "Could not parse the file \"" << filename <<"\"";
			nexusReader.NexusError(err, 0, -1, -1);
		}
		NxsToken token(inf);	
		NxsCloneBlockFactory factory;

		nexusReader.AddFactory(&factory);
		NxsCharactersBlock charsB(NULL, NULL);
		charsB.SetCreateImpliedBlock(true);
		charsB.SetImplementsLinkAPI(true);
		charsB.SetSupportMixedDatatype(true);
		
		NxsDataBlock dataB(NULL, NULL);
		dataB.SetCreateImpliedBlock(true);
		dataB.SetImplementsLinkAPI(true);
		dataB.SetSupportMixedDatatype(true);
		
		NxsTaxaBlock taxaB;
		taxaB.SetImplementsLinkAPI(false);
		
		NxsTreesBlock treesB(NULL);
		treesB.SetCreateImpliedBlock(true);
		treesB.SetImplementsLinkAPI(true);
		treesB.SetProcessAllTreesDuringParse(gTreeNum < 0);
		if (gStrictLevel < 2)
			treesB.SetAllowImplicitNames(true);
		treesB.SetWriteFromNodeEdgeDataStructure(true);
		
		NxsStoreTokensBlockReader storerB(bundle.blockName, true);

		factory.AddPrototype(&charsB, "CHARACTERS");
		factory.AddPrototype(&dataB, "DATA");
		factory.AddPrototype(&taxaB);
		factory.AddPrototype(&treesB);
		factory.AddPrototype(&storerB);
		try {
			bundle.errStream << bundle.outputPrefix << ": Parsing " << filename << endl;
			nexusReader.Execute(token);
			bundle.errStream << bundle.outputPrefix << ": Interpreting " << filename << endl;
		}
		catch(...) {
			nexusReader.RemoveFactory(&factory);
			throw;
		}
		nexusReader.RemoveFactory(&factory);
		blocks = nexusReader.GetUsedBlocksInOrder();
		for (BlockReaderList::const_iterator bIt = blocks.begin(); bIt != blocks.end(); ++bIt) {
			NxsBlock * b = *bIt;
			if (b) {
				if (b->GetID() == "TAXA") {
					NxsTaxaBlockAPI * taxaBPtr = (NxsTaxaBlockAPI *) b;
					updateKernelTaxa(bundle, taxaBPtr);
				}
				else if (b->GetID() == "CHARACTERS" || b->GetID() == "DATA" ) {
					NxsCharactersBlock * charactersBPtr = (NxsCharactersBlock *) b;
					updateKernelCharacters(bundle, charactersBPtr);
				}
				else if (b->GetID() == "TREES") {
					NxsTreesBlock * treesBPtr = (NxsTreesBlock *) b;
					updateKernelTrees(bundle, treesBPtr);
				}
				else if (b->GetID() == bundle.blockName) {
					NxsStoreTokensBlockReader * privBlockPtr = (NxsStoreTokensBlockReader*) b;
					const std::list<ProcessedNxsCommand> & cmds = privBlockPtr->GetCommands();
					for (std::list<ProcessedNxsCommand>::const_iterator cIt = cmds.begin(); cIt != cmds.end(); ++cIt) {
						NxsSimpleCommandStrings s = ProcessedNxsToken::ParseSimpleCmd(*cIt, true);
						(*bundle.commandReader)(s, bundle);
					}
				}
			}
		}
		for (BlockReaderList::const_iterator bIt = blocks.begin(); bIt != blocks.end(); ++bIt) {
			NxsBlock * b = *bIt;
			delete b;
		}
		bundle.errStream << bundle.outputPrefix << ": Finished with " << filename << endl;
	}
	catch (const NxsException &x) {
		cerr << "Error:\n " << x.msg << endl;
		if (x.line > 0)
			cerr << "at line " << x.line << ", column (approximately) " << x.col << " (and file position "<< x.pos << ")" << endl;
		exit(2);
	}
}

void readFilepathOrExit(const char *filename, InferenceKernelBundle & bundle) {
	bundle.errStream << bundle.outputPrefix << ": Reading " << filename << "	 " << std::endl;
	try {
		readFilepath(filename, bundle);
	}
	catch (const NexusShellException &x) {
			bundle.errStream << bundle.outputPrefix << ": Reading of " << filename << " failed (with a NexusShellException).\n Exception message: " << x.what() << endl;
		exit(1);
	}
	catch (...) {
		bundle.errStream << bundle.outputPrefix << ": Reading of " << filename << " failed (with an unknown exception)" << endl;
		exit(1);
	}
}	


bool doPrivateInferenceCommand(const NxsSimpleCommandStrings &s, InferenceKernelBundle & bundle)
{
	InferenceKernel * kernel = (InferenceKernel *)(&(bundle.kernel));
	std::vector<std::string> opts;
	std::vector<double> dv;
	std::vector<double> secondDV;
	double doubleVar;
	/******** Execute *********************************************************/
	if (s.cmdName == "execute") {
		opts = s.GetOptValue("file").second;
		if (opts.empty())
			throw NexusShellException("Expecting a file option in the Execute command");
		std::string fn = *(opts.rbegin());
		readFilepath(fn.c_str(), bundle);
		return true;
	}
	/******** MCMC ************************************************************/
	if (s.cmdName == "mcmc") {
		std::cerr << "Starting MCMC.\n";
		kernel->runMCMC(-1);
		std::cerr << "Finished with MCMC.\n";
		return true;
	}
	/******** Move ************************************************************/
	if (s.cmdName == "move") {
		opts = s.GetOptValue("weight").second;
		double wt = 1.0;
		if (opts.size() == 1) {
			if ((!NxsString::to_double(opts.rbegin()->c_str(), &wt)|| (wt < 0.0)))
				throw NexusShellException("move weight must be a positive number.");
		}
		opts = s.GetOptValue("modelnumber").second;
		long submdodelNumber = 1;
		if (opts.size() == 1) {
			if ((!NxsString::to_long(opts.rbegin()->c_str(), &submdodelNumber) || (submdodelNumber > INT_MAX)))
				throw NexusShellException("modelindex must be a positive number.");
		}
		const unsigned submodelIndex = ((unsigned)submdodelNumber) - 1;
		opts = s.GetOptValue("freqwindow").second;
		if (opts.size() > 0) {
			const unsigned nStates = 4;
			if (opts.size() != nStates)
				throw NexusShellException("Expecting 4 base frequency move windows.");
			dv.clear();
			for (unsigned i = 0; i < nStates; ++i) {
				if ((!NxsString::to_double(opts[i].c_str(), &doubleVar)|| (doubleVar < 0.0)))
					throw NexusShellException("each base frequency move window must be a positive number.");
				dv.push_back(doubleVar);
			}
			secondDV.assign(nStates, 1.0);
			opts = s.GetOptValue("freqpriormean").second;
			if (opts.size() > 0) {
				if (opts.size() != nStates)
					throw NexusShellException("Expecting 4 base frequency prior means in the move command.");
				secondDV.clear();
				for (unsigned i = 0; i < nStates; ++i) {
					if ((!NxsString::to_double(opts[i].c_str(), &doubleVar)|| (doubleVar < 0.0)))
						throw NexusShellException("each base prior mean must be a positive number.");
					secondDV.push_back(doubleVar);
				}
			}
			kernel->addStateFreqMove(-1, wt, submodelIndex, dv, secondDV);
		}

		opts = s.GetOptValue("rmatwindow").second;
		if (opts.size() > 0) {
			const unsigned nRates = 6;
			if (opts.size() != nRates)
				throw NexusShellException("Expecting 6 rmat rate windows.");
			dv.clear();
			for (unsigned i = 0; i < nRates; ++i) {
				if ((!NxsString::to_double(opts[i].c_str(), &doubleVar)|| (doubleVar < 0.0)))
					throw NexusShellException("each rmat rate move window must be a positive number.");
				dv.push_back(doubleVar);
			}
			secondDV.assign(nRates, 1.0);
			opts = s.GetOptValue("rmatpriormean").second;
			if (opts.size() > 0) {
				if (opts.size() != nRates)
					throw NexusShellException("Expecting 6 rmat prior means in the move command.");
				secondDV.clear();
				for (unsigned i = 0; i < nRates; ++i) {
					if ((!NxsString::to_double(opts[i].c_str(), &doubleVar)|| (doubleVar < 0.0)))
						throw NexusShellException("each base prior mean must be a positive number.");
					secondDV.push_back(doubleVar);
				}
			}
			kernel->addRateMatMove(-1, wt, submodelIndex, dv, secondDV);
		}
		
		opts = s.GetOptValue("edgelenwindow").second;
		if (opts.size() > 0) {
			if (opts.size() != 1)
				throw NexusShellException("Expecting 1 window size for the edgelen move");
			if ((!NxsString::to_double(opts[0].c_str(), &doubleVar)|| (doubleVar < 0.0)))
				throw NexusShellException("the edgelen move window must be a positive number.");
			opts = s.GetOptValue("edgelenpriormean").second;
			double priorMeanV = 0.1;
			if (opts.size() > 0) {
				if (opts.size() != 1)
					throw NexusShellException("Expecting 1 edgelen prior means in the move command.");
				if ((!NxsString::to_double(opts[0].c_str(), &priorMeanV)|| (priorMeanV < 0.0)))
					throw NexusShellException("each base prior mean must be a positive number.");
			}
			kernel->addEdgeLenMove(-1, wt, submodelIndex, doubleVar, priorMeanV);
		}
		
		NxsSimpleCommandStrings::MatFromFile matFF = s.GetMatOptValue("partitionedfreqwindow");
		const NxsSimpleCommandStrings::MatString & mat = matFF.second;
		if (!mat.empty()) {
			std::vector<double> row;
			std::list<std::vector<double> > pfwl;
			for (NxsSimpleCommandStrings::MatString::const_iterator mIt = mat.begin(); mIt != mat.end(); ++mIt) {
				const std::vector<std::string> & mrow = *mIt;
				row.clear();
				for (std::vector<std::string>::const_iterator c = mrow.begin(); c != mrow.end(); ++c) {
					double doubleVar;
					const std::string & o = *c;
					if ((!NxsString::to_double(o.c_str(), &doubleVar)|| (doubleVar < 0.0))) {
						throw NexusShellException("each partitioned move window must be a non-negative number.");
					}
					row.push_back(doubleVar);
				}
				pfwl.push_back(row);
			}
			std::list<std::vector<double> >::const_iterator wIt = pfwl.begin();

			std::list<std::vector<double> > pfml;
			matFF = s.GetMatOptValue("partitionedfreqmean");
			const NxsSimpleCommandStrings::MatString & mat = matFF.second;
			NxsString err;
			if (!mat.empty()) {
				if (mat.size() != pfwl.size()) {
					err << "Expecting the same number of prior means as window sizes in the PartitionedFreq move. Found " << mat.size() <<" vectors of prior means, and " << pfwl.size() << " vectors of window sizes.";
					throw NexusShellException(err);
				}
				unsigned subset_ind = 0;

				for (NxsSimpleCommandStrings::MatString::const_iterator mIt = mat.begin(); mIt != mat.end(); ++mIt, ++wIt, ++subset_ind) {
					const std::vector<std::string> & mrow = *mIt;
					if (mrow.size() != wIt->size()) {
						err << "Expecting the same number of prior means as window sizes in the PartitionedFreq move. Found " << mrow.size() <<" prior means, and " << wIt->size() << " window sizes for subset " << subset_ind +1 << '.';
						throw NexusShellException(err);
					}
					row.clear();
					for (std::vector<std::string>::const_iterator c = mrow.begin(); c != mrow.end(); ++c) {
						double doubleVar;
						const std::string & o = *c;
						if ((!NxsString::to_double(o.c_str(), &doubleVar)|| (doubleVar < 0.0))) {
							throw NexusShellException("each partitioned move mean must be a non-negative number.");
						}
						row.push_back(doubleVar);
					}
					pfml.push_back(row);
				}
			}
			else {
				for (; wIt != pfwl.end() ;++wIt) {
					row.assign(wIt->size(), 1.0);
					pfml.push_back(row);
				}
			}
			kernel->addPartitionedFreqMove(-1, wt, pfwl, pfml);
		}
		return true;
	}
	/******** SaveTree *********************************************************/
	if (s.cmdName == "savetree") {
		opts = s.GetOptValue("file").second;
		if (opts.empty())
			opts.push_back(std::string("nexus_shell_tree.nex"));
		for (unsigned i = 0; i < opts.size(); ++i) {
			std::string currFilename = opts[i];
			
			std::ofstream outStream(opts[i].c_str());
			std::cerr << "Saving the last tree in memory to " << currFilename  << " ...\n";
			outStream << "#NEXUS\n";
			kernel->writeTaxa(outStream);
			kernel->writeTree(outStream, -1, -1);
			outStream.close();
			std::cerr << "Finished saving tree.\n";
		}
		return true;
	}

	/******** Set *********************************************************/
	if (s.cmdName == "set") {
		long longVar;
		opts = s.GetOptValue("seed").second;
		if (!opts.empty()) {
			if ((!NxsString::to_long(opts.rbegin()->c_str(), &longVar)|| (longVar < 1)))
				throw NexusShellException("seed must be a positive number.");
			InferenceKernel::rng.setSeed((unsigned) longVar);
		}

		opts = s.GetOptValue("geneticcode").second;
		if (!opts.empty()) {
			std::string codename = *(opts.rbegin());
			NxsGeneticCodesEnum codeE = geneticCodeNameToEnum(codename);
			kernel->setGeneticCode(codeE);
		}

		opts = s.GetOptValue("sampleinterval").second;
		if (!opts.empty()) {
			if ((!NxsString::to_long(opts.rbegin()->c_str(), &longVar)|| (longVar < 1)))
				throw NexusShellException("sampleInterval must be a positive number.");
			kernel->setSampleInterval((unsigned) longVar);
		}

		opts = s.GetOptValue("reportinterval").second;
		if (!opts.empty()) {
			if ((!NxsString::to_long(opts.rbegin()->c_str(), &longVar)|| (longVar < 1)))
				throw NexusShellException("reportInterval must be a positive number.");
			kernel->setReportInterval((unsigned) longVar);
		}

		opts = s.GetOptValue("checkpointinterval").second;
		if (!opts.empty()) {
			if ((!NxsString::to_long(opts.rbegin()->c_str(), &longVar)|| (longVar < 1)))
				throw NexusShellException("checkpointInterval must be a positive number.");
			kernel->setCheckpointInterval((unsigned) longVar);
		}

		opts = s.GetOptValue("currmcmciteration").second;
		if (!opts.empty()) {
			if ((!NxsString::to_long(opts.rbegin()->c_str(), &longVar) || (longVar < 1) || (longVar >= LONG_MAX)))
				throw NexusShellException("currmcmciteration must be a positive number (but not larger than INT_MAX)");
			kernel->setCurrMCMCIterations((unsigned long) longVar);
		}

		opts = s.GetOptValue("mcmciterations").second;
		if (!opts.empty()) {
			if ((!NxsString::to_long(opts.rbegin()->c_str(), &longVar) || (longVar < 1) || (longVar >= LONG_MAX)))
				throw NexusShellException("mcmciterations must be a positive number (but not larger than INT_MAX)");
			kernel->setMCMCIterations((unsigned long) longVar);
		}

		opts = s.GetOptValue("checkpointfile").second;
		if (!opts.empty()) {
			kernel->setCheckpointFilename(opts[0]);
		}

		opts = s.GetOptValue("sampletreefile").second;
		if (!opts.empty()) {
			kernel->setSampleTreeFilename(opts[0]);
		}

		opts = s.GetOptValue("sampleparamfile").second;
		if (!opts.empty()) {
			kernel->setSampleParamFilename(opts[0]);
		}

		opts = s.GetOptValue("reportfile").second;
		if (!opts.empty()) {
			kernel->setReportFilename(opts[0]);
		}
		
		opts = s.GetOptValue("rmat").second;
		if (!opts.empty()) {
			if (opts.size() != 6)
				throw NexusShellException("expecting 6 rmat rates");
			dv.clear();
			for (unsigned i = 0; i < 6; ++i) {
				if ((!NxsString::to_double(opts[i].c_str(), &doubleVar)|| (doubleVar < 0.0)))
					throw NexusShellException("each rate in the rmat must be a positive number.");
				dv.push_back(doubleVar);
			}
			kernel->setRMat(dv);
		}
		opts = s.GetOptValue("freq").second;
		if (!opts.empty()) {
			if (opts.size() != 4)
				throw NexusShellException("expecting 4 freq rates");
			dv.clear();
			for (unsigned i = 0; i < 4; ++i) {
				if ((!NxsString::to_double(opts[i].c_str(), &doubleVar)|| (doubleVar < 0.0)))
					throw NexusShellException("each freq must be a positive number.");
				dv.push_back(doubleVar);
			}
			kernel->setStateFreq(dv);
		}
		
		opts = s.GetOptValue("shape").second;
		if (!opts.empty()) {
			if ((!NxsString::to_double(opts.rbegin()->c_str(), &doubleVar)|| (doubleVar < 0.0)))
				throw NexusShellException("shape must be a positive number.");
			kernel->setGamma((unsigned) doubleVar);
		}

		opts = s.GetOptValue("pinv").second;
		if (!opts.empty()) {
			if ((!NxsString::to_double(opts.rbegin()->c_str(), &doubleVar)|| (doubleVar < 0.0)))
				throw NexusShellException("pinv must be a positive number.");
			kernel->setPinv((unsigned) doubleVar);
		}
		return true;
	}
	/******** SeqLScore ********************************************************/
	if (s.cmdName == "seqlscore") {
		std::cerr << "Calculating the likelihood of the tree for the sequence data.\n";
		const double lnL = kernel->seqLScore(-1);
		std::cerr << "lnL = " << lnL << '\n';
		std::cerr << "seqlscore.\n";
		return true;
	}
	return false;
}

//==============================================================================
//		rand.cpp taken from the Lot class in phycas
//==============================================================================
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
|                                                                             |
|  This program is free software; you can redistribute it and/or modify       |
|  it under the terms of the GNU General Public License as published by       |
|  the Free Software Foundation; either version 2 of the License, or          |
|  (at your option) any later version.                                        |
|                                                                             |
|  This program is distributed in the hope that it will be useful,            |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of             |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              |
|  GNU General Public License for more details.                               |
|                                                                             |
|  You should have received a copy of the GNU General Public License along    |
|  with this program; if not, write to the Free Software Foundation, Inc.,    |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

const unsigned MASKSIGNBIT = 0x80000000;

/*----------------------------------------------------------------------------------------------------------------------
|	Member function returning uniform deviate. Provided by J. Monahan, Statistics Dept., North Carolina State 
|	University. Originally from Schrage, ACM Trans. Math. Software 5:132-138 (1979). Translated to C by Paul O. Lewis, 
|	Dec. 10, 1992.
*/
#if defined(LOG_RNG_UNIFORM_CALLS)
	double RNG::uniform(const char * file, const int line)
#else
	double RNG::uniform()
#endif
	{
	const unsigned a = 16807U;
	const unsigned b15 = 32768U;
	const unsigned b16 = 65536U;
	const unsigned p = 2147483647U;
	const unsigned xhi = curr_seed / b16;
	const unsigned xalo = (curr_seed - xhi * b16) * a;
	const unsigned leftlo = xalo / b16;
	const unsigned fhi = xhi * a + leftlo;
	const unsigned k = fhi / b15;
	curr_seed = (((xalo - leftlo * b16) - p) + (fhi - k * b15) * b16) + k;
	if (curr_seed & MASKSIGNBIT) 
		curr_seed = unsigned (((int)curr_seed) + p);
#	if defined(LOG_RNG_UNIFORM_CALLS)
#		pragma TODO("undef LOG_RNG_UNIFORM_CALLS in std_force_include.hpp for release version")
		double retval = curr_seed * 4.6566128575e-10;
		std::ofstream tmpf("uniform.txt", std::ios::out | std::ios::app);
		tmpf << num_seeds_generated << " -> " << retval << " (" << file << ":" << line << ")" << std::endl;
		tmpf.close();
		num_seeds_generated++;
		return retval;
#	else
		return curr_seed * 4.6566128575e-10;
#	endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns an unsigned integer in [0-`max') in which all values are equiprobable.
*/
unsigned RNG::sampleUInt(unsigned max)
	{
	NEXUS_SHELL_ASSERT(max > 0);

	unsigned samples_uint = max;
	while(samples_uint == max) 
		{
		double r = uniform(FILE_AND_LINE);
		samples_uint = (unsigned)((double)max*r);
		}

	return samples_uint;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a multinomial deviate in [0,n) given the bin probabilities in probs. The sum of the first n elements of 
|   probs is assumed to equal 1.0.
*/
unsigned RNG::multinomialDraw(const double * probs, unsigned n, double totalProb)
	{
    NEXUS_SHELL_ASSERT(probs != NULL);
    NEXUS_SHELL_ASSERT(n > 0);
    double u = totalProb*uniform(FILE_AND_LINE);
    for (unsigned i = 0; i < n; ++i)
        {
        u -= probs[i];
        if (u < 0.0)
            return i;
        }
    return n-1;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of data member `last_seed_setting', which stores the seed used to initialize generator.
*/
void RNG::useClockToSeed()
	{
	time_t timer;
	setSeed((unsigned)time(&timer));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Designed to emulate Python function random.getrandbits(k). Returns an unsigned value in which the low `nbits' bits
|	are randomly chosen. Assumes `nbits' is greater than 0 and less than 32.
|>
|	nbits   decimal             binary
|	-----------------------------------------------------------------------
|	1       0, 1                0 or 1
|	2		0, 1, 2, 3	        00, 01, 10, 11
|	3       0, 1, 2, ..., 7     000, 001, 010, 011, 100, 101, 110, 111
|	-----------------------------------------------------------------------
|>
|	In general, specifying `nbits' to n results in a random choice of unsigned values in the range [0, 1, ..., (2^n)-1].  
*/
unsigned RNG::getRandBits(unsigned nbits)
	{
	NEXUS_SHELL_ASSERT(nbits > 0);
	NEXUS_SHELL_ASSERT(nbits < 32);

	double u = uniform(FILE_AND_LINE);

	if (u == 0.0)
		{
		return 0;
		}
	else
		{
		// e.g. for nbits = 3, u = 0.99
		//   term1 = log(8)
		//   term2 = log(0.99)
		//   term3 = 0.99*8 = 7.92
		//   return floor(7.92) = 7
		double term1 = log(2.0)*(double)nbits;
		double term2 = log(u);
		double term3 = exp(term1 + term2);
		return (unsigned)floor(term3);
		}
	}


