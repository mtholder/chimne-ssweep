#if ! defined(NEXUS_SHELL_HPP)
#define NEXUS_SHELL_HPP

#include <string>
#include <vector>
#include <iostream>

#include "ncl/nxscharactersblock.h"
#include "seq_likelihood.h"
class InferenceKernel;
class InferenceKernelBundle;
class NxsSimpleCommandStrings;
class NxsFullTreeDescription;
class SeqModel;


typedef bool (*privateBlockReaderFuncPtr)(const NxsSimpleCommandStrings &, InferenceKernelBundle & );

bool doPrivateInferenceCommand(const NxsSimpleCommandStrings &s, InferenceKernelBundle & bundle);

void readFilepath(const char * filename, InferenceKernelBundle & bundle);
void readFilepathOrExit(const char * filename, InferenceKernelBundle & );

#if defined(NDEBUG)
// no asserts in release version
#	define NEXUS_SHELL_ASSERT(expr)
#else
	void nexus_shell_assert(const char *e, const char * func, const char * fileN, long);
#	define NEXUS_SHELL_ASSERT(expr)  if (!(expr)) nexus_shell_assert((const char *)#expr, (const char *)__FUNCTION__, __FILE__, __LINE__)
#endif


#if defined(LOG_RNG_UNIFORM_CALLS)
#   define FILE_AND_LINE __FILE__,__LINE__
#else
#   define FILE_AND_LINE
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	RNG class is from Paul Lewis' Lot class from phycas
*/
class RNG
	{
	public:
		RNG(): last_seed_setting(1U), curr_seed(1U), num_seeds_generated(0) {
			useClockToSeed();
		}
		RNG(unsigned rnd_seed) : last_seed_setting(1U), curr_seed(1U), num_seeds_generated(0) {
			if (rnd_seed == 0 || rnd_seed == UINT_MAX)
				useClockToSeed();
			else
				setSeed(rnd_seed);
		}

		// Accessors
		unsigned				getSeed() const {
			return curr_seed;
		}

		unsigned 				getInitSeed() const	{
			return last_seed_setting;
		}



		// Modifiers
		void 					useClockToSeed();
		void 					setSeed(unsigned s)	{
			NEXUS_SHELL_ASSERT(s > 0 && s < UINT_MAX);
			curr_seed = last_seed_setting = s;
		}



		// Utilities
        unsigned                multinomialDraw(const double * probs, unsigned n, double totalProb=1.0);
		unsigned 				sampleUInt(unsigned);
		unsigned				getRandBits(unsigned nbits);

#		if defined(LOG_RNG_UNIFORM_CALLS)
			double 					uniform(const char * file, const int line);
#		else
			double 					uniform();
#		endif

		bool					boolean() {
			return (uniform(FILE_AND_LINE) < 0.5);
		}
	private:    	

		unsigned 				last_seed_setting;
		unsigned				curr_seed;
		unsigned				num_seeds_generated;
	};



class NexusShellException: public std::exception
	{
	public:
		std::string	msg;	/* NxsString to hold message */
		virtual ~NexusShellException() throw()
			{
			}

		NexusShellException(const std::string & s):msg(s){}
		NexusShellException(const char * s):msg(s){}
		
		const char * what () const throw ()
			{
			return msg.empty() ? "Unknown NexusShellException Exception" : msg.c_str();
			}
	};


class InferenceKernelBundle
{
	public:
		InferenceKernelBundle(
			InferenceKernel & k, 
			std::string outPref, 
			std::string nexusBlockName, 
			std::ostream & outputStream,
			std::ostream & errorStream,
			privateBlockReaderFuncPtr pbr)
			:kernel(k),
			outputPrefix(outPref),
			blockName(nexusBlockName),
			outStream(outputStream),
			errStream(errorStream),
			commandReader(pbr)
			{}
			

		InferenceKernel & kernel;
		std::string outputPrefix;
		std::string blockName;
		std::ostream & outStream;
		std::ostream & errStream;
		privateBlockReaderFuncPtr  commandReader;
};

class InferenceKernel
{
	public:
		static RNG rng;
		virtual ~InferenceKernel(){}
	
		virtual void addTaxLabels(const std::vector<std::string>  & ) = 0;
		virtual const std::vector<std::string>  & getTaxLabels() const = 0;
	
		virtual void replaceData(const NxsCharactersBlock * charactersBPtr) = 0;
	
		virtual unsigned addTree(const NxsFullTreeDescription & tree) = 0;
		virtual unsigned getNumTreeDescriptions() const = 0;
		virtual unsigned getNPMats() const = 0;

		
		virtual void runMCMC(int treeNum) = 0;

		virtual double seqLScore(int treeNum) = 0;
		
		virtual void writeTaxa(std::ostream &) const = 0;
		// flags should be or'd bits from NCL's NxsFullTreeDescription::TreeDescFlags bits
		virtual void writeTree(std::ostream &, int treeIndex, int flags=-1, double edgeLenMultiplier=1.0) const = 0;

		virtual void addStateFreqMove(int treeIndex, double wt, int subsetInd, const std::vector<double> & windows, const std::vector<double> & priors) = 0;
		virtual void addRateMatMove(int treeIndex, double wt, int subsetInd, const std::vector<double> & windows, const std::vector<double> & priors) = 0;
		virtual void addEdgeLenMove(int treeIndex, double wt, int subsetInd, double window, double mean) = 0;
		virtual void addPartitionedFreqMove(int treeIndex, double wt,  const std::list< std::vector<double> > & windows, const std::list< std::vector<double> > & priorMeans) = 0;
		
		virtual void setCheckpointFilename(std::string c) = 0;
		virtual void setCheckpointInterval(unsigned long c) = 0;
		virtual void setGamma(double v) = 0;
		virtual void setGeneticCode(NxsGeneticCodesEnum) = 0;
		virtual void setMCMCIterations(unsigned long m) = 0;
		virtual void setCurrMCMCIterations(unsigned long m) = 0;

		virtual void setPinv(double v) = 0;
		virtual void setReportFilename(std::string r) = 0;
		virtual void setReportInterval(unsigned long r) = 0;
		virtual void setRMat(const std::vector<double> & r) = 0;
		virtual void setSampleInterval(unsigned long s) = 0;
		virtual void setSampleTreeFilename(std::string s) = 0;
		virtual void setSampleParamFilename(std::string s) = 0;
		virtual void setStateFreq(const std::vector<double> & f) = 0;
	protected:	
		virtual SeqModel * generateNewSeqModel(LikeStructsBundle &) = 0;
		virtual LikeStructsBundle generateNewLikeStructs(unsigned nCLAs) = 0;


};

#endif
