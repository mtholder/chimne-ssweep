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
#if ! defined (CHIMNE_SWEEP_KERNEL_HPP)
#define CHIMNE_SWEEP_KERNEL_HPP

#include <iostream>
#include <string>
#include <vector>
#include <stack>

#include "ncl/nxscxxdiscretematrix.h"
#include "ncl/nxstreesblock.h"
#include "seq_likelihood.h"
#include "nexus_shell.hpp"


class NxsCharactersBlock;


class MCMCKernel;


class Sample;
class Tree;

class Density
{
	public:
		virtual ~Density(){}
		virtual double getRelativeLnDensity(double x) const = 0;
		virtual double getLnDensityRatio(double numerator, double denominator) const {
			return getRelativeLnDensity(numerator) - getRelativeLnDensity(denominator);
		}
		virtual double getMean() const = 0;
};


/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Encapsulates the gamma probability distribution with shape [ (alpha) and scale parameter (beta).
*/
class GammaDistribution : public Density
	{
	public:
						GammaDistribution(double shape, double scale);
						virtual ~GammaDistribution() {}

		virtual double getRelativeLnDensity(double x) const;
		virtual double getMean() const {
			return alpha*beta;
		}

    protected:
		double alpha;		/* the shape parameter */
		double beta;		/* the scale parameter */
	};





class Node
	{
	public:
		Node()
			:leafData(0L),
			cla(0L),
			pMatArray(0L)
			{
			this->tree = 0L;
			reset();
			}
		~Node() {
		}
		void setPMatArray(PMatArrayObj * ld) {
			pMatArray = ld;
		}
		PMatArrayObj * getPMatArray(void) const{
			return pMatArray;
		}
		void setLeafData(LeafDataObj * ld) {
			leafData = ld;
		}
		LeafDataObj * getLeafData(void) const{
			return leafData;
		}
		CLAObj * getCLA(void) const {
			return cla;
		}
		void setCLA(CLAObj *c) const {
			cla = c;
		}
		/// resets all field except the tree field.
		void reset() {
			this->depth = 0;
			this->sample = 0L;
			this->parent = 0L;
			this->lChild = 0L;
			this->rSib = 0L;
			this->edgeLen = 0.0;
			this->index = 0L;
		}
		void setTree(Tree *t) {
			this->tree = t;
		}
		
		const Tree * setTree() const  {
			return this->tree;
		}
		
		Tree * setTree()  {
			return this->tree;
		}
		void setParent(Node *p)  {
			this->parent = p;
		}
		Node * getParent()  {
			return this->parent;
		}
		const Node * getParentC() const {
			return this->parent;
		}
		void setLChild(Node *p)  {
			this->parent = p;
		}
		Node * getLChild()  {
			return this->lChild;
		}
		const Node * getLChildC() const {
			return this->lChild;
		}
		void setRSib(Node *p)  {
			this->rSib = p;
		}
		Node * getRSib()  {
			return rSib;
		}
		const Node * getRSibC() const {
			return rSib;
		}
		void setRChild(Node *p)  {
			this->addChild(p);
		}
		Node * getRChild()  {
			Node * n = this->lChild;
			if (n == 0L)
				return 0L;
			return n->GetRightmostSibOrSelf(); 
		}
		const Node * getRChildC() const {
			const Node * n = this->lChild;
			if (n == 0L)
				return 0L;
			return n->GetRightmostSibOrSelfC(); 
		}
		Node * getLSib()  {
			if (!this->parent)
				return 0L;
			Node * n = this->parent->lChild;
			if (n == this)
				return 0L;
			while (n && n->rSib != this)
				n = n->rSib;
			return n;
		}
		const Node * getLSibC() const {
			if (!this->parent)
				return 0L;
			Node * n = this->parent->lChild;
			if (n == this)
				return 0L;
			while (n && n->rSib != this)
				n = n->rSib;
			return n;
		}
		Node * GetRightmostSibOrSelf()  {
			Node * n = this;
			while (n->rSib)
				n = n->rSib;
			return n;
		}
		const Node * GetRightmostSibOrSelfC() const {
			const Node * n = this;
			while (n->rSib)
				n = n->rSib;
			return n;
		}
		
		void addChild(Node *c) {
			c->setParent(this);
			c->setRSib(0L);
			Node * n = this->lChild;
			if (n){
				n = n->GetRightmostSibOrSelf();
				n->setRSib(c);
			}	
			else
				this->lChild = c;
		}
		void setEdgeLen(double x) const {
			this->edgeLen = x;
		}
		double getEdgeLen() const {
			return this->edgeLen;
		}
		/// Uses depth (and parent's depth to calculate the edgelengths).
		double calcEdgeLen(bool recurse=false) const;
		void setSample(const Sample *l) {
			this->sample = l;
		}
		const Sample * getSample() const {
			return this->sample;
		}
		
		bool isLeaf() const {
			return this->lChild == 0L;
		}

		bool isSampled() const {
			return this->sample == 0L;
		}

		double getDepth() const;
		void setDepth(double x) {
			NEXUS_SHELL_ASSERT(this->sample == 0L);
			this->depth = x;
		}
		void setIndex(unsigned ind) {
			this->index = ind;
		}
		unsigned getIndex() const;

		// flags should be or'd bits from NCL's NxsFullTreeDescription::TreeDescFlags bits
		void writeNewick(std::ostream &, int flags, double edgeLenMultiplier=1.0) const;

		void addSelfAndDesToPreorder(std::vector<Node *> &p);

	private:
	
		LeafDataObj * leafData;
		mutable CLAObj * cla;
		PMatArrayObj * pMatArray;
		
		double depth;
		mutable double edgeLen; // 
		const Sample * sample;
		Node * parent;
		Node * lChild;
		Node * rSib;
		Tree * tree;
		unsigned index;
	};
	
class Tree
	{
	public:
		static int defaultTreeWritingFlags;
		
		Tree(const std::vector<Sample> &, const NxsFullTreeDescription &);
		~Tree()
			{
			std::vector<Node *>::iterator ndIt = nodeArraysToDelete.begin();
			for (; ndIt != nodeArraysToDelete.end(); ++ndIt)
				delete [] *ndIt;
			}
		
		Node * getNewNode();
		Node * getRandomEdge(double p);
		
		void ultrametricize(double minEdgeLen);
		void setEdgeLenFromDepth();
		
		double getMaxTimeAllowed() {
			return this->maxTimeAllowed;
		}
		std::vector<Node *> GetPreorderTraversal() {
			std::vector<Node *> p;
			if (root)
				root->addSelfAndDesToPreorder(p);
			return p;	

		}
		void setMaxTimeAllowed(double x) {
			this->maxTimeAllowed = x;
		}
		
		// flags should be or'd bits from NCL's NxsFullTreeDescription::TreeDescFlags bits
		std::string getNewick(int flags=-1, double edgeLenMultiplier=1.0) const;
		void writeNewick(std::ostream &, int flags=-1, double edgeLenMultiplier=1.0) const;
		void setTreeLikeArray(FullLAObj *fla) {
			treeLA = fla;
		}
		
		const Node *getRoot() const {
			return root;
		}

		FullLAObj * getTreeLA() const {
			return treeLA;
		}
			
	private:
		
		void refillExtraNodeArray(unsigned);
		
		FullLAObj * treeLA;
		
		
		Node * root;
		std::vector<Node *> leaves;
		std::set<Node *> usedNodes;
		std::stack<Node *> extraNodes;
		std::vector<Node *> nodeArraysToDelete;
		double maxTimeAllowed;
		unsigned index;
	};


class Sample 
	{
	public:
		Sample()
		  :intTime(INT_MAX),
		  depth(DBL_MAX),
		  index(UINT_MAX)
			{}
		Sample(const std::string &label)
		  :intTime(INT_MAX),
		  depth(DBL_MAX),
		  index(UINT_MAX) {
			parseCodedName(label);
		}
		
		std::ostream & writeCodedName(std::ostream &o) const; 
		std::string getCodedName() const;
		void parseCodedName(const std::string & label);

		bool operator!=(const Sample & other) const {
			return ! (*this == other);
		}

		bool operator==(const Sample & other) const {
			if (accession != other.accession)
				return false;
			if (id != other.id)
				return false;
			if (strTime != other.strTime)
				return false;
			return sampleLabel == other.sampleLabel;
		}
		int getTime() const {
			return this->intTime;
		}
		double getDepth() const {
			return this->depth;
		}
		double setDepth(int maxDepth, double multiplier);

		unsigned getIndex() const {
			return this->index;
		}
		void setIndex(unsigned ind) {
			this->index = ind;
		}
		void setLabel(const std::string & l) {
			this->sampleLabel = l;
		}
		const std::string & getLabel() const {
			return this->sampleLabel;
		}
		
	private:
		int intTime; // time since beginning	
		double depth; // distance from the present timepoint.
		std::string strTime; // string from of the intTime field
		std::string id;
		std::string accession;
		std::string sampleLabel;
		unsigned index;
	};

inline std::ostream & operator<<(std::ostream & o, const Sample &l);

class CachedCalculator
{
	public:
		CachedCalculator()
			:dirty(true) 
			{}
		virtual ~CachedCalculator(){}
		virtual void setIsDirty(bool isD) const;
		virtual void setIsDirty(bool isD, int) const;
		bool isDirty() const {
			return dirty;
		}
	protected:
		mutable bool dirty;
		
};

inline void CachedCalculator::setIsDirty(bool isD) const {
	dirty = isD;
}

inline void CachedCalculator::setIsDirty(bool isD, int) const {
	dirty = isD;
}

class ParamExpression
{
	public:
		ParamExpression(){}
		virtual ~ParamExpression(){}
		virtual double getValue() const = 0;
#		if defined(DEBUGGING_LIKELIHOOD)
			virtual std::string getDesc(bool top) const = 0;
#		endif
	private:
		ParamExpression(const ParamExpression &); // not copyable (to avoid slicing)
};


class Parameter
{
	public:
		Parameter()
			:value(1.0),
			minV(0.0),
			maxV(DBL_MAX)
			{}
		void setValue(double v) {
			this->value = v;
			for (std::vector<CachedCalculator *>::iterator lIt = listeners.begin(); lIt != listeners.end(); ++lIt)
				(*lIt)->setIsDirty(true);
			for (std::vector<std::pair<CachedCalculator *, int> >::iterator lIt = intlisteners.begin(); lIt != intlisteners.end(); ++lIt)
				(lIt->first)->setIsDirty(true, lIt->second);
		}
		double getValue() const {
			return this->value;
		}
		void addListener(CachedCalculator *cc) {
			listeners.push_back(cc);
		}
		void addListener(CachedCalculator *cc, int x) {
			intlisteners.push_back(std::pair<CachedCalculator *, int>(cc,x) );
		}
#		if defined(DEBUGGING_LIKELIHOOD)
			std::string getDesc() const {
				return desc;
			}
			void setDesc(std::string d)  {
				desc = d;
			}
			std::string desc;
#		endif
	protected:
		double value;
		double minV;
		double maxV;
		std::vector<CachedCalculator *> listeners;
		std::vector<std::pair<CachedCalculator *, int> > intlisteners;
};

class ParamProduct: public ParamExpression
{
	public:
		ParamProduct(const Parameter & f, const Parameter & s)
			:fp(&f),
			sp(&s),
			fpe(0L),
			spe(0L)
			{}
		ParamProduct(const Parameter & f, const ParamExpression & s)
			:fp(&f),
			sp(0L),
			fpe(0L),
			spe(&s)
			{}

		ParamProduct(const ParamExpression & f, const ParamExpression & s)
			:fp(0L),
			sp(0L),
			fpe(&f),
			spe(&s)
			{}
		ParamProduct(const ParamExpression & f, const Parameter & s)
			:fp(0L),
			sp(&s),
			fpe(&f),
			spe(0L)
			{}
			
		virtual ~ParamProduct(){}
		virtual double getValue() const {
			const double f = (fp ? fp->getValue() : fpe->getValue());
			const double s = (sp ? sp->getValue() : spe->getValue());
			return f*s;
		}
#		if defined(DEBUGGING_LIKELIHOOD)
			virtual std::string getDesc(bool top) const {
			const std::string f = (fp ? fp->getDesc() : fpe->getDesc(false));
			const std::string s = (sp ? sp->getDesc() : spe->getDesc(false));
			std::string ret;
			if (!top)
				ret.append(1, '(');
			ret += f;
			ret.append(1, '*');
			ret += s;
			if (!top)
				ret.append(1, ')');
			return ret;
			}
#		endif
	private:
		const Parameter * fp;
		const Parameter * sp;
		const ParamExpression *fpe;
		const ParamExpression * spe;
		
};

class ParamQuotient: public ParamExpression
{
	public:
		ParamQuotient(const Parameter & f, const Parameter & s)
			:fp(&f),
			sp(&s),
			fpe(0L),
			spe(0L)
			{}
		ParamQuotient(const Parameter & f, const ParamExpression & s)
			:fp(&f),
			sp(0L),
			fpe(0L),
			spe(&s)
			{}
		ParamQuotient(const ParamExpression & f, const ParamExpression & s)
			:fp(0L),
			sp(0L),
			fpe(&f),
			spe(&s)
			{}
		ParamQuotient(const ParamExpression & f, const Parameter & s)
			:fp(0L),
			sp(&s),
			fpe(&f),
			spe(0L)
			{}
			
		virtual ~ParamQuotient(){}
		virtual double getValue() const {
			const double f = (fp ? fp->getValue() : fpe->getValue());
			const double s = (sp ? sp->getValue() : spe->getValue());
			return f/s;
		}
#		if defined(DEBUGGING_LIKELIHOOD)
			virtual std::string getDesc(bool top)  const {
			const std::string f = (fp ? fp->getDesc() : fpe->getDesc(false));
			const std::string s = (sp ? sp->getDesc() : spe->getDesc(false));
			std::string ret;
			if (!top)
				ret.append(1, '(');
			ret += f;
			ret.append(1, '/');
			ret += s;
			if (!top)
				ret.append(1, ')');
			return ret;
			}
#		endif
	private:
		const Parameter * fp;
		const Parameter * sp;
		const ParamExpression * fpe;
		const ParamExpression * spe;
		
};

class ParamSum: public ParamExpression
{
	public:
		ParamSum(const Parameter & f, const Parameter & s)
			:fp(&f),
			sp(&s),
			fpe(0L),
			spe(0L)
			{}
		ParamSum(const Parameter & f, const ParamExpression & s)
			:fp(&f),
			sp(0L),
			fpe(0L),
			spe(&s)
			{}
		ParamSum(const ParamExpression & f, const ParamExpression & s)
			:fp(0L),
			sp(0L),
			fpe(&f),
			spe(&s)
			{}
		ParamSum(const ParamExpression & f, const Parameter & s)
			:fp(0L),
			sp(&s),
			fpe(&f),
			spe(0L)
			{}
			
		virtual ~ParamSum(){}
		virtual double getValue() const {
			const double f = (fp ? fp->getValue() : fpe->getValue());
			const double s = (sp ? sp->getValue() : spe->getValue());
			return f+s;
		}
#		if defined(DEBUGGING_LIKELIHOOD)
			virtual std::string getDesc(bool top) const {
			const std::string f = (fp ? fp->getDesc() : fpe->getDesc(false));
			const std::string s = (sp ? sp->getDesc() : spe->getDesc(false));
			std::string ret;
			if (!top)
				ret.append(1, '(');
			ret += f;
			ret.append(1, '+');
			ret += s;
			if (!top)
				ret.append(1, ')');
			return ret;
			}
#		endif
	private:
		const Parameter * fp;
		const Parameter * sp;
		const ParamExpression * fpe;
		const ParamExpression * spe;
		
};

class ParamWrapper: public ParamExpression
{
	public:
		ParamWrapper(const Parameter & f)
			:fp(&f)
			{}
			
		virtual ~ParamWrapper(){}
		virtual double getValue() const {
			return fp->getValue();
		}
#		if defined(DEBUGGING_LIKELIHOOD)
			virtual std::string getDesc(bool ) const {
				return fp->getDesc();
			}
#		endif
	private:
		const Parameter * fp;		
};

class ConstParamWrapper: public ParamExpression
{
	public:
		ConstParamWrapper(double v)
			:val(v)
			{}
			
		virtual ~ConstParamWrapper(){}
		virtual double getValue() const {
			return val;
		}
#		if defined(DEBUGGING_LIKELIHOOD)
			virtual std::string getDesc(bool ) const {
				NxsString s;
				s += this->val;
				return std::string(s);
			}
#		endif
	private:
		double val;
};

class ParamDifference: public ParamExpression
{
	public:
		ParamDifference(const Parameter & f, const Parameter & s)
			:fp(&f),
			sp(&s),
			fpe(0L),
			spe(0L)
			{}
			
		virtual ~ParamDifference(){}
		virtual double getValue() const {
			const double f = (fp ? fp->getValue() : fpe->getValue());
			const double s = (sp ? sp->getValue() : spe->getValue());
			return f-s;
		}
#		if defined(DEBUGGING_LIKELIHOOD)
			virtual std::string getDesc(bool top) const {
			const std::string f = (fp ? fp->getDesc() : fpe->getDesc(false));
			const std::string s = (sp ? sp->getDesc() : spe->getDesc(false));
			std::string ret;
			if (!top)
				ret.append(1, '(');
			ret += f;
			ret.append(1, '-');
			ret += s;
			if (!top)
				ret.append(1, ')');
			return ret;
			}
#		endif
	private:
		const Parameter * fp;
		const Parameter * sp;
		const ParamExpression *fpe;
		const ParamExpression * spe;
		
};


class ParamGroup: public CachedCalculator
{
	public:
		virtual ~ParamGroup(){}
		
		ParamGroup()
			:CachedCalculator() {
		}
		ParamGroup(unsigned nParam)
			:CachedCalculator() {
			if (nParam > 0)
				resize(nParam);
		}
		virtual void resize(unsigned nParam) {
			rawParams.clear();
			rawParams.resize(nParam);
		}

		const std::vector<Parameter> & getRawParam() const {
			return rawParams;
		}

		std::vector<Parameter> & getRawParamRef() {
			return rawParams;
		}
			
		void setParam(unsigned index, double val) {
			rawParams.at(index).setValue(val);
			this->setIsDirty(true);
		}
		
		unsigned size() const {
			return rawParams.size();
		}

	protected:
		std::vector<Parameter> rawParams;
		
};

class StateFreqParamGroup: public ParamGroup
{
	public:
		StateFreqParamGroup()
			:ParamGroup()
			{}
		StateFreqParamGroup(unsigned nStates)
			:ParamGroup()
			{
			if (nStates > 0)
				resize(nStates);
			}
	

		virtual void resize(unsigned nParam) {
			ParamGroup::resize(nParam);
			for (std::vector<Parameter>::iterator pIt = rawParams.begin(); pIt != rawParams.end(); ++pIt)
				pIt->addListener(this);
		}
	
		const std::vector<double> & getFreqs() const {
			if (this->isDirty()) {
				freqs.clear();
				double t = 0.0;
				for (std::vector<Parameter>::const_iterator pIt = rawParams.begin(); pIt != rawParams.end(); ++pIt) {
					const double v = pIt->getValue();
					t += v;
					freqs.push_back(v);
				}
				for (std::vector<double>::iterator fIt = freqs.begin(); fIt != freqs.end(); ++fIt)
					*fIt /= t;
				this->setIsDirty(false);
			}
			return freqs;
		}

	protected:
		mutable std::vector<double> freqs;
};

class SeqModel: public CachedCalculator
	{
	public:
		SeqModel()
			:rescaleThreshold(200)
			{}
		virtual ~SeqModel(){}

		virtual double lnLikeFromEffectiveRoot(const Node * nd, FullLAObj * full_la, double * subsetLnL);
		virtual void updateNode(const Node &nd, const std::vector<CalcContext> &);

		virtual void writeSample(std::ostream &) {
		}
		virtual void writePrivateCommands(std::ostream & ) const = 0;
		void addTreeLikePtr(FullLAObj *tl) {
			treeLikePtrs.push_back(tl);
		}
		virtual unsigned getNSubsets() const = 0;
		virtual unsigned getNMixtures() const = 0;
		virtual StateFreqParamGroup & getStateFreqRef(int subModelIndex) = 0;
		virtual ParamGroup & getRMatParamGroup(int subModelIndex) = 0;
		virtual std::vector<Parameter> & getRMatParamRef(int subModelIndex) = 0;
		virtual Parameter & getGammaShapeParam(int subModelIndex) = 0;
		virtual Parameter & getPInvarParam(int subModelIndex) = 0;
		int getRescaleThreshold() const {
			return rescaleThreshold;
		}
	protected:
		virtual PMatArrayObj * getNodeUpdatedPMat(const Node &nd, double edgeLen) = 0;
		virtual void propagateQMat() const = 0;

		int rescaleThreshold;
		std::vector<FullLAObj *> treeLikePtrs;
		
	};
	
class CSMCMCChain;

class MoveHistory
{
	public:
		MoveHistory(unsigned nAccepted, unsigned nTried)
			:nMovesAccepted(nAccepted),
			nMovesTried(nTried)
			{}
		MoveHistory(bool accepted)
			:nMovesAccepted(accepted ? 1 : 0),
			nMovesTried(1)
			{}
		
		unsigned nMovesAccepted;
		unsigned nMovesTried;
};
class MCMCMove 
	{
		public:
			MCMCMove()
				:modelIndicator(INT_MAX)
				{}
			virtual ~MCMCMove(){}
			virtual MoveHistory update() = 0;
			virtual void writePrivateCommands(std::ostream & ) const = 0;
			void setModelIndex(int m) {
				modelIndicator = m;
			}
		protected:
			int modelIndicator;
	};
	
double windowMove(double x, double window, double minV, double maxV);
bool metropolisHastingsRule(double lnAcceptRatio);


class SampledSurface
	{
	public:
		virtual ~SampledSurface(){}
		virtual double getRelativeLnDensity() const = 0;
		virtual void revert() const {
		}
	};

class PartitionedSampledSurface
	{
	public:
		virtual ~PartitionedSampledSurface(){}
		virtual std::vector<double> getRelativeLnDensityBySubsets() const = 0;
		virtual void revertSubset(unsigned) const = 0;
	};

class ProbPartitionedSeqGivenGenealogy: public PartitionedSampledSurface
	{
	public:
		ProbPartitionedSeqGivenGenealogy(CSMCMCChain & c)
			:chain(c)
			{}
		virtual ~ProbPartitionedSeqGivenGenealogy(){}
		virtual std::vector<double> getRelativeLnDensityBySubsets() const;
		virtual void revertSubset(unsigned) const;

	protected:
		mutable CSMCMCChain & chain;
		mutable std::vector<double> prevLnLike;
	};


class PositiveParamGroupMove: public MCMCMove
{
	public:
		PositiveParamGroupMove(ParamGroup & params, const std::vector<double> & windows, const std::vector<double> & priors, SampledSurface * likelihoodSurface);
		virtual ~PositiveParamGroupMove()
			{
			for (std::vector<Density *>::iterator m = priors.begin(); m != priors.end(); ++m)
				delete (*m);
			
			}
		virtual MoveHistory update() {
			return updateParamGroup(paramGroup);
		}
		virtual void writePrivateCommands(std::ostream & ) const;
	protected:
		MoveHistory updateParamGroup(ParamGroup &pg);

		ParamGroup & paramGroup;
		std::vector<double> paramWindows;
		std::vector<Density *> priors;
		SampledSurface * likelihoodDensity;

};

class DiscreteProbMove: public PositiveParamGroupMove
{
	public:
		DiscreteProbMove(StateFreqParamGroup & params, const std::vector<double> & windows, const std::vector<double> & priors, SampledSurface * likelihoodSurface);
		virtual ~DiscreteProbMove()
			{
			for (std::vector<Density *>::iterator m = priors.begin(); m != priors.end(); ++m)
				delete (*m);
			
			}
		virtual MoveHistory update() {
			return updateParamGroup(paramGroup);
		}
		virtual void writePrivateCommands(std::ostream & ) const;
	protected:
		StateFreqParamGroup & freqParams;
};

class ChainSampler
{
	public:
		ChainSampler(MCMCKernel & k)
			:checkpointStream(0L),
			reportStream(0L),
			treeStream(0L),
			paramStream(0L),
			kernel(k),
			sampleInterval(1),
			reportInterval(1),
			checkpointInterval(1),
			checkpointTreeFilename("chkpt.tre")
			{}

		void sample(Tree * tree, SeqModel * m, unsigned long iteration, unsigned long totalIter) const;
		void report(Tree * tree, SeqModel * m, unsigned long iteration, unsigned long totalIter) const;
		void checkpoint(Tree * tree, SeqModel * m, unsigned long iteration, unsigned long totalIter, CSMCMCChain & chain) const;
		
		unsigned long getCheckpointInterval() const {
			return checkpointInterval;
		}
		unsigned long getSampleInterval() const {
			return sampleInterval;
		}
		unsigned long getReportInterval() const {
			return reportInterval;
		}
		
		void setCheckpoint(unsigned interval, std::ostream * cmdfile) {
			checkpointInterval = interval;
			if (checkpointStream)
				checkpointStream->flush();
			checkpointStream = cmdfile;
		}
		void setReport(unsigned interval, std::ostream * reportfile) {
			reportInterval = interval;
			if (reportStream)
				reportStream->flush();
			reportStream = reportfile;
		}
		void setSample(unsigned interval, std::ostream * treefile, std::ostream * paramfile) {
			sampleInterval = interval;
			if (treeStream)
				treeStream->flush();
			treeStream = treefile;
			if (paramStream)
				paramStream->flush();
			paramStream = paramfile;
		}
		
		
	protected:
		std::ostream * 	checkpointStream;	
		std::ostream * 	reportStream;	
		std::ostream * 	treeStream;	
		std::ostream * 	paramStream;
		MCMCKernel & kernel;
		unsigned long sampleInterval;
		unsigned long reportInterval;
		unsigned long checkpointInterval;
		std::string checkpointTreeFilename;
};

// Continuous state-space MCMC Chain
class CSMCMCChain
	{
	public:
		CSMCMCChain(double maxD, double timeMultiplier);
		~CSMCMCChain();
		bool isReady() const {
			return (tree && currModel && likeStructs.leafData);
		}
		SeqModel * getModel() {
			return currModel;
		}

		void setModel(SeqModel *m) {
			currModel = m;
		}

		void setCurrSeqLScoreDirty(bool v) {
			currSeqLScoreDirty = v;
		}

		double getCurrSeqLScore() const {
			return currTotalSeqLScore;
		}

		std::vector<double> getCurrSeqPartitionedLScore() const {
			return currSeqLScoreBySite;
		}

		void setCurrSeqLScore(double l) {
			currTotalSeqLScore = l;
		}
		void setCurrSeqPartitionedLScore(const std::vector<double> & l) {
			currSeqLScoreBySite = l;
		}
		
		double seqLScore() const {
			return seqLScore(ndsMoved);
		}
		std::vector<double> seqPartitionedLScore() const {
			this->seqLScore();
			return this->currSeqLScoreBySite;
		}
		
		void run(unsigned long startIteration, unsigned long nIterations, const ChainSampler & s); 
		
		void clearMoves() {
			movesToWeights.clear();
		}
		
		void addMove(MCMCMove * move, double wt) {
			movesToWeights.push_back(std::pair<MCMCMove *, double>(move, wt));
		}
		const std::list< std::pair<MCMCMove *, double> > & getMoves() const {
			return movesToWeights;
		}

		Tree * getTree() {
			return this->tree;
		}
		
		void edgeChangedLen(Node *nd, bool isRevert = false);
		
	protected:
		double seqLScore(std::set<const Node*> dirtyNodes) const;

		void setTree(Tree * t){
			if (this->tree)
				delete this->tree;
			this->tree = t;
		}
		Tree * tree;
		SeqModel * currModel;
		LikeStructsBundle likeStructs;
		std::vector<CLAObj *> spareCLAs;
		double maxDepth;
		double sampleTimeToTreeTimeMultiplier;
		std::set<const Node*> ndsMoved;
		std::list< std::pair<MCMCMove *, double> > movesToWeights;
		mutable std::vector<double> currSeqLScoreBySite;
		mutable double currTotalSeqLScore;
		mutable bool currSeqLScoreDirty;
		
		
		friend class MCMCKernel;
		friend class ChimneSweepKernel;
	};
	
class ProbSeqGivenGenealogy: public SampledSurface
	{
	public:
		ProbSeqGivenGenealogy(CSMCMCChain & c)
			:chain(c)
			{}
		virtual ~ProbSeqGivenGenealogy(){}
		virtual double getRelativeLnDensity() const;
		virtual void revert() const;
		
	protected:
		mutable CSMCMCChain & chain;
		mutable double prevLnLike;
	};



class MCMCKernel : public InferenceKernel
	{
	public:
#		if defined (DEBUGGING_LIKELIHOOD) && DEBUGGING_LIKELIHOOD
			static MCMCKernel * gKernel;
			virtual void debugLikelihoodHook(double, const FullLAObj *, Tree *, SeqModel *) {}
#		endif
		virtual ~MCMCKernel(){}
		virtual const std::vector<std::string>  & getTaxLabels() const {
			taxLabels.clear();
			for (std::vector<Sample>::const_iterator sIt = samples.begin(); sIt != samples.end(); ++sIt) {
				taxLabels.push_back(sIt->getCodedName());
			}
			return taxLabels;
		}
		virtual void addTaxLabels(const std::vector<std::string>  & );
		MCMCKernel();
		const std::vector<Sample> &getSamples() const {
			return this->samples;
		}
		
		unsigned addSamples(const std::vector<Sample> &samples);
		unsigned addTree(const NxsFullTreeDescription &tree) {
			treeDescriptions.push_back(tree);
			return this->getNumTrees();
		}
		unsigned getNumTreeDescriptions() const {
			return this->treeDescriptions.size();
		}
		unsigned getNumTrees() const {
			return this->chains.size();
		}
		void replaceData(const NxsCharactersBlock *);
		void runMCMC(int treeNum);
		double seqLScore(int treeNum);
		void setRMat(const std::vector<double> & r) {
			chains.clear();
			rMatStartingValues = r;
		}
		void setStateFreq(const std::vector<double> & f) {
			chains.clear();
			stateFreqStartingValues = f;
		}
		void setGamma(double v) {
			chains.clear();
			gammaStartingValue = v;
		}
		void setPinv(double v) {
			chains.clear();
			pinvStartingValue = v;
		}
		void ultrametricize(int treeNum, double minEdgeLen);
		void setCurrMCMCIterations(unsigned long m) {
			currMcmcIteration = m;
		}
		void setMCMCIterations(unsigned long m) {
			mcmcIterations = m;
		}
		
		void writeTaxa(std::ostream &) const;
		// flags should be or'd bits from NCL's NxsFullTreeDescription::TreeDescFlags bits
		void writeTree(std::ostream &, int treeIndex, int flags=-1, double edgeLenMultiplier=1.0) const;
		// flags should be or'd bits from NCL's NxsFullTreeDescription::TreeDescFlags bits
		void writeTree(std::ostream &o, const Tree * t, int flags=-1, double edgeLenMultiplier=1.0) const;

		void addStateFreqMove(int treeIndex, double wt, int subModelIndex, const std::vector<double> & windows, const std::vector<double> & priors);
		void addRateMatMove(int treeIndex, double wt, int subModelIndex, const std::vector<double> & windows, const std::vector<double> & priors);

		virtual void writeCheckpointPreExecSettings(std::ostream & out) const;
		virtual void writeCheckpointPostExecSettings(std::ostream & ) const {
			// pass
		}
		void writeMCMCCheckpointInfo(std::ostream & out) const;
		
		void setCheckpointInterval(unsigned long c) {
			this->checkpointInterval = c;
		}
		void setGeneticCode(NxsGeneticCodesEnum gCode) {
			this->geneticCode = gCode;
		}
		void setSampleInterval(unsigned long s) {
			this->sampleInterval = s;
		}
		void setReportInterval(unsigned long r) {
			this->reportInterval = r;
		}
		void setSampleTreeFilename(std::string s) {
			this->sampleTreeFilename = s;
		}
		void setSampleParamFilename(std::string s) {
			this->sampleParamFilename = s;
		}
		void setReportFilename(std::string r) {
			this->reportFilename = r;
		}
		void setCheckpointFilename(std::string c) {
			this->checkpointFilename = c;
		}

		virtual unsigned getNPMats() const;
	protected:
		virtual void addLikelihoodStructsToChain(CSMCMCChain &, unsigned nCLAs = 0);
		LikeStructsBundle generateNewLikeStructs(unsigned nCLAs);
		const Tree * getTreeC(int treeNum) const;
		Tree * getTree(int treeNum);
		void clearData() {
			this->charMatrix.Initialize(0L, true);
			this->charMatrixInitialized = false;
		}
		void deleteMoves() {
			for (std::list<MCMCMove *>::iterator m = movesToDelete.begin(); m != movesToDelete.end(); ++m)
				delete (*m);
			movesToDelete.clear();
		}
		void deleteSurfaces() {
			for (std::list<SampledSurface *>::iterator m = surfacesToDelete.begin(); m != surfacesToDelete.end(); ++m)
				delete (*m);
			surfacesToDelete.clear();
			for (std::list<PartitionedSampledSurface *>::iterator p = partSurfacesToDelete.begin(); p != partSurfacesToDelete.end(); ++p)
				delete (*p);
			partSurfacesToDelete.clear();
		}

		void clearTrees() {
			deleteSurfaces();
			deleteMoves();
			chains.clear();
		}
		void clearTreeDescriptions() {
			treeDescriptions.clear();
		}
		virtual void initCharMatrix() {
		}
		
		void rescaleSampleDepths();
		void taxStatus(std::ostream &) const;
		CSMCMCChain & getChainForTreeIndex(int treeIndex);

		std::vector<Sample> samples;
		NxsCXXDiscreteMatrix charMatrix;
		std::vector<NxsFullTreeDescription> treeDescriptions;
		std::vector<CSMCMCChain> chains;
		std::vector< std::vector<int> >  paramIdentityMatrix;
		int nRateCats;
		
		double maxDepth;
		double sampleTimeToTreeTimeMultiplier;
		int maxSampleTime;
		std::vector<double> stateFreqStartingValues;
		std::vector<double> rMatStartingValues;
		double gammaStartingValue;
		double pinvStartingValue;
		std::list<MCMCMove *> movesToDelete;
		std::list<SampledSurface *> surfacesToDelete;
		std::list<PartitionedSampledSurface *> partSurfacesToDelete;
		unsigned long mcmcIterations;
		unsigned long currMcmcIteration;
		unsigned long sampleInterval;
		unsigned long reportInterval;
		unsigned long checkpointInterval;
		std::string sampleTreeFilename;
		std::string sampleParamFilename;
		std::string reportFilename;
		std::string checkpointFilename;
		mutable std::vector<std::string> taxLabels;
		bool charMatrixInitialized;
		NxsGeneticCodesEnum geneticCode;

	};



inline std::ostream & operator<<(std::ostream & o, const Sample &l) {
	l.writeCodedName(o);
	return o;
}
	

inline double Node::getDepth() const 
{
	return (this->sample ? this->sample->getDepth() : this->depth);
}

inline unsigned Node::getIndex() const
{
	return (this->sample ? this->sample->getIndex() : this->index);
}

#if ! defined(EDGE_LEN_FROM_DEPTH)
	inline double Node::calcEdgeLen(bool) const 
	{
		return this->getEdgeLen();
	}
#endif //defined(EDGE_LEN_FROM_DEPTH)


namespace chim {

class ChimneSweepSeqModel: public SeqModel
	{
	public:
		ChimneSweepSeqModel(DSCTModelObj & model, ASRVObj * asrvObj, const std::vector< std::vector<int> > & paramIdentityMatrix);
		virtual ~ChimneSweepSeqModel(){}
		DSCTModelObj & getCCoreModel() const {
			return this->cCoreModel;
		}
		StateFreqParamGroup & getStateFreqRef(int) {
			return stateFreq;
		}
		ParamGroup & getRMatParamGroup(int) {
			return rMatParams;
		}
		std::vector<Parameter> & getRMatParamRef(int) {
			return rMatParams.getRawParamRef();
		}
		Parameter & getGammaShapeParam(int) {
			return gammaShape;
		}
		Parameter & getPInvarParam(int) {
			return pInvar;
		}
		void writePrivateCommands(std::ostream & ) const;
		virtual std::string getPaupLSetOptions() const;
		virtual unsigned getNMixtures() const {
			return (asrv ? asrv->n : 0);
		}
		virtual unsigned getNSubsets() const {
			return 1;
		}


	protected:
		virtual PMatArrayObj * getNodeUpdatedPMat(const Node &nd, double edgeLen);
		virtual void propagateQMat() const;
		virtual void propagateRateHet() const;
	
		mutable DSCTModelObj & cCoreModel; //alias to DSCTModelObj * object.
		ASRVObj * asrv;
		ParamGroup rMatParams;
		typedef std::vector<Parameter * > QRateRow;
		typedef std::vector<QRateRow> QRateMatrix;
		QRateMatrix qRateMatrix;
		StateFreqParamGroup stateFreq;
		Parameter gammaShape;
		Parameter pInvar;
	};

class ChimneSweepKernel : public MCMCKernel
	{
	public:
#		if defined (DEBUGGING_LIKELIHOOD) && DEBUGGING_LIKELIHOOD
			void debugLikelihoodHook(double, const FullLAObj *, Tree *, SeqModel *);
#		endif
		virtual void addTaxLabels(const std::vector<std::string>  & );
		ChimneSweepKernel():MCMCKernel(){}
		virtual ~ChimneSweepKernel(){}
		void ultrametricize(int treeNum, double minEdgeLen);
		void setEdgeLenFromDepth(int treeNum);

		virtual void addEdgeLenMove(int treeIndex, double wt, int subsetInd, double window, double mean);
		virtual void addPartitionedFreqMove(int treeIndex, double wt, const std::list< std::vector<double> > & windows, const std::list< std::vector<double> > & priorMeans);
		void writeCheckpointPostExecSettings(std::ostream & out) const;
	protected:
		void addLikelihoodStructsToChain(CSMCMCChain & chain, unsigned nCLAs);
		virtual SeqModel * generateNewSeqModel(LikeStructsBundle &);
	};

} // namespace chim 


#endif

