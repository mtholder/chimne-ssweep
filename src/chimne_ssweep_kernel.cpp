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
#include <fstream>
#include <sstream>
#include <climits>
#include <algorithm>
#include <cmath>
#include "chimne_ssweep_kernel.hpp"
#include "chimne_ssweep_input.hpp" // for exception
#if defined (DEBUGGING_LIKELIHOOD) && DEBUGGING_LIKELIHOOD
#	include <fstream>
#endif

using namespace std;


#if defined (DEBUGGING_LIKELIHOOD) && DEBUGGING_LIKELIHOOD
	MCMCKernel * MCMCKernel::gKernel = NULL;
#endif



int Tree::defaultTreeWritingFlags = NxsFullTreeDescription::NXS_EDGE_LENGTH_UNION;


std::string Sample::getCodedName() const
{
	ostringstream s;
	writeCodedName(s);
	return s.str();
}




//==============================================================================
//	density.cpp taken from the probability_distribution file in Phycas
//==============================================================================


/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the shape and scale parameters to the specified values.
*/
GammaDistribution::GammaDistribution(
  double shape,		/* the shape parameter */
  double scale)		/* the scale parameter */
  	{
	NEXUS_SHELL_ASSERT(shape > 0.0);
	NEXUS_SHELL_ASSERT(scale > 0.0);
	alpha = shape;
	beta = scale;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The probability density function of the gamma is 
|>
|	       x^(alpha-1) exp(-x/beta)
|	f(x) = ------------------------
|	       beta^alpha Gamma(alpha)
|<
|	This function returns the natural log of the density function at `x':
|>
|	ln[f(x)] = (alpha - 1)*ln(x) - x/beta - alpha*ln(beta) - lnGamma(alpha)
|<
|	The sum of the last two terms are precalculated and available as the variable `ln_const'.
|	Returns -DBL_MAX if PDF is zero, which happens if x < 0.0, or if (x = 0.0 and alpha > 1.0).
*/
double GammaDistribution::getRelativeLnDensity(
  double x) const/* the value for which the density function is to be evaluated */
	{
	if (x < 0.0 || (x == 0.0 && alpha > 1.0))
		return -DBL_MAX;
	else
		{
		double term1 = (alpha - 1.0)*std::log(x);
		double term2 = -x/beta;
		double lnpdf = term1 + term2;

		return lnpdf;
		}
	}

//==============================================================================
//		mcmc_move.cpp taken from the Lot class in phycas
//==============================================================================

double ProbSeqGivenGenealogy::getRelativeLnDensity() const 
{
	prevLnLike = chain.getCurrSeqLScore();
	const double proposedLnLike = chain.seqLScore();
	return proposedLnLike - prevLnLike;
}

void ProbSeqGivenGenealogy::revert() const 
{
	chain.setCurrSeqLScore(prevLnLike);
	chain.setCurrSeqLScoreDirty(true);	
}



double windowMove(double x, double window, double minV, double maxV)
{
	const double u = InferenceKernel::rng.uniform(FILE_AND_LINE);
	double step = window*(u - 0.5);
	x += step;
	while (x < minV || x > maxV){
		if (x < minV)
			x = 2*minV + - x;
		else
			x = 2*maxV - x;
	}
	return x;
}

bool metropolisHastingsRule(double lnAcceptRatio)
{
	if (lnAcceptRatio > 0)
		return true;
	if (lnAcceptRatio < -30.0)
		return false;
	return InferenceKernel::rng.uniform(FILE_AND_LINE) < exp(lnAcceptRatio);
}



PositiveParamGroupMove::PositiveParamGroupMove(ParamGroup & params, const std::vector<double> & windows, const std::vector<double> & priorMeans,  SampledSurface * likelihoodSurface)
  :paramGroup(params),
  paramWindows(windows),
  likelihoodDensity(likelihoodSurface)
{
	NEXUS_SHELL_ASSERT(windows.size() == priorMeans.size());
	for (unsigned i = 0; i < priorMeans.size(); ++i) {
		Density * pd = new GammaDistribution(priorMeans[i], 1.0);
		priors.push_back(pd);
	}
}

DiscreteProbMove::DiscreteProbMove(StateFreqParamGroup & params, const std::vector<double> & windows, const std::vector<double> & priorMeans,  SampledSurface * likelihoodSurface)
  :PositiveParamGroupMove(params,  windows,  priorMeans,  likelihoodSurface),
  freqParams(params)
{
	NEXUS_SHELL_ASSERT(params.getFreqs().size() == windows.size());
	NEXUS_SHELL_ASSERT(priorMeans.size() == windows.size());
}

MoveHistory PositiveParamGroupMove::updateParamGroup(ParamGroup &pg) 
{
	NEXUS_SHELL_ASSERT(likelihoodDensity);
	const std::vector<Parameter> & rawParams = pg.getRawParam();
	NEXUS_SHELL_ASSERT(rawParams.size() == this->paramWindows.size());
	NEXUS_SHELL_ASSERT(priors.size() == this->paramWindows.size());
	std::vector<double> prevValues(rawParams.size(), 0.0);
	std::vector<double> proposedValues(rawParams.size(), 0.0);
	double lnPriorRatio = 0.0;
	for (unsigned i = 0; i < rawParams.size(); ++i) {
		double x = rawParams[i].getValue();
		prevValues[i] = x;
		double y = windowMove(x, this->paramWindows.at(i), 0, DBL_MAX);
		proposedValues[i] = y;
		Density * prior = this->priors[i];
		lnPriorRatio += prior->getLnDensityRatio(y, x);
	}
	for (unsigned i = 0; i < rawParams.size(); ++i) {
		pg.setParam(i, proposedValues[i]);
	}
	double lnLikeRatio = likelihoodDensity->getRelativeLnDensity();
	const bool shouldAccept = metropolisHastingsRule(lnPriorRatio+lnLikeRatio);
	if (!shouldAccept) {
		for (unsigned i = 0; i < rawParams.size(); ++i) {
			pg.setParam(i, prevValues[i]);
		}
		likelihoodDensity->revert();
	}
	return MoveHistory(shouldAccept);
}

void DiscreteProbMove::writePrivateCommands(std::ostream & s) const 
{
	if (modelIndicator != INT_MAX) {
		s << " modelnumber = " <<  NxsString::GetEscapedInt(1 + this->modelIndicator) ;
	}
	s << " freqWindow = (";
	for (std::vector<double>::const_iterator wIt = this->paramWindows.begin(); wIt !=  this->paramWindows.end(); ++wIt)
		s << ' ' << NxsString::GetEscapedDouble(*wIt) <<' ';
	s << ") freqPriorMean = (";
	for (std::vector<Density *>::const_iterator dIt = this->priors.begin(); dIt !=  this->priors.end(); ++dIt)
		s << ' ' << NxsString::GetEscapedDouble((*dIt)->getMean()) << ' ';
	s << "); ";
	
	
}

void PositiveParamGroupMove::writePrivateCommands(std::ostream & s) const 
{
	s << " rmatWindow = (";
	for (std::vector<double>::const_iterator wIt = this->paramWindows.begin(); wIt !=  paramWindows.end(); ++wIt)
		s << ' ' << NxsString::GetEscapedDouble(*wIt) <<' ';
	s << ") rmatPriorMean = (";
	for (std::vector<Density *>::const_iterator dIt = this->priors.begin(); dIt !=  priors.end(); ++dIt)
		s << ' ' << NxsString::GetEscapedDouble((*dIt)->getMean()) << ' ';
	s << "); ";
}




//==============================================================================
//		node.cpp
//==============================================================================

#if defined(EDGE_LEN_FROM_DEPTH)
	double Node::calcEdgeLen(bool recurse) const 
	{
			if (this->parent) {
				this->edgeLen = this->parent->getDepth() - this->getDepth();
				//cout << this->parent->getDepth() << " - " <<  this->getDepth() << " = " << this->edgeLen << '\n';
			}
			else {
				NEXUS_SHELL_ASSERT(this->tree);
				this->edgeLen = this->tree->getMaxTimeAllowed() - this->getDepth();
			}
			if (recurse) {
				const Node * nd = this->getLChildC();
				while (nd) {
					nd->calcEdgeLen(true);
					nd = nd->getRSibC();
				}
			}	
		return this->getEdgeLen();
	}
#endif //defined(EDGE_LEN_FROM_DEPTH)

void Node::writeNewick(std::ostream &o, int flags, double edgeLenMultiplier) const 
{
	if (this->lChild) {
		o << '(';
		const Node * nd = this->lChild;
		while (nd) {
			if (nd != this->lChild)
				o << ',';
			nd->writeNewick(o, flags, edgeLenMultiplier);
			nd = nd->rSib;
		}
		o << ')';
		if (flags&NxsFullTreeDescription::NXS_HAS_INTERNAL_NAMES_BIT)
			o << 1 + this->getIndex() ;

	}
	else 
		o << 1 + this->getIndex() ;

	if (this->parent && (flags&NxsFullTreeDescription::NXS_EDGE_LENGTH_UNION))
		o << ':' << edgeLenMultiplier*this->getEdgeLen();
}


void Node::addSelfAndDesToPreorder(std::vector<Node *> &p)
{
	p.push_back(this);
	Node * currCh = this->lChild;
	while (currCh) {
		currCh->addSelfAndDesToPreorder(p);
		currCh = currCh->rSib;
	}
}

//==============================================================================
//		tree.cpp
//==============================================================================
Tree::Tree(
	const std::vector<Sample> &sampleVec, 
	const NxsFullTreeDescription &treeDesc)
	:treeLA(0L)
{
	this->maxTimeAllowed = DBL_MAX;
	this->refillExtraNodeArray(3*sampleVec.size());
	const std::string &newick = treeDesc.GetNewick();
	//std::cout << "tree string is " << newick << '\n';

	std::string n;
	n.reserve(2 + newick.length());
	n = newick;
	n.append(1, ';');
	istringstream newickstream(n);
	NxsToken token(newickstream);
	token.SetEOFAllowed(false);
	double lastFltEdgeLen;
	long currTaxNumber;
	unsigned nodeIndex = sampleVec.size();
	token.GetNextToken();
	NEXUS_SHELL_ASSERT(token.Equals("("));
	this->root =  this->getNewNode();
	Node * currNd = this->root;
	bool prevInternalOrLength;
	bool currInternalOrLength = false;
	std::set<unsigned> sampleInds;
	std::set<unsigned> internalInds;
	for (;;) {
		if (token.Equals(";")) {
			if (currNd != root)
				throw NexusShellException("Semicolon found before the end of the tree description.  This means that more \"(\" characters  than \")\"  were found.");
			break;
		}
		const NxsString & tstr = token.GetTokenReference();
		const char * t = tstr.c_str();
		bool handled;
		handled = false;
		prevInternalOrLength = currInternalOrLength;
		currInternalOrLength = false;
		if (tstr.length() == 1) {
			handled = true;
			if (t[0] == '(') {
				currNd->setIndex(nodeIndex++);
				Node * tmpNode = this->getNewNode();
				currNd->addChild(tmpNode);
				currNd = tmpNode;
			}
			else if (t[0] == ')') {
				currNd = currNd->getParent();
				NEXUS_SHELL_ASSERT(currNd);
				currInternalOrLength = true;
			}
			else if (t[0] == ':') {
				token.SetLabileFlagBit(NxsToken::hyphenNotPunctuation); // this allows us to deal with sci. not. in branchlengths (and negative branch lengths).
				token.GetNextToken();
				t = token.GetTokenReference().c_str();
				if (!NxsString::to_double(t, &lastFltEdgeLen)) {
					NxsString emsg;
					emsg << "Expecting a number as a branch length. Found " << tstr;
					throw NxsException(emsg, token);
				}
				currNd->setEdgeLen(lastFltEdgeLen);
				currInternalOrLength = true;
			}
			else if (t[0] == ',') {
				currNd = currNd->getParent();
				NEXUS_SHELL_ASSERT(currNd);
				Node * tmpNode = this->getNewNode();
				currNd->addChild(tmpNode);
				currNd = tmpNode;
			}
			else 
				handled = false;
		}
		if (!handled) {
			bool wasReadAsNumber = NxsString::to_long(t, &currTaxNumber);
			if (!wasReadAsNumber || currTaxNumber < 1) {
				NxsString emsg;
				emsg << "Expecting a positive number as the taxon number in the tree description. Found " << tstr;
				throw NxsException(emsg, token);
			}
			unsigned ind = (unsigned)(currTaxNumber - 1);
			if (!prevInternalOrLength) {
				if (ind >= sampleVec.size()) {
					NxsString emsg;
					emsg << "Expecting a taxon number less than < " << sampleVec.size() + 1 << ". Found " << tstr;
					throw NxsException(emsg, token);
				}
				if (sampleInds.end() != sampleInds.find(ind)) {
					NxsString emsg;
					emsg << "Taxon number " << ind + 1 << " repeated in the tree description";
					throw NxsException(emsg, token);
				}
				sampleInds.insert(ind);
				currNd->setSample(&sampleVec[ind]);
				if (currNd->getLChild() == NULL) {
					while (ind >= this->leaves.size())
						this->leaves.push_back(0L);
					this->leaves[ind] = currNd;
				}
			}
			else {
				if (ind < sampleVec.size()) {
					NxsString emsg;
					emsg << "Expecting an internal node number to greater than < " << sampleVec.size()  << ". Found " << tstr;
					throw NxsException(emsg, token);
				}
				if (internalInds.end() != internalInds.find(ind)) {
					NxsString emsg;
					emsg << "Internal node number " << ind + 1 << " repeated in the tree description";
					throw NxsException(emsg, token);
				}
				internalInds.insert(ind);
				currNd->setIndex(ind);
			}
		}
		token.GetNextToken();
	}
	while (this->leaves.size() < sampleVec.size())
		this->leaves.push_back(0L);
	if (sampleInds.size() != sampleVec.size())
		std::cerr << "chimne_sweep warning: There are " << sampleInds.size() << " taxa in the tree, but " << sampleVec.size() << " taxa in memory.\n";
	//std::cout << "root depth = " << root->getDepth() << '\n';
}

Node * Tree::getRandomEdge(double p) 
{
	NEXUS_SHELL_ASSERT(p >= 0.0);
	NEXUS_SHELL_ASSERT(p <= 1.0);
	
	std::vector<Node *> preorder = GetPreorderTraversal();
	const unsigned nNodes = preorder.size();
	NEXUS_SHELL_ASSERT(nNodes > 0);
	if (nNodes < 2 || p < 0.0)
		return 0L;
	const double nnmo = (float)(nNodes-1); // we subtract one, because we do not want to return the root (which has no edge length).
	const unsigned offset = 1 + (unsigned)(nnmo*p);// cast to unsigned rounds down except if p is exactly -- in this case we return the last node 
	if (nNodes <= offset)
		return preorder[nNodes-1];
	return preorder[offset];
}
		
void Tree::refillExtraNodeArray(unsigned numToAlloc)
{
	if (numToAlloc == 0)
		return;
	Node * ndArr = new Node[numToAlloc];
	this->nodeArraysToDelete.push_back(ndArr);
	for (unsigned i = 0; i < numToAlloc; ++i) {
		Node * nd = &(ndArr[i]);
		nd->setTree(this);
		this->extraNodes.push(nd);
	}
}


Node * Tree::getNewNode()
{
	if (this->extraNodes.empty())
		this->refillExtraNodeArray(leaves.size());
	Node * nd = this->extraNodes.top(); 
	this->extraNodes.pop();
	nd->reset();
	this->usedNodes.insert(nd);
	return nd;
}


////////////////////////////////////////////////////////////////////////////////
// Use branch lengths and depth to set reset depths for all nodes such that
//	the distance from an internal node to each of its descendants is as small as
//	possible without being shorter than the sum of branch lengths on the tree.
////////////////////////////////////////////////////////////////////////////////
void Tree::ultrametricize(double minEdgeLen)
{

	std::set<Node *>::iterator ndIt = this->usedNodes.begin();
	for (; ndIt != this->usedNodes.end(); ++ndIt)  {
		if (!(*ndIt)->isLeaf())
			(*ndIt)->setDepth(0.0);
	}

	for (std::vector<Node *>::const_iterator leafIt = leaves.begin(); leafIt != leaves.end(); ++leafIt) {
		Node * nd = *leafIt; 
		if (nd) {
			double d = nd->getDepth();
			double e = std::max(minEdgeLen, nd->getEdgeLen());
			double pd = d + e;
			nd = nd->getParent();
			while (nd) {
				if (nd->getDepth() > pd)
					break;
				//cout << "setting depth " << pd << '\n';
				nd->setDepth(pd);
				pd += nd->getEdgeLen();
				nd = nd->getParent();
			}
		}
	}
}


std::string Tree::getNewick(int flags, double edgeLenMultiplier) const
{
	std::ostringstream sstr;
	this->writeNewick(sstr, flags, edgeLenMultiplier);
	return sstr.str();
}

void Tree::writeNewick(std::ostream &o, int flags, double edgeLenMultiplier) const 
{
	if (!root)
		return;
	int f = (flags >= 0 ? flags : Tree::defaultTreeWritingFlags);
	root->writeNewick(o, f, edgeLenMultiplier);
}


////////////////////////////////////////////////////////////////////////////////
// Use node depths to set edgeLen field of all nodes.
////////////////////////////////////////////////////////////////////////////////
void Tree::setEdgeLenFromDepth()
{
	if (root)
		root->calcEdgeLen(true);
}






//==============================================================================
//		sample.cpp
//==============================================================================


double Sample::setDepth(int maxDepth, double multiplier)
{
	int diff = maxDepth - this->intTime;
	this->depth = multiplier* (double) diff;
	//std::cout << this->depth << '\n';
	return this->depth;
}


//==============================================================================
//	SeqModel.cpp
//==============================================================================
const Node * getCalcNodeAndPath(const Node *nd, double *pathLen)
{
	NEXUS_SHELL_ASSERT(nd && pathLen);
	nd->calcEdgeLen(false);
	*pathLen = nd->getEdgeLen();
	const Node * c = nd->getLChildC();
	while ((c != NULL) && (c->getRSibC() == NULL)) {
		c->calcEdgeLen(false);
		*pathLen += c->getEdgeLen();
		nd = c;
		c = nd->getLChildC();
	}
	return nd;
}		

void SeqModel::updateNode(const Node &nd, const std::vector<CalcContext> & ccs)
{
	const Node * f = nd.getLChildC();
	if (f == 0L)
		return;
	const Node * s = f->getRSibC();
	if (s == 0L)
		return; // we don't calc. cond. likelihood at nodes of degree 2.

	double fEdgeLen;
	const Node * firstCalcNd = getCalcNodeAndPath(f, &fEdgeLen);
	NEXUS_SHELL_ASSERT(firstCalcNd);
	PMatArrayObj * fPMat =  getNodeUpdatedPMat(*firstCalcNd, fEdgeLen);

	double sEdgeLen;
	const Node * secCalcNd = getCalcNodeAndPath(s, &sEdgeLen);
	NEXUS_SHELL_ASSERT(secCalcNd);
	PMatArrayObj * sPMat =  getNodeUpdatedPMat(*secCalcNd, sEdgeLen);

	CLAObj * destinationCLA = nd.getCLA();
	NEXUS_SHELL_ASSERT(destinationCLA);
	
#	if defined(DEBUGGING_LIKELIHOOD) && DEBUGGING_LIKELIHOOD
		std::cerr << "Updating node " << (long)(&nd) << '\n';
		std::cerr << " leftchild " << (long)(f) << '\n';
		std::cerr << "   leftchild calc " << (long)(firstCalcNd) << '\n';
		std::cerr << "   brlen " << fEdgeLen << '\n';
		if (firstCalcNd->isLeaf()) {
			std::cerr << "   leftchild calc leaf ";
			firstCalcNd->getSample()->writeCodedName(std::cerr);
			std::cerr << '\n';
		}
		std::cerr << " nextchild " << (long)(s) << '\n';
		std::cerr << "   nextchild calc " << (long)(secCalcNd) << '\n';
		std::cerr << "   brlen " << sEdgeLen << '\n';
		if (secCalcNd->isLeaf()) {
			std::cerr << "   nextchild calc leaf ";
			secCalcNd->getSample()->writeCodedName(std::cerr);
			std::cerr << '\n';
		}
#	endif
	if (firstCalcNd->isLeaf()) {
		LeafDataObj * fLeafData = firstCalcNd->getLeafData();
		NEXUS_SHELL_ASSERT(fLeafData);
		if (secCalcNd->isLeaf()) {
			LeafDataObj * sLeafData = secCalcNd->getLeafData();
			NEXUS_SHELL_ASSERT(sLeafData);
			for (std::vector<CalcContext>::const_iterator cc = ccs.begin(); cc != ccs.end(); ++cc)
				do_two_tip_cla(fLeafData, fPMat, sLeafData, sPMat, destinationCLA, *cc);
		}
		else {
			const CLAObj * sCLA = secCalcNd->getCLA();
			NEXUS_SHELL_ASSERT(sCLA);
			for (std::vector<CalcContext>::const_iterator cc = ccs.begin(); cc != ccs.end(); ++cc)
				do_one_tip_cla(sCLA, sPMat, fLeafData, fPMat, destinationCLA, *cc);
		}
	}
	else {
		const CLAObj * fCLA = firstCalcNd->getCLA();
		NEXUS_SHELL_ASSERT(fCLA);
		if (secCalcNd->isLeaf()) {
			LeafDataObj * sLeafData = secCalcNd->getLeafData();
			NEXUS_SHELL_ASSERT(sLeafData);
			for (std::vector<CalcContext>::const_iterator cc = ccs.begin(); cc != ccs.end(); ++cc)
				do_one_tip_cla(fCLA, fPMat, sLeafData, sPMat, destinationCLA, *cc);
		}
		else {
			const CLAObj * sCLA = secCalcNd->getCLA();
			NEXUS_SHELL_ASSERT(sCLA);
			for (std::vector<CalcContext>::const_iterator cc = ccs.begin(); cc != ccs.end(); ++cc)
				do_internal_cla(fCLA, fPMat, sCLA, sPMat, destinationCLA, *cc);
		}
	}

	const Node *other = s->getRSibC();
	while (other) {
		double oEdgeLen;
		const Node * otherCalcNd = getCalcNodeAndPath(other, &oEdgeLen);
		NEXUS_SHELL_ASSERT(otherCalcNd);

#		if defined(DEBUGGING_LIKELIHOOD) && DEBUGGING_LIKELIHOOD
			std::cerr << " anotherchild " << (long)(other) << '\n';
			std::cerr << "   anotherchild calc " << (long)(otherCalcNd) << '\n';
			std::cerr << "   brlen " << oEdgeLen << '\n';
			if (otherCalcNd->isLeaf()) {
				std::cerr << "   anotherchild calc leaf ";
				otherCalcNd->getSample()->writeCodedName(std::cerr);
				std::cerr << '\n';
			}
#		endif

		PMatArrayObj * oPMat =  getNodeUpdatedPMat(*otherCalcNd, oEdgeLen);
		if (otherCalcNd->isLeaf()) {
			LeafDataObj * oLeafData = otherCalcNd->getLeafData();
			NEXUS_SHELL_ASSERT(oLeafData);
			for (std::vector<CalcContext>::const_iterator cc = ccs.begin(); cc != ccs.end(); ++cc)
				do_add_leaf_to_cla(oLeafData, oPMat, destinationCLA, *cc);
		}
		else {
			const CLAObj * oCLA = otherCalcNd->getCLA();
			NEXUS_SHELL_ASSERT(oCLA);
			for (std::vector<CalcContext>::const_iterator cc = ccs.begin(); cc != ccs.end(); ++cc)
				do_add_internal_to_cla(oCLA, oPMat, destinationCLA, *cc);
		}
		other = other->getRSibC();
	}
}


double SeqModel::lnLikeFromEffectiveRoot(const Node * nd, FullLAObj * full_la, double * subsetLnL)
{
	double scratch;
	const Node * rootCalcNode = getCalcNodeAndPath(nd, &scratch);
	NEXUS_SHELL_ASSERT(rootCalcNode);
	const CLAObj * rootCLA = rootCalcNode->getCLA();
	NEXUS_SHELL_ASSERT(rootCLA);
	
	NEXUS_SHELL_ASSERT(full_la);
	CLAObj* full_cla = full_la->full_cla;
	copy_cla(full_cla, rootCLA);
	const double lnL = partitioned_ln_likelihood(full_la, subsetLnL);
	return lnL;
}


//==============================================================================
//		chain_sampler.cpp
//==============================================================================
void ChainSampler::sample(Tree * tree, SeqModel * m, unsigned long iteration, unsigned long /*totalIter*/) const
{
	if (treeStream && tree) {
		*treeStream << "Tree rep." << iteration << " = [&U] ";
		tree->writeNewick(*treeStream, -1, 1.0);
		*treeStream << ";\n";
		treeStream->flush();
	}
	if (paramStream && m) {
		*paramStream << iteration << "\t";
		m->writeSample(*paramStream);
		*paramStream << "\n";
		paramStream->flush();
	}

}

void ChainSampler::report(Tree * , SeqModel * , unsigned long iteration, unsigned long totalIter) const
{
	if (reportStream == 0L)
		return;
	*reportStream << iteration << " iteration(s) completed. " << totalIter-iteration << " to go...\n";
	reportStream->flush();
}

void ChainSampler::checkpoint(Tree * tree, SeqModel * m, unsigned long iteration, unsigned long totalIter, CSMCMCChain & chain) const
{
	if (!checkpointStream)
		return;
	std::ofstream treef(checkpointTreeFilename.c_str());
	kernel.writeTaxa(treef);
	kernel.writeTree(treef, tree, -1);
	treef.close();
	*checkpointStream << "Set seed = " << InferenceKernel::rng.getSeed() << "; ";
	if (totalIter > iteration)
		*checkpointStream << "Set currmcmciteration = "<< iteration << " mcmcIterations = " << totalIter << "; ";
	kernel.writeCheckpointPreExecSettings(*checkpointStream);
	m->writePrivateCommands(*checkpointStream);
	kernel.writeMCMCCheckpointInfo(*checkpointStream);
	
	*checkpointStream << " ; execute file = " << NxsString::GetEscaped(checkpointTreeFilename) << " ; ";
	kernel.writeCheckpointPostExecSettings(*checkpointStream);
	const std::list< std::pair<MCMCMove *, double> > & movesToWeights = chain.getMoves();
	std::list< std::pair<MCMCMove *, double> >::const_iterator mtwIt = movesToWeights.begin();
	for (; mtwIt != movesToWeights.end(); ++mtwIt) {
		*checkpointStream << "Move weight = " <<   mtwIt->second << ' ';
		mtwIt->first->writePrivateCommands(*checkpointStream);
		*checkpointStream << " ; ";
	}
	*checkpointStream << " mcmc ; \n";
	checkpointStream->flush();
	
}

//==============================================================================
//		csmcmc_chain.cpp
//==============================================================================
void CSMCMCChain::run(unsigned long startIteration, unsigned long nIterations, const ChainSampler & s)
{
	
	if (movesToWeights.empty())
		throw NexusShellException("MCMC cannot be run until move objects have been added");
	/* Set up vectors of the Moves and their probabilities */
	const unsigned nMoveTypes = movesToWeights.size();
	std::vector<double> moveProbs;
	moveProbs.reserve(nMoveTypes);
	std::vector<MCMCMove *> moves;
	moves.reserve(nMoveTypes);
	double total = 0.0;
	std::list< std::pair<MCMCMove *, double> >::const_iterator mtwIt = movesToWeights.begin();
	for (; mtwIt != movesToWeights.end(); ++mtwIt) {
		total += mtwIt->second;
		moveProbs.push_back(mtwIt->second);
		moves.push_back(mtwIt->first);
	}
	for (unsigned i = 0; i < nMoveTypes; ++i) {
		moveProbs[i] /= total;
	}
	
	const unsigned long chkpInterval = s.getCheckpointInterval();
	const unsigned long sampleInterval = s.getSampleInterval();
	const unsigned long reportInterval = s.getReportInterval();

	// calculate initial likelihood
	this->seqLScore();
	unsigned nAccepted = 0;
	unsigned nMovesTried = 0;
	for (unsigned long i = startIteration; i <= nIterations; ++i) {
		if (i % reportInterval == 0)
			s.report(tree, currModel, i, nIterations);
		if (i % sampleInterval == 0)
			s.sample(tree, currModel, i, nIterations);
		if (i % chkpInterval == 0)
			s.checkpoint(tree, currModel, i, nIterations, *this);
		unsigned moveIndex = InferenceKernel::rng.multinomialDraw(&moveProbs[0], nMoveTypes);
		MCMCMove *move = moves[moveIndex];
		MoveHistory mh = move->update();
		nMovesTried += mh.nMovesTried;
		nAccepted += mh.nMovesAccepted;
	}
}

double CSMCMCChain::seqLScore(std::set<const Node*> dirtyNodes) const
{
	NEXUS_SHELL_ASSERT(tree);
	NEXUS_SHELL_ASSERT(currModel);
	std::vector<Node *> preorder =  tree->GetPreorderTraversal();
	std::vector<Node *>::const_reverse_iterator postIt = preorder.rbegin();
	std::vector<Node *>::const_reverse_iterator endIt = preorder.rend();
	const Node * nd = NULL;
	const Node * root = tree->getRoot();
	NEXUS_SHELL_ASSERT(root);
	const bool onlyRecalcDirtyNodes = !(currSeqLScoreDirty || currModel->isDirty());
	FullLAObj * treeLA = tree->getTreeLA();
	
	CalcContext cc;
	const unsigned nCategs = currModel->getNMixtures();
	const unsigned nSubsets = currModel->getNSubsets();
	configure_context(0, nSubsets, nSubsets, 0, nCategs, currModel->getRescaleThreshold(), nCategs, &cc);
	std::vector<CalcContext> ccs;
	ccs.push_back(cc);
	
	for (; postIt != endIt; ++postIt) {
		nd = *postIt;
		NEXUS_SHELL_ASSERT(nd);
		bool needToUpdate = true;
		if (onlyRecalcDirtyNodes) {
			std::set<const Node*>::iterator nIt = dirtyNodes.find(nd);
			if (nIt == dirtyNodes.end())
				needToUpdate = false;
			else {
				dirtyNodes.erase(nIt);
				const Node * p = nd->getParentC();
				if (p)
					dirtyNodes.insert(p);
			}
		}
		if (needToUpdate)
			currModel->updateNode(*nd, ccs);
	}
	
	currSeqLScoreBySite.resize(nSubsets);
	currTotalSeqLScore = currModel->lnLikeFromEffectiveRoot(root,  treeLA, &currSeqLScoreBySite[0]);
	currSeqLScoreDirty = false;
#	if defined (DEBUGGING_LIKELIHOOD) && DEBUGGING_LIKELIHOOD
	if (MCMCKernel::gKernel)
		MCMCKernel::gKernel->debugLikelihoodHook(currTotalSeqLScore, treeLA, this->tree, this->currModel);
#	endif
	return currTotalSeqLScore;
}



void CSMCMCChain::edgeChangedLen(Node *nd, bool /*isRevert*/ )
{
	ndsMoved.insert(nd);
}

CSMCMCChain::CSMCMCChain(double maxD, double timeMultiplier)
	:tree(0L),
	currModel(0L),
	maxDepth(maxD),
	sampleTimeToTreeTimeMultiplier(timeMultiplier)
	,currSeqLScoreDirty(true)
{
	zeroLikeStructFields(&likeStructs);
}

CSMCMCChain::~CSMCMCChain()
{
	if (tree)
		delete tree;
	if (currModel)
		delete currModel;
	freeLikeStructFields(&likeStructs);
}





void MCMCKernel::writeTaxa(std::ostream &o) const
{
	if (this->samples.empty())
		return;
	o << "BEGIN TAXA;\n  Dimensions NTax = " << this->samples.size() << " ;\n  TaxLabels";
	for (std::vector<Sample>::const_iterator s = this->samples.begin(); s != this->samples.end(); ++s) {
		o << ' ' << NxsString::GetEscaped(s->getCodedName());
	}
	o << " ;\nEND;\n";
}

void MCMCKernel::writeTree(std::ostream &o, int treeIndex, int flags, double edgeLenMultiplier) const
{
	const Tree * t = this->getTreeC(treeIndex);	
	writeTree(o, t, flags, edgeLenMultiplier);
}

void MCMCKernel::writeTree(std::ostream &o, const Tree * t, int flags, double edgeLenMultiplier) const
{
	o << "BEGIN Trees;\n";

	o << "Translate\n  ";
	for (unsigned i = 0; i < this->samples.size(); ++i) {
		o << 1 + i << ' ' <<  NxsString::GetEscaped(this->samples[i].getCodedName());
		if (i + 1 < this->samples.size())
			o << ",\n  ";
	}
	o << ";\n";

	
	o << "Tree ChimNeSSweepTree = [&"; 
	o << (flags & NxsFullTreeDescription::NXS_IS_ROOTED_BIT ? 'R' : 'U');
	o << "] ";
	t->writeNewick(o, flags, edgeLenMultiplier);
	o<< " ;\nEND;\n";
}

unsigned MCMCKernel::addSamples(const std::vector<Sample> &samplesArg)
{
	this->clearTrees();
	this->clearTreeDescriptions();
	this->clearData();
	this->samples = samplesArg;
	this->rescaleSampleDepths();

	taxStatus(std::cerr);
	return samples.size();
}

void MCMCKernel::rescaleSampleDepths()
{
	this->maxSampleTime = 0;
	std::vector<Sample>::iterator sampleIt = this->samples.begin();
	for (; sampleIt != this->samples.end(); ++sampleIt)
		this->maxSampleTime = std::max(this->maxSampleTime, sampleIt->getTime());
	for (sampleIt = this->samples.begin(); sampleIt != this->samples.end(); ++sampleIt)
		sampleIt->setDepth(this->maxSampleTime, this->sampleTimeToTreeTimeMultiplier);
}



unsigned MCMCKernel::getNPMats() const {
	const NxsCDiscreteMatrix & mat = charMatrix.getConstNativeC();
	const unsigned nLeaves = mat.nTax;
	return 2*nLeaves;
}

LikeStructsBundle MCMCKernel::generateNewLikeStructs( unsigned nCLAs) {
	const NxsCDiscreteMatrix & mat = charMatrix.getConstNativeC();
	const unsigned nLeaves = mat.nTax;
	if (nCLAs == 0)
		nCLAs = 2*nLeaves;
	const unsigned nRates = nRateCats;
	const unsigned nPMats = this->getNPMats();
	return newLikeStructs(mat, nCLAs, nPMats,  nRates);
}
void MCMCKernel::addLikelihoodStructsToChain(CSMCMCChain & chain, unsigned nCLAs)
{
	std::cerr << "Constructing data structures for likelihood calculations...\n";

	if (!this->charMatrixInitialized)
		this->initCharMatrix();

	const NxsCDiscreteMatrix & mat = charMatrix.getConstNativeC();
	const unsigned nLeaves = mat.nTax;
	if (nCLAs == 0)
		nCLAs = 2*nLeaves;
	

	chain.likeStructs = this->generateNewLikeStructs(nCLAs);
	if (chain.likeStructs.leafData == 0L) {
		throw NexusShellException("Insufficient memory");
	}
	if (chain.tree == NULL) {
		throw NexusShellException("Chain does not have a tree.");
	}

	const unsigned nPMats = this->getNPMats();
	unsigned claCounter = 0;
	unsigned pmatCounter = 0;
	std::vector<Node *> preorder =  chain.tree->GetPreorderTraversal();
	for (std::vector<Node *>::iterator ndIt = preorder.begin(); ndIt != preorder.end(); ++ndIt) {
		Node * nd = *ndIt;
		if (nd->isLeaf()) {
			unsigned leafInd = nd->getIndex();
			NEXUS_SHELL_ASSERT(leafInd < nLeaves);
			nd->setLeafData(chain.likeStructs.leafData[leafInd]);
		}
		else {
			if (claCounter < nCLAs)
				nd->setCLA(chain.likeStructs.clas[claCounter++]);
			else {
				freeLikeStructFields(&(chain.likeStructs));
				zeroLikeStructFields(&(chain.likeStructs));
				addLikelihoodStructsToChain(chain, 2*nCLAs);
				return;
			}
		}
		if (nd->getRSibC() || nd->getLSibC()) {
			// the root and nodes of degree 2 are the only nodes without PMats
			NEXUS_SHELL_ASSERT(pmatCounter < nPMats);
			nd->setPMatArray(chain.likeStructs.pmats[pmatCounter++]);
		}
	}
	chain.tree->setTreeLikeArray(chain.likeStructs.treeLike);
	chain.spareCLAs.clear();
	for (; claCounter < nCLAs; ++claCounter) {
		chain.spareCLAs.push_back(chain.likeStructs.clas[claCounter]);
	}


	SeqModel * m = this->generateNewSeqModel(chain.likeStructs);
	chain.setModel(m);
	m->addTreeLikePtr(chain.likeStructs.treeLike);
}

MCMCKernel::MCMCKernel()
	:nRateCats(4), // /* Assuming 4 cat gamma distributed rates here!  */ 
	maxDepth(1.0e12),
	sampleTimeToTreeTimeMultiplier(0.010),
	maxSampleTime(0),
	stateFreqStartingValues(4,0.25),
	rMatStartingValues(6,1.0),
	gammaStartingValue(0.5),
	pinvStartingValue(0.0),
	mcmcIterations(2L),
	currMcmcIteration(1L),
	sampleInterval(1L),
	reportInterval(1L),
	checkpointInterval(1L),
	charMatrixInitialized(false),
	geneticCode(NXS_GCODE_NO_CODE)
	{
	vector<int> row(4,-1);
	row[1] = 0;
	row[2] = 1;
	row[3] = 2;
	paramIdentityMatrix.push_back(row);
	row[0] = 0;
	row[1] = -1;
	row[2] = 3;
	row[3] = 4;
	paramIdentityMatrix.push_back(row);
	row[0] = 1;
	row[1] = 3;
	row[2] = -1;
	row[3] = 5;
	paramIdentityMatrix.push_back(row);
	row[0] = 2;
	row[1] = 4;
	row[2] = 5;
	row[3] = -1;
	paramIdentityMatrix.push_back(row);
	}

void MCMCKernel::addStateFreqMove(int treeIndex, double wt, int subsetInd, const std::vector<double> & windows, const std::vector<double> & priors) 
{
	CSMCMCChain & chain = this->getChainForTreeIndex(treeIndex);
	if (chain.getModel() == 0L)
		this->addLikelihoodStructsToChain(chain);
	
	SeqModel * seqM = chain.getModel();
	
	StateFreqParamGroup & sfpg = seqM->getStateFreqRef(subsetInd);
	SampledSurface * surf = new ProbSeqGivenGenealogy(chain);
	MCMCMove * m = new DiscreteProbMove(sfpg, windows, priors, surf);
	m->setModelIndex(subsetInd);
	movesToDelete.push_back(m);
	surfacesToDelete.push_back(surf);
	chain.addMove(m, wt);
}

void MCMCKernel::addRateMatMove(int treeIndex, double wt, int subsetInd, const std::vector<double> & windows, const std::vector<double> & priors) 
{
	CSMCMCChain & chain = this->getChainForTreeIndex(treeIndex);
	if (chain.getModel() == 0L)
		this->addLikelihoodStructsToChain(chain);
	
	SeqModel * seqM = chain.getModel();
	
	ParamGroup & pg = seqM->getRMatParamGroup(subsetInd);
	SampledSurface * surf = new ProbSeqGivenGenealogy(chain);
	MCMCMove * m = new PositiveParamGroupMove(pg, windows, priors, surf);
	movesToDelete.push_back(m);
	surfacesToDelete.push_back(surf);
	chain.addMove(m, wt);
}


void MCMCKernel::writeCheckpointPreExecSettings(std::ostream & out) const 
{
	if (geneticCode != NXS_GCODE_NO_CODE) {
		out << " Set geneticcode = " << geneticCodeEnumToName(geneticCode) << " ; "; 
	}
}

void MCMCKernel::writeMCMCCheckpointInfo(std::ostream & out) const
{
	out << "Set sampleinterval = "<< sampleInterval;
	out << " reportInterval = " <<  reportInterval;
	out << " checkpointInterval = " <<  checkpointInterval;
	if (!checkpointFilename.empty())
		out << " checkpointfile = " << NxsString::GetEscaped(checkpointFilename);
	if (!sampleTreeFilename.empty())
		out << " sampletreefile = " << NxsString::GetEscaped(sampleTreeFilename);
	if (!sampleParamFilename.empty())
		out << " sampleparamfile = " << NxsString::GetEscaped(sampleParamFilename);
	if (!reportFilename.empty())
		out << " reportfile = " << NxsString::GetEscaped(reportFilename);
}

void MCMCKernel::runMCMC(int treeIndex)
{
	CSMCMCChain & chain = this->getChainForTreeIndex(treeIndex);
	if (chain.getModel() == 0L)
		this->addLikelihoodStructsToChain(chain);
	
	ChainSampler s(*this);

	std::ofstream checkpointFile;
	std::ofstream reportFile;
	std::ofstream paramFile;
	std::ofstream treeFile;
	
	if (checkpointFilename.empty())
		s.setCheckpoint(ULONG_MAX, NULL);
	else if (checkpointFilename == "cout")
		s.setCheckpoint(checkpointInterval, &std::cout);
	else if (checkpointFilename == "cerr")
		s.setCheckpoint(checkpointInterval, &std::cerr);
	else {
		checkpointFile.open(checkpointFilename.c_str());
		s.setCheckpoint(checkpointInterval, &checkpointFile);
	}
	
	
	if (reportFilename.empty())
		s.setReport(ULONG_MAX, NULL);
	else if (reportFilename == "cout")
		s.setReport(reportInterval, &std::cout);
	else if (reportFilename == "cerr")
		s.setReport(reportInterval, &std::cerr);
	else {
		reportFile.open(reportFilename.c_str());
		s.setReport(reportInterval, &reportFile);
	}
	
	std::ostream * tp;
	std::ostream * pp;
	if (sampleTreeFilename.empty())
		tp = NULL;
	else if (sampleTreeFilename == "cout")
		tp = &std::cout;
	else if (sampleTreeFilename == "cerr")
		tp = &std::cerr;
	else {
		treeFile.open(sampleTreeFilename.c_str());
		tp = &treeFile;
	}
	if (sampleParamFilename.empty())
		pp = NULL;
	else if (sampleParamFilename == "cout")
		pp = &std::cout;
	else if (sampleParamFilename == "cerr")
		pp = &std::cerr;
	else {
		paramFile.open(sampleParamFilename.c_str());
		pp = &paramFile;
	}
	s.setSample(sampleInterval, tp, pp);
	
	return chain.run(currMcmcIteration, mcmcIterations, s);
}

double MCMCKernel::seqLScore(int treeIndex)
{
	CSMCMCChain & chain = this->getChainForTreeIndex(treeIndex);
	if (chain.getModel() == 0L)
		this->addLikelihoodStructsToChain(chain);
	return chain.seqLScore();
}

const Tree * MCMCKernel::getTreeC(int treeIndex) const 
{
	MCMCKernel * p = const_cast<MCMCKernel * >(this);
	Tree * t = p->getTree(treeIndex);
	return t;
	
}


CSMCMCChain & MCMCKernel::getChainForTreeIndex(int treeIndex)
{
	const int  nTinMem = (int)getNumTrees();
	const int  ntd = (int)treeDescriptions.size();
	int ind = (treeIndex < 0 ? nTinMem + ntd + treeIndex :  treeIndex);
	if (ind < 0  || ind >= nTinMem + ntd) {
		if (ntd == 0)
			throw NexusShellException("No trees in memory");
		throw NexusShellException("Tree index out of range");
	}
	if (ind < ntd) {
		getTree(treeIndex);
		return *chains.rbegin();
	}
	return chains[ind - ntd];
}

Tree * MCMCKernel::getTree(int treeIndex)
{
	const int  nTinMem = (int)getNumTrees();
	const int  ntd = (int)treeDescriptions.size();
	int ind = (treeIndex < 0 ? nTinMem + ntd + treeIndex :  treeIndex);
	if (ind < 0  || ind >= nTinMem + ntd) {
		if (ntd == 0)
			throw NexusShellException("No trees in memory");
		throw NexusShellException("Tree index out of range");
	}
	if (ind < ntd) {
		const NxsFullTreeDescription & td = treeDescriptions[ind];
		Tree * t = new Tree(samples, td);
		chains.push_back(CSMCMCChain(maxDepth, sampleTimeToTreeTimeMultiplier));
		CSMCMCChain & c = *chains.rbegin();
		c.setTree(t);
		return t;
	}
	return chains[ind - ntd].getTree();
}


void MCMCKernel::taxStatus(std::ostream &o) const
{
	o << samples.size() << " taxa in memory\n";
	for (unsigned i = 0 ; i < samples.size(); ++i) {
		const Sample & sample = samples[i];
		o << i << ' ' << sample << '\n';
	}
}




std::ostream & Sample::writeCodedName(std::ostream &o) const
{
	o << this->sampleLabel;
	if (this->intTime != INT_MAX)
		o << '.' << this->strTime << "." << this->id << '.' << this->accession;
	return o;
}

using namespace chim;
#if defined (DEBUGGING_LIKELIHOOD) && DEBUGGING_LIKELIHOOD

//==============================================================================
//		chimne_sweep_kernel.cpp
//==============================================================================
void chim::ChimneSweepKernel::debugLikelihoodHook(double likeCalc, const FullLAObj *, Tree *tree, SeqModel *model)
{
	ChimneSweepSeqModel * currModel = (ChimneSweepSeqModel *)(model);
	std::ofstream likeTree("likelihoodCheck.tre");
	likeTree << "#NEXUS\n";
	MCMCKernel::gKernel->writeTaxa(likeTree);
	MCMCKernel::gKernel->writeTree(likeTree, tree, -1);
	likeTree.close();

	ofstream paupCmd("paupCmd.nex");
	paupCmd << "#NEXUS\nBegin Paup;\n\tSet storebrlens;\ngettrees file = likelihoodCheck.tre;\n";
	paupCmd << "	LScore  / UserBrlens " << currModel->getPaupLSetOptions() << ";END;\n";
	paupCmd << "[!\n****************************************************************************************\n The PAUP LSCORE SHOULD HAVE BEEN\n-ln L   " ;
	paupCmd.precision(5);
	paupCmd.setf(ios::fixed,ios::floatfield);  
	paupCmd << -likeCalc << "\n**************************************************************************************** ]\n";
	paupCmd.close();
	
	if (system("./checkWithPaup.py") != 0) {
		throw  NexusShellException("likelihood bug: disagreement with PAUP");
	}
}
#endif

std::string ChimneSweepSeqModel::getPaupLSetOptions() const
{
	const unsigned nStates = cCoreModel.dim;
	const std::vector<double> & freqs =  stateFreq.getFreqs();
	ostringstream s;
	s.precision(10);
	s.setf(ios::fixed,ios::floatfield);  

	s << " BaseFreq = (";
	for (unsigned i = 0; i < nStates-1; ++i)
		s << ' ' << freqs[i] << ' ';
		
	
	s << ") nst = 6 rmat =(";
	const double refRate = qRateMatrix[2][3]->getValue();
	for (unsigned i = 0; i < nStates-2; ++i) {
		const QRateRow & qMatParamRow = qRateMatrix[i];
		for (unsigned j = i+1; j < nStates; ++j) {
			const Parameter * p = qMatParamRow[j];
			assert(p);
			s << ' ' << p->getValue()/refRate << ' ';
		}
	}
	s << ") ";
	const double gammaShapeV = gammaShape.getValue();
	if (gammaShapeV > 0.0)
		s << " rates=gamma shape=" << gammaShapeV;
	else
		s << " rates=equal";
	const double pinv = pInvar.getValue();
	if (pinv > 0.0)
		s << " pinv=" << gammaShapeV;
	else
		s << " pinv=0.0";
	return s.str();
}

void ChimneSweepSeqModel::writePrivateCommands(std::ostream & s) const
{
	const unsigned nStates = cCoreModel.dim;
	s.precision(10);
	s.setf(ios::fixed,ios::floatfield);  
	const std::vector<Parameter> &freq = stateFreq.getRawParam();
	s << "Set  freq = (";
	for (unsigned i = 0; i < nStates; ++i)
		s << ' ' << NxsString::GetEscapedDouble(freq[i].getValue()) << ' ';
		
	
	s << ") rmat =(";
	for (unsigned i = 0; i < nStates - 1 ; ++i) {
		const QRateRow & qMatParamRow = qRateMatrix[i];
		for (unsigned j = i+1; j < nStates; ++j) {
			const Parameter * p = qMatParamRow[j];
			assert(p);
			s << ' ' << NxsString::GetEscapedDouble(p->getValue()) << ' ';
		}
	}
	s << ") ";
	const double gammaShapeV = gammaShape.getValue();
	if (gammaShapeV > 0.0)
		s << " shape =" << NxsString::GetEscapedDouble(gammaShapeV);
	s << ';' ;
}


void chim::ChimneSweepKernel::setEdgeLenFromDepth(int treeIndex)
{
	Tree * t = this->getTree(treeIndex);
	t->setEdgeLenFromDepth();
}

void chim::ChimneSweepKernel::addLikelihoodStructsToChain(CSMCMCChain & chain, unsigned nCLAs)
{
	MCMCKernel::addLikelihoodStructsToChain(chain, nCLAs);

	ChimneSweepSeqModel * m = static_cast<ChimneSweepSeqModel *>(chain.getModel());
	std::vector<Parameter> & sfp = m->getStateFreqRef(0).getRawParamRef();
	NEXUS_SHELL_ASSERT(sfp.size() == stateFreqStartingValues.size());
	std::vector<Parameter>::iterator pIt = sfp.begin();
	std::vector<double>::iterator vIt = stateFreqStartingValues.begin();
	for (; pIt != sfp.end(); ++pIt, ++vIt) {
		pIt->setValue(*vIt);
	}
	
	
	std::vector<Parameter> &rp = m->getRMatParamRef(0);
	NEXUS_SHELL_ASSERT(rp.size() == rMatStartingValues.size());
	pIt = rp.begin();
	vIt = rMatStartingValues.begin();
	for (; pIt != rp.end(); ++pIt, ++vIt) {
		pIt->setValue(*vIt);
	}
	
	m->getGammaShapeParam(0).setValue(gammaStartingValue);
	m->getPInvarParam(0).setValue(pinvStartingValue);
}

SeqModel * chim::ChimneSweepKernel::generateNewSeqModel(LikeStructsBundle & bundle)
{
	NEXUS_SHELL_ASSERT(bundle.model);
	NEXUS_SHELL_ASSERT(bundle.model[0]);
	return new ChimneSweepSeqModel(*(bundle.model[0]), bundle.asrv, paramIdentityMatrix);
}

void chim::ChimneSweepKernel::ultrametricize(int treeIndex, double minEdgeLen)
{
	Tree * t = this->getTree(treeIndex);
	t->ultrametricize(minEdgeLen);	
}


ChimneSweepSeqModel::ChimneSweepSeqModel(
	DSCTModelObj & model, 
	ASRVObj * asrvObj, 
	const std::vector< std::vector<int> > & paramIdentityMatrix)
:cCoreModel(model),
asrv(asrvObj),
rMatParams(0),
stateFreq(model.dim)
{
	std::vector<Parameter> & sf = stateFreq.getRawParamRef();
	for (std::vector<Parameter>::iterator pIt = sf.begin(); pIt != sf.end(); ++pIt)
		pIt->addListener(this);
	int maxInd = -1;
	std::vector< std::vector<int> >::const_iterator rIt = paramIdentityMatrix.begin();
	for (; rIt != paramIdentityMatrix.end(); ++rIt) {
		std::vector<int>::const_iterator cIt = rIt->begin(); 
		for (; cIt != rIt->end(); ++cIt) {
			maxInd = max(maxInd, *cIt);
		}
	}
	if (maxInd >= 0)
		rMatParams.resize(maxInd + 1);
	std::vector<Parameter> & rateMatParams =  rMatParams.getRawParamRef();
	rIt = paramIdentityMatrix.begin();
	for (; rIt != paramIdentityMatrix.end(); ++rIt) {
		QRateRow qr;
		std::vector<int>::const_iterator cIt = rIt->begin(); 
		for (; cIt != rIt->end(); ++cIt) {
			Parameter * p = (*cIt < 0 ? NULL : &(rateMatParams[*cIt]));
			qr.push_back(p);
		}
		NEXUS_SHELL_ASSERT(qr.size() == model.dim);
		qRateMatrix.push_back(qr);
	}
	NEXUS_SHELL_ASSERT(qRateMatrix.size() == model.dim);
	for (unsigned i = 0; i < model.dim; ++i) {
		if (qRateMatrix[i][i])
			throw  ChimneSweepException("Diagonal element of qMat must be 0");
		for (unsigned j = i; j < model.dim; ++j) {
			if (qRateMatrix[i][j] != qRateMatrix[j][i])
				throw  ChimneSweepException("Q-matrix is not reversible");
		}
	}
	std::vector<Parameter>::iterator pIt = rateMatParams.begin();
	for (; pIt != rateMatParams.end(); ++pIt)
		pIt->addListener(this);
	std::vector<Parameter> &sfp =  stateFreq.getRawParamRef();
	for (pIt = sfp.begin(); pIt != sfp.end(); ++pIt)
		pIt->addListener(this);
}


// use prameters to set the cCoreModel.qmat
void ChimneSweepSeqModel::propagateQMat() const {
	NEXUS_SHELL_ASSERT(cCoreModel.q_mat);
	NEXUS_SHELL_ASSERT(cCoreModel.q_mat[0]);
	NEXUS_SHELL_ASSERT(qRateMatrix.size() == cCoreModel.dim);
	const unsigned nStates = cCoreModel.dim;
	const std::vector<double> & freqs =  stateFreq.getFreqs();
	double wtdDiagSum = 0.0;
	for (unsigned i = 0; i < nStates; ++i) {
		const QRateRow & qMatParamRow = qRateMatrix[i];
		double * qMatRow = cCoreModel.q_mat[i];
		qMatRow[i] = 0;
		for (unsigned j = 0; j < nStates; ++j) {
			if (i != j) {
				const Parameter * p = qMatParamRow[j];
				assert(p);
				const double unNorm = freqs[j]*p->getValue();
				qMatRow[j] = unNorm;
				qMatRow[i] -= unNorm;
			}
		}
		wtdDiagSum -= qMatRow[i]*freqs[i];
	}
	const double mult = 1.0/wtdDiagSum;
	for (unsigned i = 0; i < nStates; ++i) {
		double * qMatRow = cCoreModel.q_mat[i];
		for (unsigned j = 0; j < nStates; ++j) {
			qMatRow[j] *= mult;
		}
	}
	cCoreModel.eigen_calc_dirty = 1;
	
	
	for (std::vector<FullLAObj *>::const_iterator tlpIt = treeLikePtrs.begin(); tlpIt != treeLikePtrs.end(); ++tlpIt) {
		FullLAObj * tl = *tlpIt;
		const unsigned nc = tl->cso.n_categ_or_subs;
		for (unsigned i = 0; i < nc; ++i) {
			double * scf = tl->state_categ_freqs[i];
			for (unsigned j = 0; j < nStates; ++j)
				scf[j] = freqs[j];
		}
	}
}

void ChimneSweepSeqModel::propagateRateHet() const {
	if (!asrv)
		return;
	internal_asrv_set_shape(asrv, gammaShape.getValue());
}

PMatArrayObj * ChimneSweepSeqModel::getNodeUpdatedPMat(const Node & nd, double edgeLen) {
	PMatArrayObj * p = nd.getPMatArray();
	NEXUS_SHELL_ASSERT(p);
	if (this->isDirty()) {
		this->propagateRateHet();
		this->propagateQMat();
	}
	setIsDirty(false);
	
	const int rc = do_asrv_pmat_array_calc(p, &cCoreModel, asrv, edgeLen);
	NEXUS_SHELL_ASSERT(rc != 0);
	return p;
}

void chim::ChimneSweepKernel::addEdgeLenMove(int , double , int , double , double ) {
	throw ChimneSweepException("Edge length moves are not supported in chimne sweep -- use node depth moves");
}

void  chim::ChimneSweepKernel::addPartitionedFreqMove(int , double , const std::list< std::vector<double> > & , const std::list< std::vector<double> > & ) {
	throw ChimneSweepException("Partitioned Freq moves are not supported in chimne sweep -- partitioned models are not available.");
}

void chim::ChimneSweepKernel::writeCheckpointPostExecSettings(std::ostream & out) const 
{
	out << "; Ultrametricize ; EdgeLenFromDepth ; ";
}

