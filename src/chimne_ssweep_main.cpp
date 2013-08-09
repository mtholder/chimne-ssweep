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

#include <iostream>

#include "chimne_ssweep_input.hpp"
#include "chimne_ssweep_kernel.hpp"


void printHelp(std::ostream & out);

void printHelp(std::ostream & out) 
{
	out << "chimne_ssweep reads a NEXUS file as its only command line argument.\n";
	out << "The correct invocation is:\n\n    chimne_ssweep <path to NEXUS file>\n";
}

int main(int argc, char *argv[])
{
	if (argc < 2) {
		printHelp(std::cerr);
		return 1;
	}
	chim::ChimneSweepKernel kernel;
#	if defined (DEBUGGING_LIKELIHOOD) && DEBUGGING_LIKELIHOOD
		MCMCKernel::gKernel = &kernel;
#	endif

	Tree::defaultTreeWritingFlags = (     NxsFullTreeDescription::NXS_EDGE_LENGTH_UNION 
										| NxsFullTreeDescription::NXS_HAS_DEG_TWO_NODES_BIT
										| NxsFullTreeDescription::NXS_HAS_POLYTOMY_BIT
										| NxsFullTreeDescription::NXS_IS_ROOTED_BIT
										| NxsFullTreeDescription::NXS_HAS_INTERNAL_NAMES_BIT
									);


	for (int i = 1; i < argc; ++i) {
		const char * curr_arg = argv[i];
		if (strlen(curr_arg) > 1 && curr_arg[0] == '-' && curr_arg[1] == 'h') {
			printHelp(std::cout);
			return 0;
		}
	}
	InferenceKernelBundle kernelWrapper(kernel, "chimne_ssweep", "CHIMNESSWEEP", std::cout, std::cerr, chim::doChimneSweepNexusCommand);
	for (int i = 1; i < argc; ++i) {
		readFilepathOrExit(argv[i], kernelWrapper);
	}
	return 0;
}

