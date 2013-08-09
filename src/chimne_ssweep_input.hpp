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
#if ! defined (CHIMNE_SWEEP_INPUT_HPP)
#define CHIMNE_SWEEP_INPUT_HPP

#include <iostream>

class NxsSimpleCommandStrings;

#include "nexus_shell.hpp"


namespace chim {

class ChimneSweepKernel;
bool doChimneSweepNexusCommand(const NxsSimpleCommandStrings &s, InferenceKernelBundle & bundle);


class ChimneSweepException: public NexusShellException
	{
	public:
		virtual ~ChimneSweepException() throw()
			{
			}

		ChimneSweepException(const std::string & s):NexusShellException(s){}
		ChimneSweepException(const char * s):NexusShellException(s){}

		const char * what () const throw ()
			{
			return msg.empty() ? "Unknown ChimneSSweep Exception" : msg.c_str();
			}
	};


} // namespace chim 


#endif
