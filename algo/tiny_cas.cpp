/*
 *  OGF/Graphite: Geometry and Graphics Programming Library + Utilities
 *  Copyright (C) 2000-2015 INRIA - Project ALICE
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact for Graphite: Bruno Levy - Bruno.Levy@inria.fr
 *  Contact for this Plugin: Bruno Levy - Bruno.Levy@inria.fr
 *
 */

#include <OGF/Experiment/algo/tiny_cas.h>

namespace OGF {
    namespace TinyCAS {

        /**********************************************************************/

	std::string Polynom::Monom::to_string(
	    const std::vector<std::string>& var_names
	) const {
	    geo_assert(index_t(var_names.size()) == nb_var());
	    std::string result;
	    if(a == 0.0) {
		return "0.0";
	    }
	    if(a == 1.0) {
		result += "+";
	    } else if(a == -1.0) {
		result += "-";
	    } else {
		if(a > 0.0) {
		    result += "+";
		}
		result += String::to_string(a);
	    }
	    bool first = true;
	    for(index_t i=0; i<nb_var(); ++i) {
		if(pow[i] != 0) {
		    if(!first) {
			result += "*";
		    }
		    result += var_names[i];
		    if(pow[i] != 1) {
			result += "^";
			result += String::to_string(pow[i]);
		    }
		    first = false;
		}
	    }
	    if(first) {
		result += "1";
	    }
	    return result;
	}

	std::string Polynom::to_string(
	    const std::vector<std::string>& var_names
	) const {
	    if(is_zero()) {
		return "0";
	    }
	    std::string result;
	    for(const Monom& M: monoms_) {
	    if(M.a == 0.0) {
		continue;
	    }
	    if(result.length() != 0) {
		result += " ";
	    }
	    result += M.to_string(var_names);
	    }
	    return result;
	}

       /***********************************************************************/

    } // namespace TinyCAS
} // namespace OGF
