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

#ifndef GEOGRAM_SYMBOLICS_TINY_CAS
#define GEOGRAM_SYMBOLICS_TINY_CAS

#include <OGF/Experiment/common/common.h>
#include <geogram/basic/memory.h>
#include <string>
#include <iostream>

namespace OGF {

    namespace TinyCAS {

	/**
	 * \brief Quick and dirty CAS system for polynoms
	 */
	class Polynom {
	public:

	    struct Monom {

	        Monom(index_t nb_var) : a(0.0), pow(nb_var, 0) {
		}

	        Monom(index_t nb_var, index_t var) : a(1.0), pow(nb_var,0) {
		    pow[var] = 1;
		}

	        Monom(index_t nb_var, double a_in) : a(a_in), pow(nb_var,0) {
		}

		index_t nb_var() const {
		    return pow.size();
		}

		Monom operator*(const Monom& rhs) const {
		    geo_assert(rhs.nb_var() == nb_var());
		    Monom result(nb_var());
		    result.a = a*rhs.a;
		    for(index_t i=0; i<nb_var(); ++i) {
			result.pow[i] = pow[i]+rhs.pow[i];
		    }
		    return result;
		}

		std::string to_string(
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
		double a;
		vector<index_t> pow;
	    };

	    index_t nb_var() const {
		return monoms_[0].nb_var();
	    }

	    Polynom& operator+=(const Monom& M1) {
		geo_assert(M1.nb_var() == nb_var());
		if(M1.a == 0.0) {
		    return *this;
		}
		for(Monom& M2: monoms_) {
		    if(M2.pow == M1.pow) {
			M2.a += M1.a;
			return *this;
		    }
		}
		monoms_.push_back(M1);
		return *this;
	    }

	    Polynom& operator-=(const Monom& M1) {
		geo_assert(M1.nb_var() == nb_var());
		Monom mM1 = M1;
		mM1.a = -mM1.a;
		return this->operator+=(mM1);
	    }

	    Polynom& operator+=(const Polynom& rhs) {
		geo_assert(rhs.nb_var() == nb_var());
		for(const Monom& M2: rhs.monoms_) {
		    *this += M2;
		}
		return *this;
	    }

	    Polynom operator+(const Polynom& rhs) const {
		Polynom result = *this;
		result += rhs;
		return result;
	    }

	    Polynom operator-=(const Polynom& rhs) {
		geo_assert(rhs.nb_var() == nb_var());
		for(const Monom& M2: rhs.monoms_) {
		    *this -= M2;
		}
		return *this;
	    }

	    Polynom operator-(const Polynom& rhs) const {
		Polynom result = *this;
		result -= rhs;
		return result;
	    }

	    Polynom operator*(const Polynom& rhs) const {
		geo_assert(rhs.nb_var() == nb_var());
		Polynom result(monoms_[0].pow.size());
		for(const Monom& M1: monoms_) {
		    for(const Monom& M2: rhs.monoms_) {
			result += (M1*M2);
		    }
		}
		return result;
	    }

	    Polynom pow(index_t p) {
		if(p == 0) {
		    return Polynom(nb_var(),1.0);
		}
		Polynom result = *this;
		while(p > 1) {
		    if((p & 1) == 0) {
			result = result * result;
			p /= 2;
		    } else {
			result = result * (*this);
			--p;
		    }
		}
		return result;
	    }

	    std::string to_string(
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

	    Polynom(index_t nb_var, double val=0.0) {
		monoms_.push_back(Monom(nb_var, val));
	    }

	    Polynom(index_t nb_var, index_t var) {
		monoms_.push_back(Monom(nb_var, var));
	    }

	    bool is_zero() const {
		for(const Monom& M: monoms_) {
		    if(M.a != 0.0) {
			return false;
		    }
		}
		return true;
	    }

	public:
	    vector<Monom> monoms_;
	};
    }
}

#endif
