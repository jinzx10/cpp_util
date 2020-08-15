#ifndef __CXX_UTIL_PHYSICS_H__
#define __CXX_UTIL_PHYSICS_H__

#include <armadillo>
#include <tuple>
#include "math_helper.h"

namespace cxut {

	// sigma_x and sigma_p of the Wigner quasiprobability distribution of harmonic oscillators
	inline std::tuple<double, double> ho_wigner(double const& mass, double const& omega, double const& kT = 0) {
		double sigma_x, sigma_p;
		if (kT < arma::datum::eps) {
			sigma_x = std::sqrt(0.5 / mass / omega);
			sigma_p = std::sqrt(0.5 * omega * mass);
		} else {
			sigma_x = std::sqrt( 0.5 / mass / omega / std::tanh(omega/2.0/kT) );
			sigma_p = std::sqrt( 0.5 * mass * omega / std::tanh(omega/2.0/kT) );
		}
		return std::make_tuple(sigma_x, sigma_p);
	}
	

	// Boltzmann weight
	inline arma::vec boltzmann(arma::vec const& E, double const& kT) {
		arma::uword imin = E.index_min();
		arma::vec v;
		if ( std::abs(kT) < arma::datum::eps ) {
			v.zeros(E.n_elem);
			v(imin) = 1.0;
			return v;
		}
		v = arma::exp(-(E-E(imin))/kT);
		return v / accu(v);
	}
	

	// Fermi function
	inline double fermi(double const& E, double const& mu, double const& kT) {
		return ( std::abs(kT) < arma::datum::eps ) ? 
			(E <= mu) : 1.0 / ( std::exp( (E - mu) / kT ) + 1.0 );
	}
	
	inline arma::vec fermi(arma::vec const& E, double const& mu, double const& kT) {
		return ( std::abs(kT) < arma::datum::eps ) ? 
			arma::conv_to<arma::vec>::from(E <= mu) : 1.0 / ( exp( (E - mu) / kT ) + 1.0 );
	}
	

	// solve for the chemical potential given E, N and kT
	inline int findmu(double& mu, arma::vec const& E, arma::uword const& N, double const& kT = 0.0) {
		if ( N > E.n_elem ) {
			std::cerr << "findmu: there are more particles than energy levels." << std::endl;
			return -1;
		}
	
		if ( std::abs(kT) < arma::datum::eps ) {
			mu = arma::sort(E, "ascend").eval()(N-1);
			return 0;
		}
	
		auto dn = [&] (double const& mu) { return arma::accu(fermi(E, mu, kT)) - N; };
		mu = E(0);
		return broydenroot(dn, mu, 1e-12);
	}
}

#endif
