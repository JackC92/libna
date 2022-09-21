#ifndef NA_LEGENDRE_H
#define NA_LEGENDRE_H
#include "Eigen/Core"

// This file contains code for constructing Guass-Legendre quadrature rules, and for
// forming and manipulating univariate Legendre expansions.  A univariate Legendre 
// expansion is a sum of the form 
//
//             n-1
//     f(x) = \sum   a_i \tilde{P}_i(x),                                                   (1)
//             i=0
//
// where \tilde{P}_i denotes the L^2 normalized Legendre polynomial of the first
// kind of degree i.  This code provides routines for real-valued expansions.
//
// Expansions of the form (1) are represented either via the vector of coefficients
//  
//   ( a_0     )
//   ( a_1     )
//   ( ...     )                                                                            (2)
//   ( a_{n-1} )
//
// or via the vector 
//
//   ( f(x_1) \sqrt{w_1} )
//   ( f(x_2) \sqrt{w_2} )
//   (      ...          )                                                                 (3)
//   ( f(x_n) \sqrt{w_n} )
//
// of their *scaled* values at the nodes of the n-point Gauss-Legendre quadrature 
// rule.  Note that the matrix taking (2) to (3) is orthogonal.
//
// The following routines should be regarded as public:
//     
//   quadrature - return the nodes and weights of the n-point Gauss-Legendre 
//     quadrature rule
//
//   expansion_coefficients - compute the vector (2) of the coefficients of one of more 
//      expansions of the form (1) given the vectors of their scaled values (3)
//
//   interpolate - use barycentric interpolation to evaluate an expansion 
//    of the form (1) given its *scaled* values at nodes of the n-point
//    Gauss-Legendre quadrature rule
//
//   evaluate_expansion - evaluate one or more expansions of the form (1) (and their
//     derivatives) at a specified point given the vector (2) of coefficients
//
//   interpolation_matrix - return an (m,n) matrix which takes the vector (3) of the
//     *scaled* values of an expansion at the nodes of the n-point Gauss-Legendre rule
//     to *scaled* values at a user-specified collection of points
//
//   coefficient_matrix - return the (n,n) matrix which takes the vector (3) of 
//     the *scaled* values of an expansion of the form (1) to the vector (2) of its
//     coefficients -- note that this matrix is orthogonal and its transpose
//     takes the vector (2) to the vector (3)

namespace na
{
	namespace legendre
	{
		// Compute the coefficients of the Legendre polynomial of degree n
		// in the monomial basis
		//
		// Input parameters:
		//   n - the degree of the Legendre polynomial
		// 
		// Output parameters:
		//   coefs - the coefficients of the Legendre polynomial
		Eigen::VectorXd polynomial_coefficients(const int n);

		// Evaluate the Legendre polynomial of degree n at the point x
		//
		// Input parameters:
		//   n - the degree of the Legendre polynomial
		//   x - the point at which to evaluate the Legendre polynomial
		// 
		// Output parameters:
		//   val - the value of the Legendre polynomial at the point x
		double evaluate_polynomial(
			const int n,
			const double x);

		// Evaluate the derivative of the Legendre polynomial of degree
		// n at the point x
		//
		// Input parameters:
		//   n - the degree of the Legendre polynomial
		//   x - the point at which to evaluate the Legendre polynomial
		// 
		// Output parameters:
		//   val - the value of the Legendre polynomial at the point x
		double evaluate_derivative(
			const int n,
			const double x);

		// Evaluate the Legendre polynomial of degree n and its derivative
		// at the point x
		//
		// Input parameters:
		//   n - the degree of the Legendre polynomial
		//   x - the point at which to evaluate the Legendre polynomial
		// 
		// Output parameters:
		//   pol - the value of the Legendre polynomial at the point x
		//   der - the value of the derivative of the Legendre polynomial at the point x
		void evaluate_polynomial(
			const int n,
			const double x,
			double& pol,
			double& der);

		// Form the n-point Gauss-Legendre quadrature rule.
		//
		// Input parameters:
		//   n - the length of the desired quadrature rule
		//
		// Output parameters:
		//   xslege - the nodes of the desired quadrature rule
		//   whtslege - the weights of the desired quadrature rule
		void quadrature(
			const int n,
			Eigen::VectorXd& xslege,
			Eigen::VectorXd& whtslege);

		// Compute the vector (2) of the coefficients of one of more 
		// expansions of the form (1) given the vectors of their scaled values (3)
		//
		// Input parameters:
		//   n - the number of terms in the expansion
		//   vals - the *scaled* values of (1) at the n Gauss-Legendre nodes
		//
		// Output parameters:
		//   coefs - an array specifying the coefficients
		Eigen::VectorXd expansion_coefficients(
			const int n,
			const Eigen::VectorXd& vals);

		// Evaluate a real-valued expansion of the form (1) given the vector
		// (2) of its expansion coefficients.
		//
		// Input parameters:
		//   n - the number of terms in the expansion
		//   coefs - an array specifying the coefficients
		//   x - the point at which to evaluate (1)
		//
		// Output parameters:
		//   val - the value of the expansion at the point x
		double evaluate_expansion(
			const int n,
			const Eigen::VectorXd& coefs,
			const double x);

		// Evaluate one or more expansion of the form (1) and their
		// derivatives at a specified point given the vector (2) of coefficients
		//
		// Input parameters:
		//   n - the number of terms in the expansion
		//   coefs - an array specifying the coefficients
		//   x - the point at which to evaluate (1)
		//
		// Output parameters:
		//   pol - the value of the expansion at the point x
		//   der - the value of the derivative of the expansion at the point x
		void evaluate_expansion(
			const int n,
			const Eigen::VectorXd& coefs,
			const double x,
			double& pol,
			double& der);

		// Use barycentric Lagrange interpolation to evaluate a real-valued expansion of
		// the form (1) at a specified points.
		//
		// Input parameters:
		//   n - the number of terms in the expansion (1)
		//   (xslege,whtslege) - the nodes and weights of the n-point Guass-Legendre
		//     quadrature rule
		//   vals - the *scaled* values of (1) at the n Gauss-Legendre nodes
		//   x - the point at which to evaluate (1)
		//
		// Output parameters:
		//   valout - the value of (1) at the point x
		double interpolate(
			const int n,
			const Eigen::VectorXd& xslege,
			const Eigen::VectorXd& whtslege,
			const Eigen::VectorXd& vals,
			const double x);

		// Construct the matrix which takes the *scaled* values of (1) at the Gauss-Legendre
		// nodes to its scaled values at the nodes of a user-specified quadrature rule.
		//
		// Input parameters:
		//   n - the number of terms in the expansion (1)
		//   (xslege,whtslege) - the nodes and weights of the n-point Guass-Legendre
		//     quadrature rule
		//   m - the number of points in the output target quadrature rule
		//   (xsout,whtsout) - the nodes and weights of the output quadrature rule
		//
		// Output parameters:
		//   ainterp - the (m,n) interpolation matrix
		Eigen::MatrixXd interpolation_matrix(
			const int n,
			const Eigen::VectorXd& xslege,
			const Eigen::VectorXd& whtslege,
			const int m,
			const Eigen::VectorXd& xsout,
			const Eigen::VectorXd& whtsout);

		// Return the (n,n) matrix which takes the *scaled* values of an expansion of the
		// form (1) at the nodes of the n-point Gauss-Legendre quadrature to the
		// expansion coefficients.
		//
		// Input parameters:
		//   n - the number of terms in the expansion
		//   xslege - the nodes of the n-point Gauss-Legendre quadrature rule
		//   whtslege - the weights of the n-point Gauss-Legendre rule
		//
		// Output parameters:
		//   umatr - the (n,n) matrix which takes values to coefficients
		Eigen::MatrixXd coefficient_matrix(
			const int n,
			const Eigen::VectorXd& xslege,
			const Eigen::VectorXd& whtslege);

		// Compute the vector of coefficients of an expansion of the form (1)
		// in the monomial basis.
		//
		// Input parameters:
		//   n - the number of terms in the expansion
		//   coefs - an array specifying the coefficients (2)
		//
		// Output parameters:
		//   coefsout - an array specifying the coefficients in the monomial basis
		Eigen::VectorXd to_monomial(
			const int n,
			const Eigen::VectorXd& coefs);

		// Compute the vector (2) of coefficients of a polynomial from the vector
		// of coefficients in the monomial basis.
		//
		// Input parameters:
		//   n - the number of terms in the expansion
		//   coefs - an array specifying the coefficients in the monomial basis
		//
		// Output parameters:
		//   coefsout - an array specifying the coefficients (2)
		Eigen::VectorXd from_monomial(
			const int n,
			const Eigen::VectorXd& coefs);
	}
}

#endif // !NA_LEGENDRE_H
