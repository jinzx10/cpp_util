#ifndef __CXX_UTIL_MATH_HELPER_H__
#define __CXX_UTIL_MATH_HELPER_H__

#include <armadillo>
#include "arma_helper.h"
#include <functional>
#include <iostream>

namespace cxut {

    inline bool null_qr(arma::mat& ns, arma::mat const& A) {
        if (A.is_empty()) {
            ns.clear();
            return true;
        }
        arma::mat q, r;
        bool status = arma::qr(q, r, A.t());
        if (status) {
            arma::vec s = arma::sum(arma::abs(r), 1);
            double tol = arma::datum::eps * std::max(A.n_rows, A.n_cols);
            ns = q.cols(arma::find(s<tol));
        }
        return status;
    }
    
    inline arma::mat null_qr(arma::mat const& A) {
        arma::mat ns;
        bool status = null_qr(ns, A);
        if (!status) {
            std::cout << "null_qr(): qr decomposition failed." << std::endl;
            exit(EXIT_FAILURE);
        }
        return ns;
    }
    
    
    inline arma::mat orth_lowdin(arma::mat const& A) {
        arma::vec eigval;
        arma::mat eigvec;
        arma::eig_sym(eigval, eigvec, A*A.t());
        return arma::solve(eigvec*arma::diagmat(arma::sqrt(eigval))*eigvec.t(), A);
    }
    

    // find the smallest/largest number
    template <typename T1, typename T2>
    auto min(T1 const& i, T2 const& j) {
        return (i < j) ? i : j;
    }
    
    template <typename T1, typename T2, typename ...Ts>
    auto min(T1 const& i, T2 const& j, Ts const& ...args) {
        auto tmp = min(j, args...);
        return (i < tmp) ? i : tmp;
    }
    
    template <typename T1, typename T2>
    auto max(T1 const& i, T2 const& j) {
        return (i > j) ? i : j;
    }
    
    template <typename T1, typename T2, typename ...Ts>
    auto max(T1 const& i, T2 const& j, Ts const& ...args) {
        auto tmp = max(j, args...);
        return (i > tmp) ? i : tmp;
    }
    
    
    // first-order numerical differentiation by finite difference
    inline std::function<double(double)> grad(std::function<double(double)> const& f, double const& delta = 1e-3) {
        return [=] (double x) -> double {
            double m1 = f(x+delta) - f(x-delta);
            double m2 = f(x+2.0*delta) - f(x-2.0*delta);
            double m3 = f(x+3.0*delta) - f(x-3.0*delta);
            return (m3 - 9.0*m2 + 45.0*m1) / (60.0*delta);
        };
    }
    
    // (d/dxi) f(x0,x1,...)
    template <typename V>
    std::function<double(V)> gradi(std::function<double(V)> const& f, size_t const& i, double const& delta = 1e-3) {
        return [=] (V const& v) -> double {
            std::function<double(double)> g = [=, v=v] (double const& x) mutable {
                v[i] = x;
                return f(v);
            };
            return grad(g, delta)(v[i]); 
        };
    }
    
    // grad f(x0,x1,...)
    template <typename V>
    std::function<V(V)> grad(std::function<double(V)> const& f, double const& delta = 1e-3) {
        return [=] (V const& x) -> V {
            V df = x;
            for (size_t i = 0; i != x.size(); ++i) {
                df[i] = gradi(f, i, delta)(x);
            }
            return df;
        };
    }
    
    
    // generate 1D grid points according to some grid density
    inline arma::vec grid1d(double const& xmin, double const& xmax, std::function<double(double)> density, double const& x0) {
        if (x0 < xmin || x0 > xmax)
            return arma::vec{};
    
        arma::vec gr = {x0};
        double gr_last = gr.back();
        while ( gr_last < xmax ) {
            gr.resize(gr.n_elem+1);
            gr.back() = gr_last + 1.0 / density(gr_last);
            gr_last = gr.back();
        }
    
        arma::vec gl = {x0};
        double gl_last = gl.back();
        while ( gl_last > xmin ) {
            gl.resize(gl.n_elem+1);
            gl.back() = gl_last - 1.0 / density(gl_last);
            gl_last = gl.back();
        }
    
        return join_cols(flipud(gl.tail(gl.n_elem-1)), gr);
    }
    
    inline arma::vec grid1d(double const& xmin, double const& xmax, std::function<double(double)> density) {
        return grid1d(xmin, xmax, density, xmin);
    }

    
    // Broyden's quasi-Newton method
    // the good and bad Broyden's methods are identical in 1D
    inline int broydenroot(std::function<double(double)> f, double& x, double const& beta = 0.7, double const& tol = 1e-12, unsigned int const& max_iter = 50) {
        double fx = f(x);
        if (std::abs(fx) < tol)
            return 0;
    
        // compute the initial Jacobian by finite difference
        double delta = 1e-6 * std::max(1.0, std::sqrt(std::abs(x)));
        double J = (f(x+delta) - fx) / delta; 
    
        double dx = 0.0;
        double fx_new = 0.0;
        for (unsigned int counter = 1; counter != max_iter; ++counter) {
            if (std::abs(J) < 1e-14) {
                std::cout << "broydenroot: the Jacobian appears to be singular." << std::endl;
                return 2;
            }
    
            dx = -fx / J * beta;
            x += dx;
            fx_new = f(x);
            if (std::abs(fx_new) < tol)
                return 0;
    
            J = (fx_new - fx) / dx;
            fx = fx_new;
        }
    
        std::cout << "broydenroot: fails to find the root." << std::endl;
        return 1;
    }
    
    inline int broydenroot(std::function<arma::vec(arma::vec)> f, arma::vec& x, double const& beta = 0.7, double const& tol = 1e-12, unsigned int const& max_iter = 50, std::string const& method = "good") {
        arma::vec fx = f(x);
        if (arma::norm(fx) < tol)
            return 0;
    
        int md = -1;
        if (!method.compare("good"))
            md = 0;
        if (!method.compare("inv"))
            md = 1;
        if (!method.compare("bad"))
            md = 2;
    
        if ( md < 0 ) {
            std::cerr << "broydenroot: invalid method \"" << method << "\", use default instead. " << std::endl;
            md = 0;
        }
    
        arma::uword len_x = x.n_elem;
        arma::uword len_f = fx.n_elem;
    
        // compute the initial Jacobian by finite difference
        arma::mat J(len_f, len_x);
        double delta = 1e-6 * max(1.0, std::sqrt(arma::norm(x)));
        arma::vec dxi(len_x);
        arma::vec df(len_f);
        for (arma::uword i = 0; i != len_x; ++i) {
            dxi.zeros();
            dxi(i) = delta;
            df = f(x+dxi) - fx;
            J.col(i)  = df / delta;
        }
    
        arma::mat invJ;
        if ( md > 0 ) {
            if (len_x == len_f) {
                bool info = arma::inv(invJ, J);
                if (!info) {
                    std::cout << "broydenroot: the Jacobian appears to be singular." << std::endl;
                    return 2;
                }
            } else {
                std::cerr << "broydenroot: inverse update requires the number of equations equals to the number of variables." << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    
        arma::vec dx(len_x);
        arma::vec fx_new(len_f);
        for (unsigned int counter = 1; counter != max_iter; ++counter) {
            if (md > 0) {
                dx = -invJ * fx * beta;
            } else {
                dx = -arma::solve(J, fx) * beta;
            }
    
            x += dx;
            fx_new = f(x);
    
            if (arma::norm(fx_new) < tol)
                return 0;
    
            df = fx_new - fx;
            fx = fx_new;
    
            // Broyden's update
            switch (md) {
                case 2:
                    invJ += ( dx - invJ * df ) / arma::dot(df, df) * df.as_row();
                    break;
                case 1:
                    invJ += ( dx - invJ * df ) / arma::dot(dx, invJ * df) * dx.as_row() * invJ;
                    break;
                default:
                    J += ( df - J * dx ) / arma::dot(dx, dx) * dx.as_row();
            }
        }
    
        std::cout << "broydenroot: fails to find the root." << std::endl;
        return 1;
    }


    inline int diis(std::function<double(double)> iter, double& x, double tol = 1e-8, size_t const& max_iter = 50, size_t const& max_subspace = 20) {
        double xdiis = iter(x);
        double r = xdiis - x;
        x = xdiis;
        arma::rowvec xs = {x};
        arma::rowvec rs = {r};

        arma::mat B = rs.t() * rs;

        auto diismat = [&B] () -> arma::mat { 
            return arma::join_cols(
                    join_rows(B, arma::ones(B.n_rows, 1)),
                    join_rows(arma::ones(1, B.n_cols), arma::mat{0.0})
            ); 
        };

        auto diisvec = [&B] () -> arma::vec { 
            return join_cols(arma::zeros(B.n_cols), arma::vec{1.0}); 
        };

        size_t counter = 0;
        while (counter < max_iter) {
            if ( B.n_cols > max_subspace || arma::rcond(diismat()) < 1e-14 ) {
                rs.shed_col(0);
                xs.shed_col(0);
                B.shed_col(0);
                B.shed_row(0);
                continue;
            }
            ++counter;
            
            xdiis = arma::dot( xs, arma::solve(diismat(), diisvec()).eval().head_rows(B.n_cols) );

            x = iter(xdiis);
            r = x - xdiis;

            if (std::abs(r) < tol) {
                return 0;
            }

            xs.insert_cols(xs.n_cols, arma::vec{x});
            rs.insert_cols(rs.n_cols, arma::vec{r});

            B.resize(B.n_rows+1, B.n_cols+1);
            B.row(B.n_rows-1) = rs.col(rs.n_cols-1).t() * rs;
            B.col(B.n_cols-1) = B.row(B.n_rows-1).t();
        }

        std::cerr << "DIIS fails to converge." << std::endl;
        return 1;
    } 

    inline int diis(std::function< std::tuple<double, double>(double) > iter_err, double& x, double tol = 1e-8, size_t const& max_iter = 50, size_t const& max_subspace = 20) {
        double r;
        std::tie(x, r) = iter_err(x);
        arma::rowvec xs = {x};
        arma::rowvec rs = {r};

        arma::mat B = rs.t() * rs;

        auto diismat = [&B] () -> arma::mat { 
            return arma::join_cols(
                    join_rows(B, arma::ones(B.n_rows, 1)),
                    join_rows(arma::ones(1, B.n_cols), arma::mat{0.0})
            ); 
        };

        auto diisvec = [&B] () -> arma::vec { 
            return join_cols(arma::zeros(B.n_cols), arma::vec{1.0}); 
        };

        size_t counter = 0;
        while (counter < max_iter) {
            if ( B.n_cols > max_subspace || arma::rcond(diismat()) < 1e-14 ) {
                rs.shed_col(0);
                xs.shed_col(0);
                B.shed_col(0);
                B.shed_row(0);
                continue;
            }
            ++counter;
            
            std::tie(x, r) = iter_err( arma::dot( xs, arma::solve(diismat(), diisvec()).eval().head_rows(B.n_cols) ) );

            if (std::abs(r) < tol) {
                return 0;
            }

            xs.insert_cols(xs.n_cols, arma::vec{x});
            rs.insert_cols(rs.n_cols, arma::vec{r});

            B.resize(B.n_rows+1, B.n_cols+1);
            B.row(B.n_rows-1) = rs.col(rs.n_cols-1).t() * rs;
            B.col(B.n_cols-1) = B.row(B.n_rows-1).t();
        }

        std::cerr << "DIIS fails to converge." << std::endl;
        return 1;
    } 

    inline int diis(std::function<arma::vec(arma::vec)> iter, arma::vec& x, double tol = 1e-8, size_t const& max_iter = 50, size_t const& max_subspace = 20) {
        arma::vec xdiis = iter(x);
        arma::vec r = xdiis - x;
        x = xdiis;
        arma::mat xs = x;
        arma::mat rs = r;

        arma::mat B = rs.t() * rs;

        auto diismat = [&B] () -> arma::mat { 
            return arma::join_cols(
                    join_rows(B, arma::ones(B.n_rows, 1)),
                    join_rows(arma::ones(1, B.n_cols), arma::mat{0.0})
            ); 
        };

        auto diisvec = [&B] () -> arma::vec { 
            return join_cols(arma::zeros(B.n_cols), arma::vec{1.0}); 
        };

        size_t counter = 0;
        while (counter < max_iter) {
            if ( B.n_cols > max_subspace || arma::rcond(diismat()) < 1e-14 ) {
                rs.shed_col(0);
                xs.shed_col(0);
                B.shed_col(0);
                B.shed_row(0);
                continue;
            }
            ++counter;
            
            xdiis = xs * arma::solve(diismat(), diisvec()).eval().head_rows(B.n_cols);

            x = iter(xdiis);
            r = x - xdiis;

            if (arma::norm(r) < tol) {
                return 0;
            }

            xs.insert_cols(xs.n_cols, x);
            rs.insert_cols(rs.n_cols, r);

            B.resize(B.n_rows+1, B.n_cols+1);
            B.row(B.n_rows-1) = rs.col(rs.n_cols-1).t() * rs;
            B.col(B.n_cols-1) = B.row(B.n_rows-1).t();
        }

        std::cerr << "DIIS fails to converge." << std::endl;
        return 1;
    } 

    inline int diis(std::function< std::tuple<arma::vec, arma::vec>(arma::vec) > iter_err, arma::vec& x, double tol = 1e-8, size_t const& max_iter = 50, size_t const& max_subspace = 20) {
        arma::vec r;
        std::tie(x, r) = iter_err(x);
        arma::mat xs = x;
        arma::mat rs = r;

        arma::mat B = rs.t() * rs;

        auto diismat = [&B] () -> arma::mat { 
            return arma::join_cols(
                    join_rows(B, arma::ones(B.n_rows, 1)),
                    join_rows(arma::ones(1, B.n_cols), arma::mat{0.0})
            ); 
        };

        auto diisvec = [&B] () -> arma::vec { 
            return join_cols(arma::zeros(B.n_cols), arma::vec{1.0}); 
        };

        size_t counter = 0;
        while (counter < max_iter) {
            if ( B.n_cols > max_subspace || arma::rcond(diismat()) < 1e-14 ) {
                rs.shed_col(0);
                xs.shed_col(0);
                B.shed_col(0);
                B.shed_row(0);
                continue;
            }
            ++counter;
            
            std::tie(x, r) = iter_err(xs * arma::solve(diismat(), diisvec()).eval().head_rows(B.n_cols));

            if (arma::norm(r) < tol) {
                return 0;
            }

            xs.insert_cols(xs.n_cols, x);
            rs.insert_cols(rs.n_cols, r);

            B.resize(B.n_rows+1, B.n_cols+1);
            B.row(B.n_rows-1) = rs.col(rs.n_cols-1).t() * rs;
            B.col(B.n_cols-1) = B.row(B.n_rows-1).t();
        }

        std::cerr << "DIIS fails to converge." << std::endl;
        return 1;
    } 

}


#endif
