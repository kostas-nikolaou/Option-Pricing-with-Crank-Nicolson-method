#ifndef CN_H
#define CN_H
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <utility>
#include <numeric>
#include <sstream>
#include <iomanip>

/**
* @class CN
* @brief This class performs the Crank Nicolson finite differencing scheme for pricing European and American Option, for non-constant risk free interest rate.
*
* This class is able to price Vanilla Options with the Crank Nicolson finite differencing numerical method, for a non constant interest rate,
* using linear interpolation over the points in time given. It does for no dividents, meaning the American and european values shoudl be the same and can be compared
* to the Black Scholes explicit model. The class provides the ability to compare the two types of options, and validate performs using the explicit method.
*/
class CN {
    std::string option; /// Option contract, European or American
    std::string type; ///Type of exercise, Call or Put
    double S_0;   /// Initial asset price
    double K;        /// Strike price
    double Tf;       /// Time to maturity in years
    double T0;       /// Initial time
    std::vector<std::pair<double, double>> r; /// Risk-free interest rate
    double sigma;   /// Volatility
    int Nx;         /// Number of spatial points (asset prices)
    int Nt;        /// Number of time steps


public:
    /**
    * @brief Constructor for the class
    *  @param op String, The type of option, American or European
    *  @param type_ String, The exercise type Call or Put
    * @param s0 Double, The underlying asset (Spot) price today
    * @param k Double, The strike price
    * @param r_ A vector of pairs containing a double representing a date between T0 and Maturity, and a double representing the risk-free interest rate at that time.
    * @param s_ Double, The volatility of the underlying asset
    * @param nx Integer for the spot grid points. Default price is 100.
    * @param nt Integer for the time grid points. Default price is 100.
    * @param tf Double, Maturity time. Default time is 1.0.
    * @param t0 Double, Computation time. Default time is 0.0 (today).
    * */
    CN(std::string op, std::string type_, double s0, double k, std::vector<std::pair<double, double>> r_, double s_, int nx = 100, int nt = 100, double tf = 1.0, double t0 = 0.0) : option(op),
        type(type_), S_0(s0), K(k), Tf(tf), T0(t0), r(r_), sigma(s_), Nx(nx), Nt(nt) {
        if (type_ != "Call" && type_ != "Put") {
            throw std::invalid_argument("Option type must be 'Call' or 'Put'.");
        }
        if (op != "European" && op != "American") {
            throw std::invalid_argument("Option must be 'European' or 'American'.");
        }
        if (s0 < 0 || k < 0 || tf < 0 || t0 < 0) {
            throw std::invalid_argument("Variables S0, K, T0 and Tf must non negative.");
        }
        if (t0 >= tf) {
            throw std::invalid_argument("Final time must be higher than initial time.");
        }

        if (nx == 0 || nt == 0) {
            throw std::invalid_argument("Number of steps of both price and time must be positive.");
        }
        if (s_ <= 0) {
            throw std::invalid_argument("Volatility of underlying asset must be positive.");

        }
        if (std::max_element(r.begin(), r.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; })->first > Tf) {
            throw std::invalid_argument("Risk Free interest rate must correspond to maturity.");
        }
        if (std::min_element(r.begin(), r.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; })->first < T0) {
            throw std::invalid_argument("Risk Free interest rate must correspond to computation date.");
        }

    }
    /** 
    * @brief Destructor of the class
    *
    * */
    ~CN() { std::cout << "deleted data" << std::endl; }
    /**
    * @brief Copy Constructor of the class
    */
    CN(CN& b) {
        option = b.option;
        type = b.type;
        S_0 = b.S_0;
        K = b.K;
        Tf = b.Tf;
        T0 = b.T0;
        r = b.r;
        sigma = b.sigma;
        Nx = b.Nx;
        Nt = b.Nt;
    }
    /**
    * @brief Takes a vector of pairs from the constructor and linearly interpolates for the entire time span, creating a varying risk-free interest rate.
    * @param dr_ Double, represents a pertrubation in the interest rate. Used for deriving the Rho of the option, default value is zero.
    * @return The interpolated vector with an interest rate for every time step
    */
    std::vector<double> defineRate(double dr_ = 0) {
        std::vector<double> result(Nt, 0.0);
        double dt = (Tf - T0) / Nt;

        for (int i = 0; i < Nt; i++) {
            double t = T0 + i * dt;
            auto it = std::find_if(r.begin(), r.end(),
                [t](const std::pair<double, double>& pair) { return std::abs(pair.first - t) < 1e-9; });

            if (it != r.end()) {
                result[i] = it->second + dr_;
                continue;
            }
            auto upper = std::upper_bound(r.begin(), r.end(), t,
                [](double val, const auto& pair) { return val < pair.first; });

            if (upper == r.end()) {
                result[i] = std::prev(r.end())->second + dr_;;
                continue;
            }
            if (upper == r.begin()) {
                result[i] = r.begin()->second + dr_;
                continue;
            }

            auto lower = std::prev(upper);
            double x1 = lower->first, v1 = lower->second;
            double x2 = upper->first, v2 = upper->second;
            result[i] = dr_ + v1 + (v2 - v1) * (t - x1) / (x2 - x1);

        }
        return result;
    }
    /** @brief Solves a tridiagonal matrix using the Thomas Algorithm
    * @param lower The lower diagonal
    * @param diag The middle diagonal
    * @param upper The upper diagonal
    * @param b The right hand part vector
    * @param result An empty vector which is populated with the result
    */
    void solveTridiagonal(const std::vector<double>& lower, const std::vector<double>& diag,
        const std::vector<double>& upper, const std::vector<double>& b,
        std::vector<double>& result) {
        int n = diag.size();
        std::vector<double> c_prime(n, 0.0);
        std::vector<double> d_prime(n, 0.0);
        result.resize(n);

        //forward step
        c_prime[0] = upper[0] / diag[0];
        d_prime[0] = b[0] / diag[0];
        for (int i = 1; i < n; ++i) {
            double m = diag[i] - lower[i] * c_prime[i - 1];
            c_prime[i] = upper[i] / m;
            d_prime[i] = (b[i] - lower[i] * d_prime[i - 1]) / m;
        }

        //backsubstitution
        result[n - 1] = d_prime[n - 1];
        for (int i = n - 2; i >= 0; --i) {
            result[i] = d_prime[i] - c_prime[i] * result[i + 1];
        }
    }
    /**
    * @brief Computes the mean of the second part of a vector of pairs
    * @param vec The vector of pairs containing doubles
    * @return Double, the arithmetic mean
    */
    double meanSecondElement(const std::vector<std::pair<double, double>>& vec) { //calculate the mean of the input of interest rates
        if (vec.empty()) {
            throw std::invalid_argument("Vector is empty. Cannot compute mean.");
        }
        double sum = 0.0;
        for (int i = 0; i < vec.size(); ++i) {
            sum += vec[i].second;
        }
        
        return sum / vec.size();
    }

    /** @brief Solves the Black Schole PDE using the Crank Nicholson finite differencing method
    * @param dP Pertrubation in price, used for Delta computation
    * @param ds_ Pertrubation in volatility, used for Vega computation
    * @param dr_ Pertrubation in interest rate, passed to DefineRate, used for Rho computation
    * @return Vector of vectors with the solution of the option pricing
    */
    std::vector<std::vector<double>> solve(double dP = 0, double ds_ = 0, double dr_ = 0) {
        double T = Tf - T0; //Time Frame
        double dS = 2.0*S_0 / Nx; //Price step
        double dt = T / Nt;  //Time step
        bool isCall = (type == "Call");
        bool isAmerican = (option == "American");
        std::vector<double> R = defineRate(dr_);
        double Sigma = sigma + ds_;
        
        //price grid
        std::vector<double> S(Nx + 1, 0.0);
        for (int i = 0; i <= Nx; i++) {
            S[i] =  i * dS + dP;
            S[i] *= std::exp(R[0] * T0); //Assuming rate is constant from today to computation date 
        }
        if (dt > (dS * dS) / (2 * Sigma * Sigma * S[Nx] * S[Nx])) {
            std::cerr << "Warning: Time step too large, stability may be comprimised." << std::endl;
        }
        
        std::vector<std::vector<double>> V(Nx + 1, std::vector<double>(Nt + 1, 0.0));

        //maturity payoff
        for (int i = 0; i <= Nx; i++) {
            V[i][Nt] = isCall ? std::max(S[i] - K, 0.0) : std::max(K - S[i], 0.0);
            V[i][Nt] *= std::exp(-trapezoidal(R, T0, Tf));  
        };
        //Boundary conditions- Dirichlet
        std::vector<double> subvec; // je comprends pas pourquoi vous aimez autant cette langue 
        for (int t = 0; t < Nt; t++) {
            subvec.clear();
            std::copy(R.begin(), R.begin() + t, std::back_inserter(subvec));
            V[0][t] = isCall ? 0.0 : std::max(K - S[0], 0.0) * std::exp(-trapezoidal(subvec, 0, t * dt));
            V[Nx][t] = isCall ? std::max(S[Nx] - K, 0.0) * std::exp(-trapezoidal(subvec, 0, t * dt)) : 0.0;    
        }
        
        //Coefficients
        std::vector<double> alpha(Nx + 1, 0.0), beta(Nx + 1, 0.0), gamma(Nx + 1, 0.0);

        //Tridiagonal matrix components
        std::vector<double> lower(Nx - 1, 0.0), main(Nx - 1, 0.0), upper(Nx - 1, 0.0);

        //Time step
        for (int t = Nt - 1; t >= 0; --t) {
            for (int i = 0; i <= Nx; ++i) {
                alpha[i] = 0.25 * dt * ((Sigma * Sigma * S[i] * S[i]) / (dS * dS) - (R[t] * S[i]) / dS);
                beta[i] = -0.5 * dt * ((Sigma * Sigma * S[i] * S[i]) / (dS * dS) + R[t]);
                gamma[i] = 0.25 * dt * ((Sigma * Sigma * S[i] * S[i]) / (dS * dS) + (R[t] * S[i]) / dS);
            }

            //Boundary adjustments
            std::vector<double> b(Nx - 1, 0.0);

            for (int i = 1; i < Nx; ++i) {
                lower[i - 1] = -alpha[i];
                main[i - 1] = 1.0 - beta[i];
                upper[i - 1] = -gamma[i];
            }
            for (int i = 1; i < Nx - 2; ++i) {
                b[i] = alpha[i + 1] * V[i][t + 1] + (1.0 + beta[i + 1]) * V[i + 1][t + 1] + gamma[i + 1] * V[i + 2][t + 1];
            }

            // Apply Dirichlet boundary conditions
            b[0] = (1.0 + beta[1]) * V[1][t + 1] + gamma[1] * V[2][t + 1] + alpha[1] * (V[0][t + 1] + V[0][t]) - alpha[0] * V[0][t + 1];  // Left boundary
            b[Nx - 2] = alpha[Nx - 1] * V[Nx - 2][t + 1] + (1 + beta[Nx - 1]) * V[Nx - 1][t + 1] + gamma[Nx - 1] * (V[Nx][t + 1]+V[Nx][t]);  // Right boundary

            //Solution
            std::vector<double> solution;
            solveTridiagonal(lower, main, upper, b, solution);

            //Update
            for (int i = 1; i < Nx; i++) {
                V[i][t] = solution[i - 1];
                
            }
        }
        
         if (isAmerican) {
              for (int t = Nt - 1; t >= 0; --t) {
                subvec.clear();
                std::copy(R.begin(), R.begin() + t, std::back_inserter(subvec));
                for (int i = 1; i < Nx; i++) {
                    V[i][t] = std::max(V[i][t], isCall ? std::max(S[i] - K, 0.0) * std::exp(-trapezoidal(subvec, 0, t * dt)) : std::max(K - S[i], 0.0) * std::exp(-trapezoidal(subvec, 0, t * dt)));

                  }
              }
         }
         
        return V;
    }
    
    /** @brief Calls the solve function and saves the result on a .csv file
    */
    void solve_and_save() {
        double dS = 2.0*S_0 / Nx; //price step
        std::vector<double> R = defineRate();
        std::vector<std::vector<double>> V = solve();
        std::vector<double> S(Nx + 1, 0.0);
        for (int i = 0; i <= Nx; ++i) {
            S[i] =  i * dS;
            S[i] *= std::exp(R[0] * T0); //assume rate is constant from today to computation date
        }
        std::ofstream output("black_scholes_solution.csv");
        if (!output.is_open()) {
            std::cerr << "Error: Could not open file for writing." << std::endl;
            return;
        }
        output << "Stock Price";
        for (size_t t = 0; t < V[0].size(); ++t) {
            output << ",Timestep_" << t;
        }
        output << "\n";

        //option price for each stock price
        for (size_t i = 0; i < S.size(); ++i) {
            output << S[Nx - i];
            for (size_t t = 0; t < V[i].size(); ++t) {
                output << "," << V[i][t];
            }
            output << "\n";
        }

        output.close();
        std::cout << "Results saved." << std::endl;

    }
    /** @brief Extract the computation's date option price at time t.
    * @param t Double, time to extract the target. Must be between T0 and Tf. It is rounded to the index corresponding to the nearest time step 
    */
    std::vector<double> OptionPrice_At(double t = 0) { 
        if (t < 0 || t>Tf) {
            throw "Date cannot be before computation date or after maturity.";
        }
        if (t == 0) {}
        double dt = (Tf - T0) / Nt;
        int index = std::floor((t - T0) / dt);
        std::vector<double> result(Nx + 1, 0.0);
        std::vector<std::vector<double>> target = solve();
        if (t == Tf) {
            for (int i = 0; i <= Nx; i++) {
                result[i] = target[i][index];
            }
            return result;
        }
        for (int i = 0; i <= Nx; i++) {
            result[i] = target[i][index] + ((t - index * dt) / ((index + 1) * dt - index * dt)) * (target[i][index + 1] - target[i][index]); //linear interpolation
        }
        return result;
    }
    /** @brief Extract the computation's date option Delta at time t. Uses central differencing.
    * @param t Double, time to extract the target. Must be between T0 and Tf. It is rounded to the index corresponding to the nearest time step. Default Value is 0.0.
    * @param dP Double, the pertrubation for the differencing scheme. Default Value is 0.5.
    */
    std::vector<double> Delta_At(double t = 0, double dP = 0.5) { //use this for delta at t
        if (t < 0 || t>Tf) {
            throw "Date cannot be before computation date or after maturity.";
        }
        if (t == 0) {}
        double dt = (Tf - T0) / Nt;
        int index = std::floor((t - T0) / dt);
        std::vector<double> result(Nx + 1, 0.0);
        std::vector<std::vector<double>> target = delta(dP);
        if (t == Tf) {
            for (int i = 0; i <= Nx; i++) {
                result[i] = target[i][index];
            }
            return result;
        }
        for (int i = 0; i <= Nx; i++) {
            result[i] = target[i][index] + ((t - index * dt) / ((index + 1) * dt - index * dt)) * (target[i][index + 1] - target[i][index]); //linear interpolation
        }
    }
    /** @brief Extract the computation's date option Gamma at time t. Uses central differencing.
    * @param t Double, time to extract the target. Must be between T0 and Tf. It is rounded to the index corresponding to the nearest time step
    * @param dP Double, the pertrubation for the differencing scheme. Default Value is 0.5.
    */
    std::vector<double> Gamma_At(double t = 0, double dP = 0.5) { //use this for gamma at t
        if (t < 0 || t>Tf) {
            throw "Date cannot be before computation date or after maturity.";
        }
        if (t == 0) {}
        double dt = (Tf - T0) / Nt;
        int index = std::floor((t - T0) / dt);
        std::vector<double> result(Nx + 1, 0.0);
        std::vector<std::vector<double>> target = gamma(dP);
        if (t == Tf) {
            for (int i = 0; i <= Nx; i++) {
                result[i] = target[i][index];
            }
            return result;
        }
        for (int i = 0; i <= Nx; i++) {
            result[i] = target[i][index] + ((t - index * dt) / ((index + 1) * dt - index * dt)) * (target[i][index + 1] - target[i][index]); //linear interpolation
        }
        return result;
    }
    /** @brief Extract the computation's date option Rho at time t. Uses central differencing.
    * @param t Double, time to extract the target. Must be between T0 and Tf. It is rounded to the index corresponding to the nearest time step
    * @param dr_ Double, the pertrubation for the differencing scheme. Default Value is 0.01 (or 1%).
    */
    std::vector<double> Rho_At(double t = 0, double dr_ = 0.01) { //use this for rho at t
        if (t < 0 || t>Tf) {
            throw "Date cannot be before computation date or after maturity.";
        }
        if (t == 0) {}
        double dt = (Tf - T0) / Nt;
        int index = std::floor((t - T0) / dt);
        std::vector<double> result(Nx + 1, 0.0);
        std::vector<std::vector<double>> target = rho(dr_);
        if (t == Tf) {
            for (int i = 0; i <= Nx; i++) {
                result[i] = target[i][index];
            }
            return result;
        }
        for (int i = 0; i <= Nx; i++) {
            result[i] = target[i][index] + ((t - index * dt) / ((index + 1) * dt - index * dt)) * (target[i][index + 1] - target[i][index]); //linear interpolation
        }
        return result;
    }
    /** @brief Extract the computation's date option Vega at time t. Uses central differencing.
    * @param t Double, time to extract the target. Must be between T0 and Tf. It is rounded to the index corresponding to the nearest time step
    * @param ds_ Double, the pertrubation for the differencing scheme. Default Value is 0.1.
    */
    std::vector<double> Vega_At(double t = 0, double ds_ = 0.1) { //use this for vega at t
        if (t < 0 || t>Tf) {
            throw "Date cannot be before computation date or after maturity.";
        }
        if (t == 0) {}
        double dt = (Tf - T0) / Nt;
        int index = std::floor((t - T0) / dt);
        std::vector<double> result(Nx + 1, 0.0);
        std::vector<std::vector<double>> target = vega(ds_);
        if (t == Tf) {
            for (int i = 0; i <= Nx; i++) {
                result[i] = target[i][index];
            }
            return result;
        }
        for (int i = 0; i <= Nx; i++) {
            result[i] = target[i][index] + ((t-index*dt)/((index+1)*dt -index*dt))*(target[i][index+1]-target[i][index]); //linear interpolation
        }
        return result;
    }
    /** @brief Extract the computation's date option Theta at time t. Uses central differencing. It uses the Time Steps (Tf-T0)/Nt for differencing.
    * @param t Double, time to extract the target. Must be between T0 and Tf. It is rounded to the index corresponding to the nearest time step. Default Value is 0.0.
    * 
    */
    std::vector<double> Theta_At(double t = 0) { //use this for theta at t
        if (t < 0 || t>Tf) {
            throw "Date cannot be before computation date or after maturity.";
        }
        if (t == 0) {}
        double dt = (Tf - T0) / Nt;
        int index = std::floor((t - T0) / dt);
        std::vector<double> result(Nx + 1, 0.0);
        std::vector<std::vector<double>> target = theta();
        if (t >= Tf - dt && t <= Tf) {
            for (int i = 0; i <= Nx; i++) {
                result[i] = target[i][index];
            }
            return result;
        }
        for (int i = 0; i <= Nx; i++) {
            result[i] = target[i][index] + ((t - index * dt) / ((index + 1) * dt - index * dt)) * (target[i][index + 1] - target[i][index]); //linear interpolation
        }
        return result;
    }
    /** @brief Return the Cumulative Distribution Function(CDF) for a Standard Normal Distribution of double x.
    * @param x Double
    */
    double normCdf(double x) {
        return 0.5 * (1.0 + erf(x / sqrt(2.0)));
    }
    /** @brief Uses the trapezoidal approach for approximating an integral.
    * @param f Vector of arbitary length, the function over which to integrate
    * @param a Double, lower limit of integration
    * @param b Double, upper limit if integration
    */
    double trapezoidal(const std::vector<double>& f, double a, double b) {
        int n = f.size() - 1;
        if (n <= 0) {
            return 0.0;

        }
        double h = (b - a) / n; //size
        double integral = 0.0;
        for (int i = 0; i < n; i++) {
            integral += (f[i] + f[i + 1]) * h / 2.0;
        }

        return integral;
    }
    /** @brief Uses the Black Scholes explicit model to price the option.
    * 
    * Due to the lemma that states that how an european and an american option call/put with no dividents are equal when r is positive/negative, we just compute the european option using the black-scholes model(you can confirm the lemma by suing method compare us_and_eu()).
    */
    std::pair<double,int> validate() { 
        std::vector<double> d1(Nx + 1, 0.0);
        std::vector<double> d2(Nx + 1, 0.0);
        double T = (Tf - T0); //in years
        double dS = 2.0 * S_0 / Nx; //price step
        double dt = T / Nt;  //time step
        bool isCall = (type == "Call");
        bool isAmerican = (option == "American");
        std::vector<double> R = defineRate();
        std::vector<std::vector<double>> target = solve();
        double mean = 0; //For MAE
        int flags = 0; //For finding how mnay points have high error
        std::vector<double> S(Nx + 1, 0.0); //spot price
        for (int i = 0; i <= Nx; ++i) {
            S[i] = i * dS;
            S[i] *= std::exp(R[0] * T0); //compound for computation date
        }
        std::vector<std::vector<double>> V(Nx + 1, std::vector<double>(Nt + 1, 0.0));
        //maturity payoff
        for (int i = 0; i <= Nx; ++i) {
            V[i][Nt] = isCall ? std::max(S[i] - K, 0.0) : std::max(K - S[i], 0.0);
            V[i][Nt] *= std::exp(-trapezoidal(R, T0, Tf));
            mean += std::fabs(V[i][Nt] - target[i][Nt]);
            flags = (std::fabs(V[i][Nt] - target[i][Nt]) > std::max(0.05*V[i][Nt], 1e-6)) ? flags + 1 : flags;
        }
        
        
        std::vector<double> subvec_1;
        std::vector<double> subvec;
        for (int t = Nt - 1; t >= 0; --t) {
            double time = t * dt;
            subvec.clear();
            std::copy(R.begin() + t, R.end(), std::back_inserter(subvec));
            subvec_1.clear();
            std::copy(R.begin(), R.begin() + t, std::back_inserter(subvec_1));
            for (int i = 0; i <= Nx; i++) {
                d1[i] = (1 / (sigma * sqrt(T-time))) * (std::log(S[i] / K) + (sigma*sigma / 2) * (T-time));
                d2[i] = d1[i] - sigma * sqrt(T-time);
                V[i][t] = isCall ? normCdf(d1[i]) * S[i] - normCdf(d2[i]) * K  : normCdf(-d2[i]) * K  - normCdf(-d1[i]) * S[i];
                V[i][t] *= std::exp(-trapezoidal(subvec_1, 0, time));
                
                mean += std::fabs(V[i][t] - target[i][t]);
                flags = (std::fabs(V[i][t] - target[i][t]) > std::max(0.05*V[i][t], 1e-6)) ? flags + 1 : flags;
                
            }
        }
        subvec_1.clear();
        mean /= (Nx + 1) * (Nt + 1);
        std::cout << "Validation finished" << std::endl;
        return { mean,flags };
    }
    /** @brief Perform time stepping for european call/put and american call/put to prove they are the same for r>0/r<0. 
    *  We assume that the risk free interest rate remains strictly positive/negative for a call/put.
    */
    double compare_eu_us_nc() { 
        if (std::max_element(r.begin(), r.end(),
            [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
                return a.second < b.second; //Compare based on the second value
            })->second * std::min_element(r.begin(), r.end(),
                [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
                    return a.second < b.second; //Compare based on the second value
                })->second < 0) {
            std::cerr << "Warning: Risk free interest rate must have the same sign in order for the lemma to hold,and for European and American option of the same type to have identical prices." << std::endl;

        }
        if (meanSecondElement(r) >= 0.0) {
            double mean = 0.0;
            //int flags = 0;
            CN c1("European", "Call", S_0, K, r, sigma, Nx, Nt, Tf, T0);
            CN c2("American", "Call", S_0, K, r, sigma, Nx, Nt, Tf, T0);
            std::vector<std::vector<double>> v1 = c1.solve();
            std::vector<std::vector<double>> v2 = c2.solve();
            for (int t = Nt; t >= 0; --t) {
                for (int i = 0; i <= Nx; i++) {
                    mean += std::fabs(v1[i][t] - v2[i][t]);
                }

            }
            mean /= (Nx + 1) * (Nt + 1);
            std::cout << "Comparison finished." << std::endl;
            return mean;
        }
        else {
            double mean = 0;
            CN c1("European", "Put", S_0, K, r, sigma, Nx, Nt, Tf, T0);
            CN c2("American", "Put", S_0, K, r, sigma, Nx, Nt, Tf, T0);
            std::vector<std::vector<double>> v1 = c1.solve();
            std::vector<std::vector<double>> v2 = c2.solve();
            for (int t = Nt; t >= 0; --t) {
                for (int i = 0; i <= Nx; i++) {
                    mean += std::fabs(v1[i][t] - v2[i][t]);
                }

            }
            mean /= (Nx + 1) * (Nt + 1);
            std::cout << "Comparison finished." << std::endl;
            return mean;
        }
        
    }
    /** @brief Extract the computation's date option Delta. Uses central differencing.
    * @param dP Double, the pertrubation for the differencing scheme. Default Value is 0.5.
    */
    std::vector<std::vector<double>> delta(double dP=0.5) {
        double dS = dP;
        std::vector<std::vector<double>> optionPriceUp = solve(dP, 0.0, 0.0);
        std::vector<std::vector<double>> optionPriceDown = solve(-dP, 0, 0);
        size_t n = optionPriceUp.size();
        size_t m = optionPriceUp[0].size();
        std::vector<std::vector<double>> result(n, std::vector<double>(m, 0.0));

        for (size_t i = 0; i < (n); i++) {
            for (size_t j = 0; j < (m); j++) {
                result[i][j] = (optionPriceUp[i][j] - optionPriceDown[i][j]) / 2 * dS;
            }
        }

        return result;
    }
    /** @brief Extract the computation's date option Gamma. Uses central differencing.
    * @param dP Double, the pertrubation for the differencing scheme. Default Value is 0.5.
    */
    std::vector<std::vector<double>>  gamma(double dP=0.5) {
        double dS = dP;
        std::vector<std::vector<double>> optionPrice = solve(0.0, 0.0, 0.0);
        std::vector<std::vector<double>> optionPriceUp = solve(dP, 0.0, 0.0);
        std::vector<std::vector<double>> optionPriceDown = solve(-dP, 0.0, 0.0);
        size_t n = optionPrice.size();
        size_t m = optionPrice[0].size();
        std::vector<std::vector<double>> result(n, std::vector<double>(m, 0.0));
        for (size_t i = 1; i < (n); i++) {
            for (size_t j = 0; j < (m); j++) {
                result[i][j] = (optionPriceUp[i][j] - 2 * optionPrice[i][j] + optionPriceDown[i][j]) / (dS * dS);
            }
        }
        return result;
    }
    /** @brief Extract the computation's date option Theta. Uses forward differencing.
    *
    */
    std::vector<std::vector<double>> theta() {
        double dt = (Tf - T0) / Nt;
        std::vector<std::vector<double>> optionPrice = solve();
        size_t n = optionPrice.size();
        size_t m = optionPrice[0].size();
        std::vector<std::vector<double>> result(n, std::vector<double>(m - 1, 0.0));
        for (size_t i = 0; i < (n); i++) {
            for (size_t j = 0; j < (m - 1); j++) {
                result[i][j] = (optionPrice[i][j + 1] - optionPrice[i][j]) / (dt);
            }
        }
        return result;
    }
    /** @brief Extract the computation's date option Rho. Uses central differencing.
    * @param dr Double, the pertrubation for the differencing scheme. Default Value is 0.01 (or 1%).
    */
    std::vector<std::vector<double>> rho(double dr = 0.01) {
        if (dr <= 0.0) {
            throw std::invalid_argument("differencing must be positive to define rho, eg. 1 basis point (0.0001).");
        }
        std::vector<std::vector<double>> optionPriceUp = solve(0.0, 0.0, dr);
        std::vector<std::vector<double>> optionPriceDown = solve(0.0, 0.0, -dr);
        size_t n = optionPriceUp.size();
        size_t m = optionPriceUp[0].size();
        std::vector<std::vector<double>> result(n, std::vector<double>(m, 0.0));
        for (size_t i = 0; i < (n); i++) {
            for (size_t j = 0; j < m; j++) {
                result[i][j] = (optionPriceUp[i][j] - optionPriceDown[i][j]) / (2 * dr); //central differencing
            }
        }
        return result;
    }
    /** @brief Extract the computation's date option Vega. Uses central differencing.
    * @param dSigma Double, the pertrubation for the differencing scheme. Default Value is 0.1.
    */
    std::vector<std::vector<double>> vega(double dSigma = 0.1) {
        if (dSigma <= 0.0) {
            throw std::invalid_argument("differencing must be positive to define vega, eg. 1 basis point (0.0001).");
        }
        std::vector<std::vector<double>> optionPriceUp = solve(0.0, dSigma, 0.0);
        std::vector<std::vector<double>> optionPriceDown = solve(0.0, -dSigma, 0.0);
        size_t n = optionPriceUp.size();
        size_t m = optionPriceUp[0].size();
        std::vector<std::vector<double>> result(n, std::vector<double>(m, 0.0));
        for (size_t i = 0; i < (n); i++) {
            for (size_t j = 0; j < m; j++) {
                result[i][j] = (optionPriceUp[i][j] - optionPriceDown[i][j]) / (2 * dSigma); //central differencing
            }
        }
        return result;
    }
   
};

#endif
