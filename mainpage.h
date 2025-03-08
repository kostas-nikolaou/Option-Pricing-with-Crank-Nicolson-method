/**
* 
* @mainpage Crank Nicolson scheme for pricing Vanilla Options
* 
* This class concerns the simple case of a vanilla option with no dividents. We seperate our options into European or American, Put or Call,
* and use the Crank Nicolson finite differnecing scheme to approximate the solution the Black Scholes equation. The class can also compute the five basic Greeks of the option
* (Delta, Gamma, Theta, Rho, and Vega), compare the difference in value between European and American option ceteris paribus, and validate the accuracy of the model using the explicit Black Scholes model.
* Additionally the class can handle time variant risk free interest rate. By adding at leat one value of the rate in any time point between computation date and maturity, the class will linearly interpolate the rest of the time points.
* 
* ## Getting Started 
* Simply include the class header, plus some common libraries, string.h, vector.h and utility.h (for using std::pairs). 
* 
* ```cpp
* #include<string>
* #include<utility>
* #include<vector>
* #include "CN.h"
* CN object(...);
* std::vector<std::vector<double>> optionPrices = object.solve();
* ```
* 
* ## Crank Nicolson finite differencing method
* @image html grid.png
*
* The Crank-Nicholson method is widely used for solving time-dependent partial differential equations, particularly the heat equation and
* diffusion problems. It is an implicit, finite difference method that combines the features of both the forward Euler and backward Euler methods, balancing
* between stability and accuracy. Its key features are:
*
* - **Unconditional Stability:** The method is stable regardless of time step size.
* - **Second-Order Accuracy:** The scheme achieves second-order accuracy in both space and time.
* - **Implicit Time-Stepping:** Unlike explicit methods, it allows for larger time steps without instability.
* 
* It transforms the partial differential equation into a system of linear equations that can be
* efficiently solved using a tridiagonal matrix algorithm.
*
* The Black-Scholes equation is given by:
* \f[
* \frac{\partial V}{\partial t} + \frac{1}{2} \sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + rS \frac{\partial V}{\partial S} - rV = 0,
* \f]
* where:
* - \f$ V(S, t) \f$ is the option price as a function of stock price \f$ S \f$ and time \f$ t \f$.
* - \f$ \sigma \f$ is the volatility.
* - \f$ r \f$ is the risk-free interest rate.
*
* Using the Crank-Nicolson method, the discretized equation at grid point \f$ i \f$ becomes:
* \f[
* a_i V_{i-1}^{n+1} + b_i V_i^{n+1} + c_i V_{i+1}^{n+1} = d_i,
* \f]
* where:
* \f{align*}{
* a_i &= -\frac{\Delta t}{4} \left( \sigma^2 i^2 - ri \right), \\
* b_i &= 1 + \frac{\Delta t}{2} \left( \sigma^2 i^2 + r \right), \\
* c_i &= -\frac{\Delta t}{4} \left( \sigma^2 i^2 + ri \right), \\
* d_i &= \frac{\Delta t}{4} \left( \sigma^2 i^2 - ri \right) V_{i-1}^n \\
*      &\quad + \left( 1 - \frac{\Delta t}{2} \left( \sigma^2 i^2 + r \right) \right) V_i^n \\
*      &\quad + \frac{\Delta t}{4} \left( \sigma^2 i^2 + ri \right) V_{i+1}^n.
* \f}
*
* This leads to a tridiagonal matrix system that can be solved efficiently using the Thomas algorithm. To ensure stability, the following constraint should be met:
* \f[
* dt < \frac{dS^{2}}{2* \sigma^{2} * S_{max} ^{2}}
* \f]
* 
* The Option Greeks are calculated using simple forward or central differencing schemes of the options's value over the respective difference or pertrubation.
*
*
* Advantages:
* - Coverges in Error faster than Implicit/Explicit method as number of timesteps increases
* - Suitable for large time step simulations without instability.
* - Efficient when solved using methods like the Thomas algorithm for tridiagonal systems.
* 
* @image html convergence.png
*
* Reference for both images [here](https://quintus-zhang.github.io/post/on_pricing_options_with_finite_difference_methods/).
* 
* See Also:
* - [Black-Scholes Equation](https://en.wikipedia.org/wiki/Black%E2%80%93Scholes_equation)
* - [Numerical Methods Overview](https://en.wikipedia.org/wiki/Numerical_methods)
* - [Crank-Nicolson Method](https://en.wikipedia.org/wiki/Crank%E2%80%93Nicolson_method)
* - [Option Greeks](https://shorturl.at/x9uWg)
* - [On pricing options with finite difference methods](https://quintus-zhang.github.io/post/on_pricing_options_with_finite_difference_methods/)
*/

