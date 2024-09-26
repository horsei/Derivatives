import math
import numpy as np

# Function to generate normally distributed random numbers
def generate_gaussian_noise(mean, stddev):
    return np.random.normal(mean, stddev)

# Function to calculate the payoff of a European call option
def call_option_payoff(S, K):
    return max(S - K, 0.0)

# Function to calculate the payoff of a European put option
def put_option_payoff(S, K):
    return max(K - S, 0.0)

# Monte Carlo Simulation for European option pricing
def monte_carlo_option_pricing(S0, K, r, sigma, T, num_simulations, is_call_option):
    payoff_sum = 0.0

    for i in range(num_simulations):
        # Generate a random price path
        ST = S0 * math.exp((r - 0.5 * sigma ** 2) * T + sigma * math.sqrt(T) * generate_gaussian_noise(0.0, 1.0))

        if is_call_option:
            payoff = call_option_payoff(ST, K)
        else:
            payoff = put_option_payoff(ST, K)

        payoff_sum += payoff

    # Calculate the average payoff and discount it to present value
    average_payoff = payoff_sum / num_simulations
    current_value = math.exp(-r * T) * average_payoff
    estimated_future_value = average_payoff
    
    return current_value, estimated_future_value

if __name__ == "__main__":
    # Option parameters
    S0 = 100.0   # Initial stock price
    K = 110.0    # Strike price
    r = 0.05     # Risk-free rate
    sigma = 0.2  # Volatility
    T = 1        # Time to maturity (1 year)
    num_simulations = 1000000  # Number of simulations
    
    # Calculate option prices
    call_current_value, call_estimated_future_value = monte_carlo_option_pricing(S0, K, r, sigma, T, num_simulations, True)
    put_current_value, put_estimated_future_value = monte_carlo_option_pricing(S0, K, r, sigma, T, num_simulations, False)
    
    # Calculate returns
    call_returns = call_estimated_future_value - call_current_value
    put_returns = put_estimated_future_value - put_current_value

    # Output the results
    print(f"European Call Option Price current worth: {call_current_value}")
    print(f"European Call Option Price estimated future value: {call_estimated_future_value}")
    print(f"European Call Option Returns: {call_returns}")
    print(f"European Put Option Price current worth: {put_current_value}")
    print(f"European Put Option Price estimated future value: {put_estimated_future_value}")
    print(f"European Put Option Returns: {put_returns}")