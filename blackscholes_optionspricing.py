import math
from scipy.stats import norm

# Define the Contract structure
class Contract:
    def __init__(self):
        self.premium = 0.0
        self.dte = 0
        self.delta = 0.0
        self.gamma = 0.0
        self.theta = 0.0
        self.vega = 0.0
        self.rho = 0.0
        self.implied_volatility = 0.0
        self.intrinsic_value = 0.0

# Black-Scholes option pricing model
def blackScholesOptionPricing(S0, K, r, sigma, T, isCallOption):
    con = Contract()
    days_till_expiry = int(T * 365.2425)
    con.dte = days_till_expiry

    d1 = (math.log(S0/K) + (r + 0.5 * sigma**2) * T) / (sigma * math.sqrt(T))
    d2 = d1 - sigma * math.sqrt(T)

    if isCallOption:
        # Call option price
        con.premium = S0 * norm.cdf(d1) - K * math.exp(-r * T) * norm.cdf(d2)
        # Greeks
        con.delta = norm.cdf(d1)
        con.gamma = norm.pdf(d1) / (S0 * sigma * math.sqrt(T))
        con.theta = -(S0 * norm.pdf(d1) * sigma) / (2 * math.sqrt(T)) - r * K * math.exp(-r * T) * norm.cdf(d2)
        con.vega = S0 * norm.pdf(d1) * math.sqrt(T)
        con.rho = K * T * math.exp(-r * T) * norm.cdf(d2)
        con.implied_volatility = sigma - ((con.premium - (con.premium - 0.01)) / con.vega)
        con.intrinsic_value = max(S0 - K, 0.0)

    else:
        # Put option price
        con.premium = K * math.exp(-r * T) * norm.cdf(-d2) - S0 * norm.cdf(-d1)
        # Greeks
        con.delta = norm.cdf(d1) - 1
        con.gamma = norm.pdf(d1) / (S0 * sigma * math.sqrt(T))
        con.theta = -(S0 * norm.pdf(d1) * sigma) / (2 * math.sqrt(T)) + r * K * math.exp(-r * T) * norm.cdf(-d2)
        con.vega = S0 * norm.pdf(d1) * math.sqrt(T)
        con.rho = -K * T * math.exp(-r * T) * norm.cdf(-d2)
        con.implied_volatility = sigma - (((con.premium - 0.01) - con.premium) / con.vega)
        con.intrinsic_value = max(K - S0, 0.0)

    return con

if __name__ == "__main__":
    # Option parameters
    S0 = 100.0   # Initial stock price
    K = 110.0    # Strike price
    r = 0.05     # Risk-free rate
    sigma = 0.2  # Volatility
    T = 1        # Time to maturity (in years)
    
    # Calculate option prices
callContract = blackScholesOptionPricing(S0, K, r, sigma, T, True)
putContract = blackScholesOptionPricing(S0, K, r, sigma, T, False)

# Output the results
print("European Call Option:")
print(f"  Price: {callContract.premium}")
print(f"  Days to Expiry (DTE): {callContract.dte}")
print(f"  Delta: {callContract.delta}")
print(f"  Gamma: {callContract.gamma}")
print(f"  Theta: {callContract.theta}")
print(f"  Vega: {callContract.vega}")
print(f"  Rho: {callContract.rho}")
print(f"  Implied Volatility: {callContract.implied_volatility}")
print(f"  Intrinsic Value: {callContract.intrinsic_value}")

print("\nEuropean Put Option:")
print(f"  Price: {putContract.premium}")
print(f"  Days to Expiry (DTE): {putContract.dte}")
print(f"  Delta: {putContract.delta}")
print(f"  Gamma: {putContract.gamma}")
print(f"  Theta: {putContract.theta}")
print(f"  Vega: {putContract.vega}")
print(f"  Rho: {putContract.rho}")
print(f"  Implied Volatility: {putContract.implied_volatility}")
print(f"  Intrinsic Value: {putContract.intrinsic_value}") 
