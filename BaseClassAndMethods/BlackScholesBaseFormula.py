from abc import ABC, abstractmethod
from math import log, sqrt, exp
from scipy.stats.norm import cdf

class BlackScholeNonDividend(ABC):
    """The most basic Black-Scholes Formula that values a European Option on a stock that does not pay dividends before the option's expiry date

    Args:
        S: Spot Price of the Stock
        X: Strike Price of the Option, usually denoted in textbooks as K
        r: Annualized Risk Free Interest Rate, usually the 10-Y Treasury Yield
        T: Time to Maturity in Years
        sigma: Annualized Historical Volatility (Standard Deviation) of the Stock Price 
        call_put: 'c' for call, 'p' for put
    
    Returns:
        Option Price for a non-dividend paying stock if .price is called
        
    """

    def __init__(self, S: float, X: float, r: float, T: float, sigma: float, call_put: str) -> float:
        super().__init__()
        assert isinstance(S, (float, int)) and (S>=0), "Spot Price - S shld be a non negative float"
        assert isinstance(X, (float, int)) and (X>=0), "Strike/Exercise Price - X shld be a non negative float"
        assert isinstance(r, (float, int)), "Risk Free Rate - r shld be a float" # rfr can be negative
        assert isinstance(T, (float, int)) and (T>=0), "Time to maturity - T shld be a non negative float"
        assert isinstance(sigma, (float, int)) and (sigma>=0), "Historical Vol - sigma shld be a non negative float"
        assert isinstance(call_put, str) and (sigma in ['c', 'p']), r"call_put should be either 'c' or 'p'"

        self.S = S
        self.X = X
        self.r = r
        self.T = T
        self.sigma = sigma
        self.call_put = call_put
    
    @property
    def d1(self):
        return (log(self.S/self.X) + (self.r+(self.sigma**2)/2)*self.T)/(self.sigma*sqrt(self.T))

    @property
    def d2(self):
        return self.d1 - (self.sigma * self.T)

    @property
    def price(self):
        if self.call_put == 'c':
            return self.S*cdf(self.d1) - self.X*(exp**(-self.r*self.T))*cdf(self.d2)
        else:
            return self.X*(exp**(-self.r*self.T))*cdf(-self.d2) - self.S*cdf(-self.d1)
    


