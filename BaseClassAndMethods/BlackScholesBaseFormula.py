from abc import ABC, abstractmethod
from math import log, sqrt, exp
from scipy.stats import norm

class BlackScholeSNoDividend(ABC):
    """The most basic Black-Scholes(1973) Formula that values a European Option on a stock that does not pay dividends before the option's expiry date

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
        assert isinstance(call_put, str) and (call_put in ['c', 'p']), r"call_put should be either 'c' or 'p'"

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
        return self.d1 - (self.sigma * sqrt(self.T))

    @property
    def price(self):
        if self.call_put == 'c':
            return self.S*norm.cdf(self.d1) - self.X*(exp(-self.r*self.T))*norm.cdf(self.d2)
        else:
            return self.X*(exp(-self.r*self.T))*norm.cdf(-self.d2) - self.S*norm.cdf(-self.d1)
    

class BlackScholesMertonWithDividend(BlackScholeSNoDividend):
    """The Merton Extended(1973) Black-Scholes Formula that values a European Option on a stock/index that pays a known dividend yield

    Args:
        S: Spot Price of the Stock/Index
        X: Strike Price of the Option, usually denoted in textbooks as K
        r: Annualized Risk Free Interest Rate, usually the 10-Y Treasury Yield
        q: Annualized Dividend Yield of the stock or index
        T: Time to Maturity in Years
        sigma: Annualized Historical Volatility (Standard Deviation) of the Stock Price 
        call_put: 'c' for call, 'p' for put
    
    Returns:
        Option Price for a dividend paying stock/index if .price is called
        
    """

    def __init__(self, S: float, X: float, r: float, q: float, T: float, sigma: float, call_put: str) -> float:
        assert isinstance(q, (float, int)) and (q>=0), "Dividend Yield - q shld be a non negative float"
        super().__init__(S, X, r, T, sigma, call_put)
        self.q = q

    @property
    def d1(self):
        return (log(self.S/self.X) + (self.r-self.q+(self.sigma**2)/2)*self.T)/(self.sigma*sqrt(self.T))

    ## We don't have to write d2 because it is automatically invoked by the @property d2

    @property
    def price(self):
        if self.call_put == 'c':
            return self.S*exp(-self.q*self.T)*norm.cdf(self.d1) - self.X*(exp(-self.r*self.T))*norm.cdf(self.d2)
        else:
            return self.X*(exp(-self.r*self.T))*norm.cdf(-self.d2) - self.S*exp(-self.q*self.T)*norm.cdf(-self.d1)


class Black76ModelOptionsOnFutures(BlackScholeSNoDividend):
    """The Black(1976) Formula that values a European Option on a Futures/Forward contract with Initial Price F

    Args:
        F: Initial Price of the Futures/Forward Contract
        X: Strike Price of the Option, usually denoted in textbooks as K
        r: Annualized Risk Free Interest Rate, usually the 10-Y Treasury Yield
        T: Time to Maturity in Years
        sigma: Annualized Historical Volatility (Standard Deviation) of the Futures/Forward Contract Price 
        call_put: 'c' for call, 'p' for put
    
    Returns:
        Option Price for a Futures/Forward contract if .price is called
    """

    def __init__(self, F: float, X: float, r: float, T: float, sigma: float, call_put: str) -> float:
        super().__init__(F, X, r, T, sigma, call_put)
        self.F = self.S

    @property
    def d1(self):
        return (log(self.F/self.X) + ((self.sigma**2)/2)*self.T)/(self.sigma*sqrt(self.T))

    ## We don't have to write d2 because it is automatically invoked by the @property d2

    @property
    def price(self):
        if self.call_put == 'c':
            return exp(-self.r*self.T)*(self.F*norm.cdf(self.d1) - self.X*norm.cdf(self.d2))
        else:
            return exp(-self.r*self.T)*(self.F*norm.cdf(-self.d1) - self.X*norm.cdf(-self.d2))*-1


class GeneralBlackScholesMertonFormula(ABC):##for checking?
    """The General Black-Scholes-Merton Formula that can be used for different option pricing settings (see b in Args)

    Args:
        S: Spot Price of the Stock
        X: Strike Price of the Option, usually denoted in textbooks as K
        b: Annualized Cost-of-Carry rate, \
            if b == r :     BlackScholes(1973) stock option model
            elif b == r-q : Merton(1973) Model with dividend
            elif b == 0 :   Black76 Futures options model
            elif b == 0 and r == 0 : Asay(1982) Margined Futures Options Model (e.g. option premium is reinvested)
            elif b == r_local - r_foreign : Garman and Kohlhagen(1983) Currency Option Model
        r: Annualized Risk Free Interest Rate, usually the 10-Y Treasury Yield
        T: Time to Maturity in Years
        sigma: Annualized Historical Volatility (Standard Deviation) of the Stock Price 
        call_put: 'c' for call, 'p' for put
    
    Returns:
        Option Price for a non-dividend paying stock if .price is called
        
    """
    def __init__(self, S: float, X: float, b:float, r: float, T: float, sigma: float, call_put: str) -> float:
        super().__init__()
        assert isinstance(S, (float, int)) and (S>=0), "Spot Price - S shld be a non negative float"
        assert isinstance(X, (float, int)) and (X>=0), "Strike/Exercise Price - X shld be a non negative float"
        assert isinstance(b, (float, int)) and (b>=0), "Cost of Carry Rate - b shld be a non negative float"
        assert isinstance(r, (float, int)), "Risk Free Rate - r shld be a float" # rfr can be negative
        assert isinstance(T, (float, int)) and (T>=0), "Time to maturity - T shld be a non negative float"
        assert isinstance(sigma, (float, int)) and (sigma>=0), "Historical Vol - sigma shld be a non negative float"
        assert isinstance(call_put, str) and (call_put in ['c', 'p']), r"call_put should be either 'c' or 'p'"

        self.S = S
        self.X = X
        self.b = b
        self.r = r
        self.T = T
        self.sigma = sigma
        self.call_put = call_put

    @property
    def df(self): ##discount Factor
        return exp((self.b-self.r)*self.T)
    
    @property
    def d1(self):
        return (log(self.S/self.X) + (self.b+(self.sigma**2)/2)*self.T)/(self.sigma*sqrt(self.T))

    @property
    def d2(self):
        return self.d1 - (self.sigma * sqrt(self.T))

    @property
    def price(self):
        if self.call_put == 'c':
            return self.S*self.df*norm.cdf(self.d1) - self.X*(exp(-self.r*self.T))*norm.cdf(self.d2)
        
        return self.X*(exp(-self.r*self.T))*norm.cdf(-self.d2) - self.S*self.df*norm.cdf(-self.d1)
        
    @property
    def delta(self):
        if self.call_put == 'c':
            return self.df*norm.cdf(self.d1)
        
        return self.df*(norm.cdf(self.d1)-1)
    
    @property
    def gamma(self):
        return (norm.pdf(self.d1) * self.df)/(self.S*self.sigma*sqrt(self.T))
    
    @property
    def vega(self):
        return self.S*self.df*norm.pdf(self.d1)*sqrt(self.T)
    
    @property
    def vanna(self):
        #d(Delta)/d(vol), or, d(Vega)/d(spot)
        return (-self.df*self.d2/self.sigma)*norm.pdf(self.d1)
    
    @property
    def theta(self):
        #expected bleed / time decay
        if self.call_put == 'c':
            return -(self.S*self.df*norm.pdf(self.d1)*self.sigma/(2*sqrt(self.T))) - (self.b-self.r)*self.S*self.df*norm.cdf(self.d1) - self.r*self.X*exp(-self.r*self.T)*norm.cdf(self.d2)
        return -(self.S*self.df*norm.pdf(self.d1)*self.sigma/(2*sqrt(self.T))) + (self.b-self.r)*self.S*self.df*norm.cdf(self.d1) + self.r*self.X*exp(-self.r*self.T)*norm.cdf(self.d2)

    @property
    def rho(self):
        if self.b != 0 : #not futures option(???) == stocks(???) page 25
            if self.call_put == 'c':
                return self.T*self.X*exp(-self.r*self.T)*norm.cdf(self.d2)
            return self.T*self.X*exp(-self.r*self.T)*norm.cdf(-self.d2)
        else: # futures options
            return -self.T*self.price