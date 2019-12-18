import numpy as np
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
import re

class Uncertainty:
    def __init__(self, s, symbols=None, verbose=True):
        """
        s: the input expression (only right part of '=', eg. for `y=x^2`, s = "x**2" )
        symbols: the variables contained in the expression, if None, the program will find it automatically.
        verbose: output additional informations if True.
        """
        self.s = s
        self.verbose = verbose
        self.forbidden_symbols = dir(sp)
        self.local_dict = {}
        self.subs_dict = {}

        if symbols is not None:
            self.symbols = symbols
        else:
            self.symbols = self.find_symbols()
        
        if self.verbose:
            print("Symbols: ", self.symbols)

        self.build_dict()
        
        if verbose:
            print("Local Dict: ", self.local_dict)

        self.expr = self.get_expr()
        
        if self.verbose:
            print("Expr: ", self.expr)

        self.sigma_expr = self.get_sigma_expr()

        if self.verbose:
            print("Sigma Expr: ", self.sigma_expr)

    def build_dict(self):
        """
        Register variable symbols and their uncertainty symbols (eg. variable "a" will also register "sigma_a")
        """
        for symbol in self.symbols:
            self.local_dict[symbol] = sp.symbols(symbol)
            self.local_dict["sigma_"+symbol] = sp.symbols("sigma_"+symbol)
    
    def find_symbols(self):
        """
        Automatically find symbols in the expression.
        Symbol definition: continued alphas [+ continued digits]
        Special names (sin, cos, pi, e, ...) are forbidden.
        """
        symbols = []
        x = ''
        for i, c in enumerate(self.s):
            if not str.isalnum(c):
                if x != '' and x not in self.forbidden_symbols:
                    symbols.append(x)
                x = ''
            elif str.isalpha(c) or x != '' and str.isdigit(c):
                x += c

        if x != '' and x not in self.forbidden_symbols:
            symbols.append(x)

        return np.unique(symbols)

    def get_expr(self):
        """
        Convert string to sympy expression.
        """
        expr = parse_expr(self.s, local_dict=self.local_dict)
        return expr

    def get_sigma_expr(self):
        """
        Automatically calculate the expression of uncertainty by "square-sum-root" rule.
        """
        sigma_expr = None
        for symbol in self.symbols:
            partial = sp.diff(self.expr, self.local_dict[symbol])
            if sigma_expr is None:
                sigma_expr = (partial * self.local_dict["sigma_"+symbol])**2
            else:
                sigma_expr += (partial * self.local_dict["sigma_"+symbol])**2
        sigma_expr = sp.sqrt(sigma_expr)

        return sigma_expr


    def set_value(self, symbol, values, error=0, scale=1, std=None):
        """
        Assign values to variables(symbols), and calculate substitution dictionary.
        symbol: the target variable
        values: data
            Only accept list of scalar or numpy.ndarray. 
            If the variable is a scalar, pass in with a racket. (eg. g = [9.8] instead of g = 9.8)
            If the variable is a vector, the mean and std will be calculated automatically.
        error: device error
        scale: scale of values
        std: forcely use already known std, instead of calculate it from data.
        """
        if not isinstance(values, np.ndarray):
            values = np.array(values)

        values = scale * values
        error = scale * error
        if std is not None:
            std = scale * std

        if len(values) == 1:
            mean = values[0]
            std = error/np.sqrt(3) if std is None else std
        else:
            mean = np.mean(values)
            std = np.sqrt(np.std(values, ddof=1)**2 + error**2/3) if std is None else std

        self.subs_dict[self.local_dict[symbol]] = mean
        self.subs_dict[self.local_dict["sigma_"+symbol]] = std

        if self.verbose:
            print(f"{symbol}: mean={mean:e}, std={std:e}")

    def evaluate(self):
        """
        Using subs_dict to eval expr and sigma_expr.
        """
        assert len(self.subs_dict) == len(self.local_dict)
        mean, std = self.expr.evalf(subs=self.subs_dict), self.sigma_expr.evalf(subs=self.subs_dict)
        print(f"Final result: mean={mean:e}, std={std:e}")
        return mean, std

def parse(s):
    return np.array([float(x) for x in re.sub("\s+", ",", s.strip()).split(',')])

if __name__ == "__main__":

    un = Uncertainty("m0*g*l**3/(4*lambda*a*h**3)", verbose=True)

    un.set_value('m0', [1000.57,999.51,1001.24,1001.9,1002.08], 0.02, 0.001)
    un.set_value('g', [9.8])
    un.set_value('l', [24.99], 0.01, 0.01)
    un.set_value('lambda', [0.32075,0.31985,0.31975,0.3191,0.31885], 0.001, 0.01)
    un.set_value('a', parse("1.504	1.504	1.506	1.504	1.508	1.510	1.512	1.510	1.508	1.510"), 0.002, 0.01)
    un.set_value('h', parse("1.526	1.534	1.537	1.539	1.533	1.545	1.520	1.542	1.539	1.526")-0.001, 0.004, 0.001)

    print(un.evaluate())