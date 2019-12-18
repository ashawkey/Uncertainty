# Uncertainty



```python
from uncertainty import Uncertainty, parse

un = Uncertainty("m0*g*l**3/(4*lambda*a*h**3)", verbose=True)

un.set_value('m0', [1000.57,999.51,1001.24,1001.9,1002.08], 0.02, 0.001)
un.set_value('g', [9.8])
un.set_value('l', [24.99], 0.01, 0.01)
un.set_value('lambda', [0.32075,0.31985,0.31975,0.3191,0.31885], 0.001, 0.01)
un.set_value('a', parse("1.504	1.504	1.506	1.504	1.508	1.510	1.512	1.510	1.508	1.510"), 0.002, 0.01)
un.set_value('h', parse("1.526	1.534	1.537	1.539	1.533	1.545	1.520	1.542	1.539	1.526")-0.001, 0.004, 0.001)

print(un.evaluate())
```