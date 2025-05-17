
"""

This is a the result of a university numerical methods class.

It contains methods numerically solving
	- Solving for x in an equation (roots)
	- Solving for X0 and xn in a system of equations
	- Interpolation
	- Extrapolation
	- Differentiation
	- Integration
	- Differential equations

"""


from terminalTrueColors import *
import numpy as np
from sympy import *
from spb import *
import itertools
from itertools import permutations
from statistics import *
import numpy as np 

# x is the only necesary one, y and z are only for 
x, z, y = symbols('x z y')

#used for systems of equations
x1, x2, x3, x4, x5 = symbols('x1 x2 x3 x4 x5')
xSymbols = symbols('x1 x2 x3 x4 x5')





#--------------------------------------------some printer helper functions-----------------------------
def printHeader(row, size=10):
	output = '|'
	for cell in row:
		if cell == 'i'or cell == ' ':
			output += cadd(colors['orange'],str(cell).center(3,"-")) + "|"
		else:
			output += cadd(colors['orange'],str(cell).center(size,"-")) + "|"
	print(output)
def printRow(row, size=10):
	output = '|'
	for cell in row:
		if type(cell)==str:
			output += cell.center(3," ") + "|"
		elif type(cell)==int:
			output += cadd(colors['orange'], str(cell).center(3," ")) + "|"
		else:
			if cell < 0.00004 and cell>0:
				output += cadd(rgb(9, 164, 81),str("{:.6f}".format(N(cell))).ljust(size," ")) + "|"
			else:
				output += str("{:.6f}".format(N(cell))).ljust(size," ") + "|"
			
	print(output)

# --------------------------------------------------------Numerical solver-------------------
class NumSolver:
	def __init__(self,f,row,maxError=0.00004, maxIter=100, lastRow=[0], solution= False):
		
		self.f = f
		self.f_diff = f.diff(x)
		cprint(colors["orange"], "Ecuation: ", end='')
		print(self.f,'= 0')

		if type(self) == Secant: 
			cprint(colors["orange"], 'X_0: ', end='')
			print(lastRow)
			cprint(colors["orange"], ' X_1: ', end='')
			print(row)
		else:
			cprint(colors["orange"], 'X_0: ', end='')
			print(row)
		
		if solution:
			if type(self) == FixedPoint:
				self.trueValue = solve(f-x,x)
			else:
				self.trueValue = solve(f,x)
			cprint(colors["orange"], "Solutions: ", end='')
			print(self.trueValue)
		else:
			self.trueValue=[9999]

		self.lastRow = lastRow
		self.row = row
		self.maxError = maxError
		self.maxIter = maxIter
		self.stopping = False
		self.whatTrueValue = None
		self.printTable()
		

	def getWhatTrueValue(self, xr):
		if not self.whatTrueValue:
			dif = 99999
			for i, value in enumerate(self.trueValue):
				vt = self.trueValue[i]
				possibleDif=abs(xr-vt)
				if dif > possibleDif:
					dif = possibleDif
					self.whatTrueValue = i

	def getError(self, i):
		
		xr = self.row[i]
		last_xr = self.lastRow[i]
		#self.getWhatTrueValue(xr)
		#vt = self.trueValue[self.whatTrueValue]

		#stopping criteria
		ea = abs((xr-last_xr)/xr)*100
		#et = abs((vt-xr)/vt)*100
		if type(self) == FixedPoint:
			ea = abs(xr-last_xr)
			#et = abs(vt-xr)

		self.stopping = ea < self.maxError
		return [ea]

	def printTable(self):
		i = 1
		while not self.stopping:
			
			if type(self) == Bisection or type(self) == FalsePosition:
				header = ["i"," "]+["xi","xr","xu","εa"]
			else:
				header = ["i"," "]+["xi","εa"]
			if i==1:
				printHeader(header)
			row = [i, " "] + self.getRow(i) + self.getError()
			printRow(row)
			if i > 50:
				self.stopping = True
			i+=1

class FalsePosition(NumSolver):
	def __init__(self, *args, **kwargs):
		print()
		cprint(colors["orange"], "Numerical Solver | ", end='')
		cprint(colors["yellow"],'False position')
		super().__init__(*args,**kwargs)
	def getError(self):
			return super().getError(1)

	def getRow(self, i):
		f = self.f
		xi, xr, xu = self.row
		self.lastRow = self.row.copy()

		xr = N(xu - (f.subs(x,xu)*(xi-xu))/(f.subs(x,xi)-f.subs(x,xu)))
		if N(f.subs(x,xi)*f.subs(x,xr))<0:
			xu = xr
		else:
			xi = xr

		self.row = [xi, xr, xu]
		return self.row

class Bisection(NumSolver):
	def __init__(self, *args, **kwargs):
		print()
		cprint(colors["orange"], "Numerical Solver | ", end='')
		cprint(colors["lightGreen"],'Bisection')
		super().__init__(*args,**kwargs)
	def getError(self):
			return super().getError(1)

	def getRow(self, i):
		f = self.f
		xi, xr, xu = self.row
		self.lastRow = self.row.copy()

		
		xr = N((xi + xu)/2)
		if N(f.subs(x, xi)*f.subs(x, xr))<0:
			xu = xr
		else:
			xi = xr

		self.row = [xi, xr, xu]
		return self.row

class NewtonRaphson(NumSolver):
	def __init__(self, *args, **kwargs):
		print()
		cprint(colors["orange"], "Numerical Solver | ", end='')
		cprint(rgb(9, 164, 81),'Newton-Raphson')
		super().__init__(*args,**kwargs)
	def getError(self):
			return super().getError(0)

	def getRow(self, i):
		f = self.f
		f_diff = self.f_diff
		xi = self.row[0]
		self.lastRow = self.row.copy()

		xi = N(xi-f.subs(x, xi)/f_diff.subs(x, xi))
		
		self.row = [xi]
		return self.row

class FixedPoint(NumSolver):
	# clear x, before inputing equation
	def __init__(self, *args, **kwargs):
		print()
		cprint(colors["orange"], "Numerical Solver | ", end='')
		cprint(colors["darkGreen"],'Fixed point')
		super().__init__(*args,**kwargs)
	def getError(self):
			return super().getError(0)

	def getRow(self, i):
		f = self.f
		xi = self.row[0]
		self.lastRow = self.row.copy()
		last_xi = self.lastRow[0]

		xi = N(f.subs(x, last_xi))
		
		self.row = [xi]
		return self.row

class Secant(NumSolver):
	def __init__(self, *args, **kwargs):
		print()
		cprint(colors["orange"], "Numerical Solver | ", end='')
		cprint(colors["darkBlue"],'Secant')
		super().__init__(*args,**kwargs)
	def getError(self):
			return super().getError(0)

	def getRow(self, i):
		f = self.f
		xi = self.row[0]
		last_xi = self.lastRow[0]
		self.lastRow = self.row.copy()


		xi = N(xi - (f.subs(x, xi)*(xi-last_xi))/(f.subs(x, xi)-f.subs(x, last_xi)))
		
		self.row = [xi]
		return self.row


# --------------------------------------Numerical Numsolvers for system of equations-------------------


class NumeSolverSE:
	def __init__(self,f,row,maxError=0.00004, maxIter=100, lastRow=[0], solution= False):
		self.f = f
		self.F = Matrix(f)
		self.f_diff = self.F.jacobian(Matrix(list(self.F.free_symbols)))
		cprint(colors["orange"], "Ecuation: ", end='')
		print([str(i)+' = 0 \n' for i in self.F])

		cprint(colors["orange"], "X_0: ", end='')
		print(row)
		
		self.lastRow = lastRow
		self.row = Matrix(row)
		self.maxError = maxError
		self.maxIter = maxIter
		self.stopping = False
		self.whatTrueValue = None
		self.printTable()
		

	def getError(self, i):
		xr = self.row
		last_xr = self.lastRow
		ea = ((last_xr[0]-xr[0])**2 + (last_xr[1]-xr[1])**2)**(1/2)
		et = 100.0
		
		self.stopping = ea < self.maxError or ea > 500
		return [ea]
		

	def printTable(self):
		i = 1
		while not self.stopping:
			header = ["i"," "]+["x1","x2"]+["εa"]
			if i==1:
				printHeader(header)
			row = [i, " "] + list(self.getRow(i)) + self.getError()
			printRow(row)
			if i > 50:
				self.stopping = True
			i+=1



class FixedPointSE(NumeSolverSE):
	# clear x and y, before inputing equation
	def __init__(self,f , *args, **kwargs):
		print()
		cprint(colors["orange"], "Numerical Solver | ", end='')
		cprint(colors["lightBlue"],'FixedPointLE')

		super().__init__(f, *args,**kwargs)
	def getError(self):
			return super().getError(0)

	def getRow(self, i):
		f = self.f
		self.lastRow = self.row.copy()
		xi = self.row
		last_xi = self.lastRow
		
		subsValues = {xSymbols[i]: last_xi[i] for i in range(len(xi))}
		xi[0] = [N(f[0].subs(subsValues))]
		subsValues = {xSymbols[i]: xi[i] for i in range(len(xi))}
		xi[1] = [N(f[1].subs(subsValues))]
		self.row = xi
		
		return self.row

class NewtonRaphsonSE(NumeSolverSE):
	def __init__(self, *args, **kwargs):
		print()
		cprint(colors["orange"], "Numerical Solver | ", end='')
		cprint(colors["blue"],'Newton-Raphson LE')
		super().__init__(*args,**kwargs)
	def getError(self):
			return super().getError(0)

	def getRow(self, i):
		F = self.F
		J = self.f_diff
		xi = self.row
		self.lastRow = self.row.copy()
		
		subsValues = {xSymbols[i]: xi[i] for i in range(len(xi))}
		mat = []
		JinvSubs = J.inv().subs(subsValues)
		xi = N(xi-(JinvSubs*F.subs(subsValues)))
		self.row = xi
		return self.row

# --------------------------------------------------------Interpolation-------------------
class Interpolation:
	#xydata -> eq
	def subs(self, xValue):
		print(f"f({xValue}) =",self.eq.subs(x,xValue))
		return N(self.eq.subs(x,xValue))

class Lagrange(Interpolation):
	def __init__(self,xData,yData):
		print()
		cprint(colors["orange"], "Interpolation | ", end='')
		cprint(colors["pink"],'Lagrange')

		self.x = xData
		self.y = yData

		temp2 = 0
		for i in range(len(xData)):
			temp = 1
			for j in range(len(xData)):
				if j != i:
					temp *= (x-self.x[j])/(self.x[i]-self.x[j])
			temp2 += temp * self.y[i]

		self.eq = simplify(temp2)
		print("f(x) =",self.eq)


class Newton(Interpolation):

	def b(self,i,n):
		if n == 0:
			return self.y[i]
		else:
			return ( self.b(i, n-1)-self.b(i-1, n-1) ) / (self.x[i] - self.x[i-n])

	def n(self, i):
		mult = 1
		for j in range(i):
			mult *= (x-self.x[j])
		return mult

	def __init__(self,xData,yData):
		
		print()
		cprint(colors["orange"], "Interpolation | ", end='')
		cprint(colors["purple"],'Newton')
		self.x = xData
		self.y = yData

		temp = 0
		for j in range(len(self.x)):
			#print('b',j,':',self.b(j,j))
			#print('n',j,':',self.b(j,j))
			temp += self.n(j)*self.b(j,j)
		
		self.eq = simplify(temp)
		print("f(x) =",expand(self.eq))


# --------------------------------------------------------Regression-------------------


class leastSquares:
	def subs(self, xValue):
		print(f"f({xValue}) =",self.eq.subs(x,xValue))
		return float(self.eq.subs(x,xValue))
	def __init__(self,xData,yData, model='linear'):
		
		#Least Squares
		#data: list
		#f: equation
		print()
		cprint(colors["orange"], "Regression | ", end='')
		cprint(colors["white"],'Least Squares')

		
		if model == 'exp':
			yData = [ln(i) for i in yData]
		elif model == 'pow':
			xData = [log(i,10) for i in xData]
			yData = [log(i,10) for i in yData]
		elif model == 'base10':
			yData = [log(i,10) for i in yData]
		elif model == 'frac':
			xData = [1/i for i in xData]
			yData = [1/i for i in yData]


		n = len(yData)
		xSum = sum(xData)
		ySum = sum(yData)
		xySum = sum(x*y for x,y in zip(xData, yData))
		xxSum = sum(x**2 for x in xData)
		yySum = sum(y**2 for y in yData)

		numerator = n*xySum - xSum*ySum
		denominator = (n*xxSum - xSum**2)
		b1 = numerator / denominator
		b0 = (ySum - b1*xSum)/n
		
		#r = [n(Σxy) - (Σx)(Σy)] / sqrt([nΣx² - (Σx)²][nΣy² - (Σy)²])
		r = numerator / ((n*xxSum - (xSum**2)) * (n*yySum - (ySum**2)))**0.5
		

		SSE = sum((y - (b0 + b1*x))**2 for x,y in zip(xData, yData))
		stde = (SSE/(n-2))**0.5
		if model == 'linear':
			print('alpha:',N(b0))
			print('beta:',N(b1))
			eq = b0 + b1 * x
		elif model == 'exp':
			print('alpha:',N(b0))
			print('beta:',N(b1))
			eq = exp(b0)*exp(b1 * x)
		elif model == 'base10':
			print('alpha:',N(b0))
			print('beta:',N(b1))
			eq = (10**b0)*(10**(b1 * x))
		elif model == 'pow':
			print('b0:',N(b0))
			print('b1:',N(b1))
			eq = N((10**b0))*x**N(b1)
		elif model == 'frac':
			print('b0:',N(b0))
			print('b1:',N(b1))
			eq = (1/b0)*(x/(b1/b0+x))

		print("r =",N(r))
		print("r^2 =",N(r**2))
		print("stde =",N(stde))

		print("f(x) =",eq)
		self.r = N(r)
		self.stde = N(stde)
		self.eq = eq
		




# --------------------------------------------------------Statistics-------------------
def statisticsInfo(data):
	print()
	cprint(colors["orange"], "Statistics | ", end='')
	cprint(colors["white"],'simple data')

	cprint(colors["orange"], "Statistics information", end='')
	cprint(colors["white"],'Least Squares')

	cprint(colors["orange"], "Data:", end='')
	cprint(colors["white"],data)

	cprint(colors["orange"], "Mean:", end='')
	cprint(colors["white"],mean(data))

	cprint(colors["orange"], "Median:", end='')
	cprint(colors["white"],median(data))

	cprint(colors["orange"], "Modes:", end='')
	cprint(colors["white"],multimode(data))

	cprint(colors["orange"], "Range:", end='')
	cprint(colors["white"],max(data)-min(data))

	cprint(colors["orange"], "Standard deviation:", end='')
	cprint(colors["white"],stdev(data))

	cprint(colors["orange"], "Variance:", end='')
	cprint(colors["white"],variance(data))

	cprint(colors["orange"], "Variation coefficient:", end='')
	cprint(colors["white"],pstdev(data)/mean(data)*100)



# --------------------------------------------------------Differentiation-------------------
def unequalDifferential(xData, yData, xValue):
	#return the aproximations using three points
	i=1
	return (yData[i-1]*((2*x-xData[i]-xData[i+1])/((xData[i-1]-xData[i])*(xData[i-1]-xData[i+1]))) + yData[i]*((2*x-xData[i-1]-xData[i+1])/((xData[i]-xData[i-1])*(xData[i]-xData[i+1])))).subs(x,xValue)

class numericalDerivative:
	def __init__(self, xData=[], yData=[], degree=1, quality=0):

		print()
		cprint(colors["orange"], "Derivation | ", end='')
		cprint(colors["pink"],f'Equally spaced data table')
		cprint(colors["orange"], "Degree:", end='')
		cprint(colors["white"],degree)
		cprint(colors["orange"], 'Quality:', end='')
		cprint(colors["white"],quality)
		
		self.degree = degree
		self.xData = xData 
		self.yData = yData 
		self.quality = quality
		self.printTable()

	def printTable(self):
		for i,temp1 in enumerate(zip(self.xData,self.yData)):
			xData,yData = temp1
			
			header = ["i"," "]+["x","y(x)"," "]+["diff","Debug info"]
			if i==0:
				printHeader(header)

			temp=0
			if self.degree%2 == 0:
				temp = 1
			
			if i-self.quality-temp <= 0:
				tipo = 'forward'
			else:
				tipo = 'center'

			if self.degree == 1:
				temp = 1

			if i+1+self.quality+self.degree-temp >= len(self.xData):
				tipo = 'backward'
			row = [i, " "] + [N(xData), N(yData), " "] + [self.getRow(i),tipo]
			printRow(row)
			i+=1

	def getRow(self, i):
		xData = self.xData
		yData = self.yData
		degree = self.degree
		quality = self.quality

		# this is for equally spaced data
		temp=0
		if degree%2 == 0:
			temp = 1
		
		if i-quality-temp <= 0:
			tipo = 'forward'
		else:
			tipo = 'center'

		if degree == 1:
			temp = 1

		if i+1+quality+degree-temp >= len(xData):
			tipo = 'backward'
		
		if tipo == 'forward':
			h = xData[i+1] - xData[i] 
			if degree == 1:
				if not quality:
					return (yData[i+1] - yData[i])/h
				else:
					return (-yData[i+2] + 4*yData[i+1] - 3*yData[i])/ (2*h)
			elif degree == 2:
				if not quality:
					return (yData[i+2] - 2*yData[i+1] + yData[i])/ (h**2)
				else:
					return (- yData[i+3] + 4* yData[i+2] - 5*yData[i+1] + 2*yData[i])/ (h**2)
			elif degree == 3:
				if not quality:
					return (yData[i+3] - 3* yData[i+2] + 3*yData[i+1] - yData[i])/ (h**3)
				else:
					return (-3*yData[i+4], + 14*yData[i+3] - 24* yData[i+2] + 18*yData[i+1] - 5*yData[i])/ (2*h**3)
		elif tipo == 'center':
			h = xData[i+1] - xData[i] 
			if degree == 1:
				if not quality:
					return (yData[i+1] - yData[i-1])/(2*h)
				else:
					return (-yData[i+2] + 8*yData[i+1] - 8*yData[i-1] +yData[i-2])/ (12*h)
			elif degree == 2:
				if not quality:
					return (yData[i+1] - 2*yData[i] + yData[i-1])/ (h**2)
				else:
					return (- yData[i+2] + 16* yData[i+1] - 30*yData[i] + 16*yData[i-1] -yData[i-2])/ (12*h**2)
			elif degree == 3:
				if not quality:
					return (yData[i+2] - 2* yData[i+1] + 2*yData[i-1] - yData[i-2])/ (2*h**3)
				else:
					return (-yData[i+3], + 8*yData[i+2] - 13* yData[i+1] + 13*yData[i-1] - 8*yData[i-2]+ yData[i-3])/ (8*h**3)
		elif tipo == 'backward':
			h = xData[i] - xData[i-1] 
			if degree == 1:
				if not quality:
					return (yData[i] - yData[i-1])/h
				else:
					return (3*yData[i] - 4*yData[i-1] + yData[i-2])/ (2*h)
			elif degree == 2:
				if not quality:
					return (yData[i] - 2*yData[i-1] + yData[i-2])/ (h**2)
				else:
					return (2*yData[i] - 5* yData[i-1] + 4*yData[i-2] - yData[i-3])/ (h**2)
			elif degree == 3:
				if not quality:
					return (yData[i] - 3* yData[i-1] + 3*yData[i-2] - yData[i-3])/ (h**3)
				else:
					return (5*yData[i], -18*yData[i-1] + 24* yData[i-2] -14*yData[i-3] +3*yData[i-4])/ (2*h**3)


class unequalNumericalDerivative:
	def __init__(self, xData=[], yData=[], degree=1, quality=0):

		print()
		cprint(colors["orange"], "Derivation | ", end='')
		cprint(colors["pink"],f'Unequally spaced data table')
		
		self.degree = degree
		self.xData = xData 
		self.yData = yData 
		self.quality = quality
		self.diffData=[]
		self.printTable()

	def printTable(self):
		self.diffData=[]
		for i,temp1 in enumerate(zip(self.xData,self.yData)):
			xValue,yValue = temp1
			
			header = ["i"," "]+["x","y(x)"," "]+["diff"]
			if i==0:
				printHeader(header)

			temp=0
			if self.degree%2 == 0:
				temp = 1
			
			if i-self.quality-temp <= 0:
				tipo = 'forward'
			else:
				tipo = 'center'

			if self.degree == 1:
				temp = 1

			if i+1+self.quality+self.degree-temp >= len(self.xData):
				tipo = 'backward'
			temp = self.getRow(i,xValue)
			self.diffData.append(temp)
			row = [i, " "] + [N(xValue), N(yValue), " "] + [temp]
			printRow(row)


	def getRow(self, i,xValue):
		xData = self.xData
		yData = self.yData
		degree = self.degree
		quality = self.quality

		i = xData.index(xValue)
		
		if i == 0:
			i += 1
		if i == len(xData)-1:
			i-=1
		
		return N((
			yData[i-1]*((2*x-xData[i]-xData[i+1])/((xData[i-1]-xData[i])*(xData[i-1]-xData[i+1]))) + 
			yData[i]*((2*x-xData[i-1]-xData[i+1])/((xData[i]-xData[i-1])*(xData[i]-xData[i+1]))) +
			yData[i+1]*((2*x-xData[i-1]-xData[i])/((xData[i+1]-xData[i-1])*(xData[i+1]-xData[i])))
			).subs(x,xValue))
	
def simpleNumericalDerivative(eq, h, xValue):

	print()
	cprint(colors["orange"], "Derivation | ", end='')
	cprint(colors["pink"],f'All posible derivatives')
	
	tValue = eq.diff().subs(x,xValue)

	print(cadd(colors['green'],'trueValue:'),N(tValue))
	
	xData = list(np.arange(xValue-2*h,xValue+3*h,h))
	yData = [eq.subs(x,xValue) for xValue in xData]
	print(cadd(colors['blue'],'xData:'),xData)
	print(cadd(colors['blue'],'yData:'),yData)
	
	i = 2
	print()

	temp = N((yData[i+1] - yData[i])/h)
	print(f"{cadd(colors['orange'],'Forward')} O(h):", temp)
	print(cadd(colors['red'],'εt:'),N(abs(temp-tValue)/tValue*100),'%')

	temp = N((-yData[i+2] + 4*yData[i+1] - 3*yData[i])/ (2*h))
	print(f"{cadd(colors['orange'],'Forward ')}O(h^2):", temp)
	print(cadd(colors['red'],'εt:'),N(abs(temp-tValue)/tValue*100),'%')
	print()

	temp = N((yData[i] - yData[i-1])/h)
	print(f"{cadd(colors['orange'],'Backward')} O(h):", temp)
	print(cadd(colors['red'],'εt:'),N(abs(temp-tValue)/tValue*100),'%')

	temp = N((3*yData[i] - 4*yData[i-1] + yData[i-2])/ (2*h))
	print(f"{cadd(colors['orange'],'Backward ')}O(h^2):", temp)
	print(cadd(colors['red'],'εt:'),N(abs(temp-tValue)/tValue*100),'%')

	print()

	temp = N((yData[i+1] - yData[i-1])/(2*h))
	print(f"{cadd(colors['orange'],'Center ')}O(h^2):", temp)
	print(cadd(colors['red'],'εt:'),N(abs(temp-tValue)/tValue*100),'%')

	temp = N((-yData[i+2] + 8*yData[i+1] - 8*yData[i-1] +yData[i-2])/ (12*h))
	print(f"{cadd(colors['orange'],'Center ')}O(h^4):", temp)
	print(cadd(colors['red'],'εt:'),N(abs(temp-tValue)/tValue*100),'%')

	#numericalDifferential(xData=[], yData=[])


# --------------------------------------------------------Differential equations-------------------



from sympy import *

class RungeKutta:
	def __init__(self, xData=[], y0=0, eq=None, method="euler"):
		self.xData = xData
		self.y0 = y0
		self.eq = eq
		self.method = method
		self.yData = []

		print()
		methodName = {
			"euler": "First-Order Runge-Kutta (Euler's Method)",
			"heun": "Second-Order Runge-Kutta (Heun's Method)",
			"ralston": "Second-Order Runge-Kutta (Ralston's Method)",
			"midpoint": "Second-Order Runge-Kutta (Midpoint Method)",
			"rk4": "Fourth-Order Runge-Kutta (Classic RK4)"
		}.get(self.method, "Unknown Method")

		cprint(colors["orange"], "Integration | ", end='')
		cprint(colors["pink"], methodName)

		cprint(colors["orange"], "Ecuation: ", end='')
		print("dy/dx =",self.eq)

		self.printTable()

	def printTable(self):
		x, y = symbols('x y')
		yCurrent = self.y0
		self.yData = [yCurrent]

		printHeader(["i", " ", "x", "y_i","dy"])

		for i in range(len(self.xData)):
			xi = self.xData[i]
			if i > 0:
				h = xi - self.xData[i - 1]
				xiPrev = self.xData[i - 1]
				yPrev = self.yData[-1]

				if self.method == "euler":
					k1 = self.eq.subs({x: xiPrev, y: yPrev})
					yCurrent = yPrev + h * k1

				elif self.method == "heun":
					k1 = self.eq.subs({x: xiPrev, y: yPrev})
					k2 = self.eq.subs({x: xiPrev + h, y: yPrev + h * k1})
					yCurrent = yPrev + h * (k1 + k2) / 2

				elif self.method == "ralston":
					k1 = self.eq.subs({x: xiPrev, y: yPrev})
					k2 = self.eq.subs({x: xiPrev + 3*h/4, y: yPrev + 3*h/4 * k1})
					yCurrent = yPrev + h * (1/3 * k1 + 2/3 * k2)

				elif self.method == "midpoint":
					k1 = self.eq.subs({x: xiPrev, y: yPrev})
					k2 = self.eq.subs({x: xiPrev + h/2, y: yPrev + h/2 * k1})
					yCurrent = yPrev + h * k2

				elif self.method == "rk4":
					k1 = self.eq.subs({x: xiPrev, y: yPrev})
					k2 = self.eq.subs({x: xiPrev + h/2, y: yPrev + h/2 * k1})
					k3 = self.eq.subs({x: xiPrev + h/2, y: yPrev + h/2 * k2})
					k4 = self.eq.subs({x: xiPrev + h, y: yPrev + h * k3})
					yCurrent = yPrev + h * (k1 + 2*k2 + 2*k3 + k4) / 6

				self.yData.append(yCurrent)

			dy = self.eq.subs({x: xi, y: yCurrent})
			printRow([i, " ", xi, yCurrent,dy])






# --------------------------------------------------------Integration-------------------


# list(np.arange(xValue-2*h,xValue+3*h,h))
def rectangleRule(eq, xData):

	print()
	cprint(colors["orange"], "Numerical Integration | ", end='')
	cprint(colors["green"],'Rectangule rule')

	cprint(colors['orange'],'Equation: ',end='')
	print(eq)
	cprint(colors['orange'],'xData: ',end='')
	print(xData)
		
	a = xData[0]
	n = len(xData)-1
	b = xData[-1]

	h = (b-a)/n
	val = N(sum([eq.subs(x,(xData[i]+xData[i+1])/2)*h for i in range(0,n)]))

	print('Estimated value:',val)
	tval = N(integrate(eq, (x, a, b)))
	print('True value:',tval)
	print('True error %:',abs(tval-val)/tval*100)

	return val

	


def trapezoidalRule(eq, xData):

	print()
	cprint(colors["orange"], "Numerical Integration | ", end='')
	cprint(colors["lightBlue"],'Trapezoidal rule')

	cprint(colors['orange'],'Equation: ',end='')
	print(eq)
	cprint(colors['orange'],'xData: ',end='')
	print(xData)
		
	a = xData[0]
	n = len(xData)-1
	b = xData[-1]

	h = (b-a)/n


	simpleSum = sum(eq.subs(x, xData[i]) for i in range(1, n))

	val = N((h/2) * (eq.subs(x, a) + 2 * simpleSum + eq.subs(x, b)))

	print('Estimated value:',val)
	tval = N(integrate(eq, (x, a, b)))
	print('True value:',tval)
	print('True error %:',abs(tval-val)/tval*100)

	return val

def simpsons13Rule(eq, xData):

	print()
	cprint(colors["orange"], "Numerical Integration | ", end='')
	cprint(colors["blue"],'Simpson 1/3 rule')

	cprint(colors['orange'],'Equation: ',end='')
	print(eq)
	cprint(colors['orange'],'xData: ',end='')
	print(xData)
		
	a = xData[0]
	n = len(xData)-1
	b = xData[-1]

	h = (b-a)/n

	if n % 2 != 0:
		cprint(colors['red'],'Error: ',end='')

		print("Requires an even number of intervals, use Simpson's 3/8 rule instead")
	
	sumOdd = sum(eq.subs(x, xData[i]) for i in range(1, n, 2))
	sumEven = sum(eq.subs(x, xData[i]) for i in range(2, n-1, 2))

	val = N((h/3) * (eq.subs(x, a) + 4*sumOdd + 2*sumEven + eq.subs(x, b)))

	print('Estimated value:',val)
	tval = N(integrate(eq, (x, a, b)))
	print('True value:',tval)
	print('True error %:',abs(tval-val)/tval*100)

	return val

def simpsons38Rule(eq, xData):

	print()
	cprint(colors["orange"], "Numerical Integration | ", end='')
	cprint(colors["darkGreen"],'Simpson 3/8 rule')

	cprint(colors['orange'],'Equation: ',end='')
	print(eq)
	cprint(colors['orange'],'xData: ',end='')
	print(xData)
		
	a = xData[0]
	n = len(xData)-1
	b = xData[-1]

	h = (b-a)/n

	if n % 3 != 0:
		cprint(colors['red'],'Error: ',end='')

		print("Requires an 3divisible number of intervals, use Simpson's 1/3 rule instead")
	
	sum3k = sum(eq.subs(x, xData[i]) for i in range(3, n-2, 3))
	sumOther = sum(eq.subs(x, xData[i]) for i in range(1, n) if i % 3 != 0)

	val = N((3*h/8) * (eq.subs(x, a) + 3*sumOther + 2*sum3k + eq.subs(x, b)))

	print('Estimated value:',val)
	tval = N(integrate(eq, (x, a, b)))
	print('True value:',tval)
	print('True error %:',abs(tval-val)/tval*100)

	return val

def gaussQuadrature(eq, xData, points):
	oldEq = eq

	print()
	cprint(colors["orange"], "Numerical Integration | ", end='')
	cprint(colors["pink"],'Gaussian quadrature')

	a = xData[0]
	n = len(xData)-1
	b = xData[-1]

	h = (b-a)/n

	if a < -1 or b > 1:
		varChange = ((b+a)+(b-a)*x)/2
		eq = eq.subs(x,varChange)
		eq *= diff(varChange, x)

	cprint(colors['orange'],'Equation: ',end='')
	print(eq)

	cprint(colors['orange'],'xData: ',end='')
	print(xData)

	c = []
	xSubs = []

	if points == 1:
		c = [2]
		xSubs = [0]
	elif points == 2:
		c = [1,1]
		xSubs = [-1/sqrt(3),1/sqrt(3)]
	elif points == 3:
		c = [5/9,8/9,5/9]
		xSubs = [-sqrt(3/5),0,sqrt(3/5)]
	elif points == 4:
		c = [
			(18-sqrt(30))/36,
			(18+sqrt(30))/36,
			(18+sqrt(30))/36,
			(18-sqrt(30))/36
		]
		xSubs = [
			-sqrt(525+70*sqrt(30/35)),
			-sqrt(525-70*sqrt(30/35)),
			sqrt(525-70*sqrt(30/35)),
			sqrt(525+70*sqrt(30/35))
		]
	elif points == 5:
		c = [
			(322-13*sqrt(70))/900,
			(322+13*sqrt(70))/900,
			128/225,
			(322+13*sqrt(70))/900,
			(322-13*sqrt(70))/900,
		]
		xSubs = [
			-sqrt(245+14*sqrt(70/21)),
			-sqrt(245-14*sqrt(70/21)),
			0,
			sqrt(245-14*sqrt(70/21)),
			sqrt(245+14*sqrt(70/21)),
		]
	elif points == 6:
		c = [
			0.171324492379170,
			0.360761573048139,
			0.467913934572691,
			0.467913934572691,
			0.360761573048131,
			0.171324492379170,
		]
		xSubs = [
			-0.932469514203152,
			-0.661209386466265,
			-0.238619186083197,
			0.238619186083197,
			0.661209386466265,
			0.932469514203152,
		]

	val=0
	for i in range(points):
		val += c[i]*eq.subs(x,xSubs[i])

	val = N(val)
	print('Estimated value:',val)
	tval = N(integrate(oldEq, (x, a, b)))
	print('True value:',tval)
	print('True error %:',abs(tval-val)/tval*100)

	return val



# --------------------------------------------------------Execution-------------------



def main():
	colorPallete = [
		colors["yellow"],
		colors["red"],
		colors["orange"],
		colors["orange"],
		colors["orange"]
	]
	header(colorPallete, "─", banner("NumSolvers",minus=-2), "Azcué", 1.0, "Numerical Numeolvers")

	# Insert code here, from examples
	
if __name__ == '__main__':
	main()
