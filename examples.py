
	#Some of this examples may be calling some outdated functions
	"""
	r=12
	L=5
	V = ((r**2) * acos((r-x)/r)-(r-x)*sqrt(2*r*x-x**2))*L-8
	
	FalsePosition(V, [0, 1, 4.0])
	#plot(V,0, (x, 0,4),rendering_kw=[{}, {"linestyle": "--"}])


	eq1 = ((1000+(x1*(x2**2))-((x1-1)**3))**(1/3))/x1
	eq2 = ((((8**3)+((x1**2)/x2))**(1/3))-12)/x2
	FixedPointSE([eq1,eq2],[3,2])


	hours = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
	population = [2.32, 2.70, 2.72, 3, 3.30, 3.65, 4, 4.5, 5.2, 5.44, 5.90, 6.64, 7.40, 8.3, 8.96]
	regr = leastSquares(hours, population, model='exp')
	regr.subs(20)



	time = [0,3,5,7]
	temp = [150,120,100,85]

	inter = Newton(time, temp)


	target = 95
	eq = inter.eq - target

	Bisection(eq, [5.0, 5.5, 7.0], maxError=0.00004, solution=False)
	plot(eq,0, (x, 5,6),rendering_kw=[{}, {"linestyle": "--"}])

	u = 1850
	m0=160000
	q=2500
	g=9.81
	v = u*ln(m0/(m0-q*x))-g*x
	n = 6
	xData = list(np.linspace(0,30,n+1))
	simpsons38Rule(v, xData)

	R=5
	L=0.5
	di = (20*exp(-3*x)-R*y)/L

	y0=0
	xData = list(np.arange(0.0,0.4,0.02))
	RungeKutta(xData=xData, y0=y0, eq=di, method="rk4")
	"""
	""" integration HW
	eq = 1-exp(-x)

	xData = [0,4]

	gaussQuadrature(eq,xData,1)
	rectangleRule(eq, xData)
	trapezoidalRule(eq, xData)

	xData = [0,2,4]
	trapezoidalRule(eq, xData)
	xData = [0,1,2,3,4]
	trapezoidalRule(eq, xData)

	xData = [0,2,4]
	simpsons13Rule(eq, xData)
	xData = [0,1,2,3,4]
	simpsons13Rule(eq, xData)
	
	n = 3
	xData = list(np.linspace(0,4,n+1))
	
	simpsons38Rule(eq, xData)

	n = 5
	xData = list(np.linspace(0,4,n+1))

	v1 = simpsons13Rule(eq, xData[:3]) # first 2 segments
	v2 = simpsons38Rule(eq, xData[-4:]) # last 3 segments
	val = v1+v2

	print()
	cprint(colors["orange"], "Numerical Integration | ", end='')
	cprint(colors["pink"],'Mixed Simpson rule')

	print('Estimated value:',val)
	tval = N(integrate(eq, (x, xData[0], xData[-1])))
	print('True value:',tval)
	print('True error %:',abs(tval-val)/tval*100)
	"""



	#b not included
	#list(np.arange(a,b,h))

	#n spaces
	#list(np.linspace(a,b,n+1))

	"""ecuaciones diferenciales 
	k=(1/3)*ln(13/23)
	Tm=70
	dT = k*(y-70)

	y0=300.0
	xData = list(np.arange(0,26,1))

	RungeKutta(xData=xData, y0=y0, eq=dT, method="rk4")


	
	dT = -y+x**2

	y0=1.0
	xData = list(np.arange(0,5,0.5))
	RungeKutta(xData=xData, y0=y0, eq=dT, method="midpoint")


	eq = y*(2*x-1.5)

	y0=2.0
	xData = list(np.arange(0,3,0.5))
	RungeKutta(xData=xData, y0=y0, eq=eq, method="euler")


	a=0.0002
	eq = y*(1-a*y)
	y0=1.0
	xData = list(np.arange(0,21,1))
	RungeKutta(xData=xData, y0=y0, eq=eq, method="ralston")
	"""
	
	"""
	eq = x + y  # dy/dx = x + y
	xData = [0.0, 0.1, 0.2, 0.3, 0.4]
	y0 = 1.0
	RungeKutta(xData=xData, y0=y0, eq=eq, method="euler")
	RungeKutta(xData=xData, y0=y0, eq=eq, method="heun")
	RungeKutta(xData=xData, y0=y0, eq=eq, method="ralston")
	RungeKutta(xData=xData, y0=y0, eq=eq, method="midpoint")
	RungeKutta(xData=xData, y0=y0, eq=eq, method="rk4")


	eq = y*x**2-1.1*y
	xData = list(np.arange(0,2.5,0.5))
	RungeKutta(xData=xData, y0=1.0, eq=eq, method="midpoint")

	RungeKutta(xData=xData, y0=1.0, eq=eq, method="rk4")
	"""




	"""
	eq = sin(x)
	h = (pi/12)
	xValue = (pi/4)
	simpleNumericalDerivative(eq,h,xValue)
	#-----------------------------------------------------------------
	print()
	xData = [0,25,50,75,100,125]
	yData = [0,32,58,78,92,100]

	numericalDerivative(xData=xData, yData=yData, quality=1)
	print()
	numericalDerivative(xData=xData, yData=yData, degree=2, quality=1)
	print()
	#-----------------------------------------------------------------

	xData = [0, 0.52, 1.04, 1.75, 2.37, 3.25, 3.83]
	yData = [153, 185, 208, 249, 261, 271, 273]
	diff = unequalNumericalDerivative(xData=xData, yData=yData).diffData
	print()
	unequalNumericalDerivative(xData=xData, yData=diff)
	print()
	"""
	"""

	xData = [0,19,32,47,50,58]
	yData = [0,50,72,85,102,130]
	diff = unequalNumericalDerivative(xData=xData, yData=yData).diffData

	unequalNumericalDerivative(xData=xData, yData=diff)
	eq = (cos(x)**2)/(x**2)

	test = 1+3/2*z
	print(N(integrate(eq, (x, 1, 4))))

	plot(eq,0, (x, 1,4),rendering_kw=[{}, {"linestyle": "--"}])
	plot(eq.subs(x,test),0, (z, 1,4),rendering_kw=[{}, {"linestyle": "--"}])
	print(N(integrate(eq.subs(x,test), (z, 1, 4))))
	"""

	
	


	"""
	data = [0.90, 1.42, 1.30, 1.55, 1.63,
			1.32, 1.35, 1.47, 1.95, 1.66,
			1.96, 1.47, 1.92, 1.35, 1.05,
			1.85, 1.74, 1.65, 1.78, 1.71,
			2.38, 1.82, 2.06, 2.14, 1.27]
	statisticsInfo(data)
	

	

	yData=[0, 2, 4, 6, 9, 11, 12, 15, 17, 19]
	xData=[5, 6, 7, 6, 9, 8, 8, 10, 12, 12]
	leastSquares(xData,yData)

	xData=[2.3026,2.9957,3.4012,3.6889,3.9120,4.0943,4.2485,4.3820]
	yData=[3.2189,4.2485,5.9402,6.3099,6.4135,7.1066,6.7214,7.2793]
	leastSquares(xData, yData)

	
	numInt([3,3.2,3.4,3.6,3.8,4], cos(x)/(x**2))
	fx=0.2+25*x-200*x**2+675*x**3-900*x**4+400*x**5
	numInt([0,0.8], fx,degree=1)
	print(integrate(fx,x))
	"""

	#determinación r^2
	#
	"""
	# Data points (time, temperature)
	time_data = [0, 2, 4, 6]
	temp_data = [80, 75.6, 71.2, 67.8]

	inter = Newton(time_data, temp_data)


	target_temp = 73
	f = inter.eq - target_temp

	Bisection(f, [2.0, 2.1, 4.0], maxError=0.00004, solution=False)

	#-----------------------------------------------------------------------------------------
	time_data = [0.4, 1.0, 1.6, 2.4]
	conc_data = [4.5, 3.7, 2.3, 1.9]

	inter = Newton(time_data, conc_data)

	target_conc = 3.2
	f = inter.eq - target_conc

	NewtonRaphson(f, [1.3], maxError=0.00004, solution=False)


	#-----------------------------------------------------------------------------------------
	f_data = [3.6,1.8,1.2,0.9]
	x_data = [1,2,3,4]

	inter = Lagrange(f_data, x_data)

	target_f = 1.7
	f_equation = inter.eq - target_f


	Bisection(f_equation, [1.9, 1.99, 2.3], maxError=0.00004, solution=False)
	"""
	"""Examen intersemestral
	#-----------------------------------------------Ejercicio 1 | Rocket
	u = 1800
	m0 = 160000
	q=2600
	v=750
	g=9.81
	fx = u*ln(m0/(m0-q*x))-g*x-v
	#plot(fx,0, (x, 10,50),rendering_kw=[{}, {"linestyle": "--"}])
	Bisection(fx, [10.0, 1, 50.0])
	FalsePosition(fx, [10, 1, 50])

	#-----------------------------------------------Ejercicio 2 | Throw
	v0= 30
	xf = 90
	yf=1
	y0=1.8
	fx = tan(x*2*pi/360)*xf - (g/(2*(v0**2)*(cos(x*2*pi/360)**2)))*(xf**2)+y0-yf
	NewtonRaphson(fx, [35])
	NewtonRaphson(fx, [60])
	Secant(fx, [36], lastRow=[30])
	Secant(fx, [60], lastRow=[44])
	#plot(fx,0, (x, 1,45),rendering_kw=[{}, {"linestyle": "--"}])
	#plot(fx,0, (x, 51,51.6),rendering_kw=[{}, {"linestyle": "--"}])

	
	#-----------------------------------------------Ejercicio 3 | Simple bisection
	fx = 4-log(3-x,10)-3
	Bisection(fx, [-10.0, 1, 1.0])
	#plot(fx,0, (x, -10,1),rendering_kw=[{}, {"linestyle": "--"}])


	#-----------------------------------------------Ejercicio 4 | Lagrange inter
	inter = LagrangeInter(
		[8.1,8.3,8.6,8.7],
		[16.94410,17.56492,18.50515,18.82091]
	)
	inter.subs(8.4)
	#plot(inter.eq,0, (x, 8,9),rendering_kw=[{}, {"linestyle": "--"}])

	#-----------------------------------------------Ejercicio 5 | Table
	FixedPointLE([
			(cos(x1)+x2*sin(x1))/2, 
			(4-x1)**(1/2)
			],
			[1,1]
		)

	#-----------------------------------------------Ejercicio 5 | Inter of Interpolations
	inter = Newton(
		[5,7],
		[59.8763748,55.3946902]
	)
	k10t6 = inter.subs(6)

	inter = Newton(
			[5,7],
			[62.5107146,60.0027463]
		)
	k25t6 = inter.subs(6)

	inter = Newton(
			[10,25],
			[k10t6,k25t6]
		)
	k15t6 = inter.subs(15)
	"""
	"""Differentiation
	eq = sin(x)
	h = (pi/12)
	xValue = (pi/4)
	simpleNumericalDerivative(eq,h,xValue)
	
	xData = [0,25,50,75,100,125]
	yData = [0,32,58,78,92,100]
	numericalDerivative(xData=xData, yData=yData)

	numericalDerivative(xData=xData, yData=yData, degree=2)


	xData = [0, 0.52, 1.04, 1.75, 2.37, 3.25, 3.83]
	yData = [153, 185, 208, 249, 261, 271, 273]
	diff = unequalNumericalDerivative(xData=xData, yData=yData).diffData

	unequalNumericalDerivative(xData=xData, yData=diff)

	diff = exp((x**3)/3-1.1*x)

	eq = y*x**2-1.1*y
	#y(o)=1
	initialValue = 1
	h=0.25
	start = 0
	end = 2

	datos = []
	i = start
	while i != end:
		
		datos.append((i,eq.subs(x,i)))
		i += h

	print(datos)

	plot(diff,0, (x, 0,2),rendering_kw=[{}, {"linestyle": "--"}])
	"""
	"""HW
	points = [0,1,1.6,2]
	a = points[0]
	b = points[-1]
	fx= exp(-x**2)
	temp2 = 0
	for i, A in enumerate(points):
		if i != len(points)-1:
			
			B = points[i+1]
			temp = (B-A)*((fx.subs(x,A)+fx.subs(x,B))/2)
			print(f'{A}->{B}:',N(temp))
			temp2 += temp
	
	val = N(temp2)
	print('Estimated value:',val)
	tval = N(integrate(fx, (x, a, b)))
	print('True value:',tval)
	print('Error %:',abs(tval-val)/tval*100)
	
	print('----------------')
	#compound trapezoidal rule
	#h must be equal at each point
	a = 0
	b = 2
	n = 5 # num total de puntos menos 1
	h=(b-a)/n
	print(h)
	val = N((h/2)*(fx.subs(x,a)+2*sum([fx.subs(x,i) for i in [i*h+a for i in range(1,n)]])+fx.subs(x,b)))
	print('Estimated value:',val)
	tval = N(integrate(fx, (x, a, b)))
	print('True value:',tval)
	print('Error %:',abs(tval-val)/tval*100)


	print('----------------')
	# Simpson 1/3
	fx=0.2+25*x-200*x**2+675*x**3-900*x**4+400*x**5
	a = 0
	b = 0.8
	n = 4 # num total de puntos menos 1
	h=(b-a)/n
	print(h)
	sumaPar = sum([fx.subs(x,i) for i in [i*h+a for i in range(1,n) if i%2==0 ]])
	sumaImpar = sum([fx.subs(x,i) for i in [i*h+a for i in range(1,n) if i%2==1 ]])
	
	val = N((h/3)*(fx.subs(x,a) + 4*sumaImpar + 2*sumaPar + fx.subs(x,b)))
	print('Estimated value:',val)
	tval = N(integrate(fx, (x, a, b)))
	print('True value:',tval)
	print('Error %:',abs(tval-val)/tval*100)



	print('-----------')
	"""

	""" 3
	age = [31, 32, 33, 34, 35, 36, 37, 38, 39, 40]
	income = [23500, 24000, 25000, 26700, 27500, 29200, 33000, 35100, 37400, 39500]
	regr = leastSquares(age, income)
	regr.subs(48)

	hours = [3, 2, 0, 0.5, 2.5, 1.5, 1.8, 1, 0.3, 2.2]
	grade = [10, 9.5, 5, 6.5, 10, 9, 8, 8.5, 6, 8.5]
	regr = leastSquares(hours, grade)
	regr.subs(1.9)

	age = [31, 32, 33, 34, 35, 36, 37, 38, 39, 40]
	income = [23500, 24000, 25000, 26700, 27500, 29200, 33000, 35100, 37400, 39500]
	regr = leastSquares(age, income, model='pow')
	regr.subs(52)

	xData=[0, 2, 4, 6, 9, 11, 12, 15, 17, 19]
	yData=[5, 6, 7, 6, 9, 8, 8, 10, 12, 12]
	leastSquares(xData,yData)

	cigs = [5,23,25,48,17,8,4,26,11,19,14,35,29,4,23]
	life = [80,78,60,53,85,84,73,79,81,75,68,72,58,92,65]
	regr=leastSquares(cigs, life, model='pow')
	regr.subs(47)

	time = [1,2,4,6]
	speed = [5,12,14,15]
	inter = Newton(time,speed)
	yTarget = 13
	f = inter.eq - yTarget
	Bisection(f, [2, 3.9, 4])
	#plot(f,0, (x, 0,15),rendering_kw=[{}, {"linestyle": "--"}])
	"""
	"""HM7
	data = [0.90, 1.42, 1.30, 1.55, 1.63,
			1.32, 1.35, 1.47, 1.95, 1.66,
			1.96, 1.47, 1.92, 1.35, 1.05,
			1.85, 1.74, 1.65, 1.78, 1.71,
			2.38, 1.82, 2.06, 2.14, 1.27]
	statisticsInfo(data)
	

	

	yData=[0, 2, 4, 6, 9, 11, 12, 15, 17, 19]
	xData=[5, 6, 7, 6, 9, 8, 8, 10, 12, 12]
	leastSquares(xData,yData)

	xData=[2.3026,2.9957,3.4012,3.6889,3.9120,4.0943,4.2485,4.3820]
	yData=[3.2189,4.2485,5.9402,6.3099,6.4135,7.1066,6.7214,7.2793]
	leastSquares(xData, yData)

	
	numInt([3,3.2,3.4,3.6,3.8,4], cos(x)/(x**2))
	fx=0.2+25*x-200*x**2+675*x**3-900*x**4+400*x**5
	numInt([0,0.8], fx,degree=1)
	print(integrate(fx,x))
	"""

	#determinación r^2
	#
	"""
	# Data points (time, temperature)
	time_data = [0, 2, 4, 6]
	temp_data = [80, 75.6, 71.2, 67.8]

	inter = Newton(time_data, temp_data)


	target_temp = 73
	f = inter.eq - target_temp

	Bisection(f, [2.0, 2.1, 4.0], maxError=0.00004, solution=False)

	#-----------------------------------------------------------------------------------------
	time_data = [0.4, 1.0, 1.6, 2.4]
	conc_data = [4.5, 3.7, 2.3, 1.9]

	inter = Newton(time_data, conc_data)

	target_conc = 3.2
	f = inter.eq - target_conc

	NewtonRaphson(f, [1.3], maxError=0.00004, solution=False)


	#-----------------------------------------------------------------------------------------
	f_data = [3.6,1.8,1.2,0.9]
	x_data = [1,2,3,4]

	inter = Lagrange(f_data, x_data)

	target_f = 1.7
	f_equation = inter.eq - target_f


	Bisection(f_equation, [1.9, 1.99, 2.3], maxError=0.00004, solution=False)
	"""
	"""Examen intersemestral
	#-----------------------------------------------Ejercicio 1 | Rocket
	u = 1800
	m0 = 160000
	q=2600
	v=750
	g=9.81
	fx = u*ln(m0/(m0-q*x))-g*x-v
	#plot(fx,0, (x, 10,50),rendering_kw=[{}, {"linestyle": "--"}])
	Bisection(fx, [10.0, 1, 50.0])
	FalsePosition(fx, [10, 1, 50])

	#-----------------------------------------------Ejercicio 2 | Throw
	v0= 30
	xf = 90
	yf=1
	y0=1.8
	fx = tan(x*2*pi/360)*xf - (g/(2*(v0**2)*(cos(x*2*pi/360)**2)))*(xf**2)+y0-yf
	NewtonRaphson(fx, [35])
	NewtonRaphson(fx, [60])
	Secant(fx, [36], lastRow=[30])
	Secant(fx, [60], lastRow=[44])
	#plot(fx,0, (x, 1,45),rendering_kw=[{}, {"linestyle": "--"}])
	#plot(fx,0, (x, 51,51.6),rendering_kw=[{}, {"linestyle": "--"}])

	
	#-----------------------------------------------Ejercicio 3 | Simple bisection
	fx = 4-log(3-x,10)-3
	Bisection(fx, [-10.0, 1, 1.0])
	#plot(fx,0, (x, -10,1),rendering_kw=[{}, {"linestyle": "--"}])


	#-----------------------------------------------Ejercicio 4 | Lagrange inter
	inter = LagrangeInter(
		[8.1,8.3,8.6,8.7],
		[16.94410,17.56492,18.50515,18.82091]
	)
	inter.subs(8.4)
	#plot(inter.eq,0, (x, 8,9),rendering_kw=[{}, {"linestyle": "--"}])

	#-----------------------------------------------Ejercicio 5 | Table
	FixedPointLE([
			(cos(x1)+x2*sin(x1))/2, 
			(4-x1)**(1/2)
			],
			[1,1]
		)

	#-----------------------------------------------Ejercicio 5 | Inter of Interpolations
	inter = Newton(
		[5,7],
		[59.8763748,55.3946902]
	)
	k10t6 = inter.subs(6)

	inter = Newton(
			[5,7],
			[62.5107146,60.0027463]
		)
	k25t6 = inter.subs(6)

	inter = Newton(
			[10,25],
			[k10t6,k25t6]
		)
	k15t6 = inter.subs(15)
	"""

