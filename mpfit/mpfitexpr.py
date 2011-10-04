""" 
	Copyright (C) 2009 Sergey Koposov

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import mpfit 
import re
import numpy

def mpfitexpr(func, x, y, err , start_params, check=True, full_output=False, **kw):
	"""Fit the used defined expression to the data
	Input:
	- func: string with the function definition 
	- x: x vector
	- y: y vector
	- err: vector with the errors of y
	- start_params: the starting parameters for the fit
	Output:
	- The tuple (params, yfit) with best-fit params and the values of func evaluated at x
	Keywords:
	- check: boolean parameter. If true(default) the function will be checked for sanity
	- full_output: boolean parameter. If True(default is False) then instead of best-fit parameters the mpfit object is returned
	Example:
	params,yfit=mpfitexpr('p[0]+p[2]*(x-p[1])',x,y,err,[0,10,1])
	
	If you need to use numpy functions in your function, then
		you must to use the full names of these functions, e.g.:
		numpy.sin, numpy.cos etc.
	
	This function is motivated by mpfitexpr() from wonderful MPFIT IDL package
		written by Craig Markwardt	
	
	"""

	def myfunc(p,fjac=None,x=None, y=None, err=None):
		return [0, eval('(y-(%s))/err'%func)]

	myre = "[^a-zA-Z]p\[(\d+)\]"
	r = re.compile(myre)
	maxp = -1
	for m in re.finditer(r,func):
		curp = int(m.group(1))
		maxp = curp if curp > maxp else maxp	
	if check:
		if maxp == -1: 
			raise Exception("wrong function format")
		if maxp + 1 != len(start_params):
			raise Exception("the length of the start_params != the length of the parameter verctor of the function")
	fa={'x' : x, 'y' : y,'err' : err}
	res = mpfit.mpfit(myfunc,start_params,functkw=fa,**kw)
	yfit = eval(func, globals(), {'x':x, 'p': res.params})
	if full_output:
		return (res, yfit)
	else:
		return (res.params, yfit)
