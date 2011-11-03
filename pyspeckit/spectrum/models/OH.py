"""
OH line fitter
"""
import redshiftedgroup


freq_dict={
'OH12':1.61223e9,
'OH11':1.66540e9,
'OH22':1.66736e9,
'OH21':1.72053e9,
}


OH = redshiftedgroup.redshiftedgroup(freq_dict)
OHfitter = OH.fitter
OHvheightfitter = OH.vheight_fitter
