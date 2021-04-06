from numpy import sin, cos, arcsin, arccos, arctan
from telescope_parameters import lat, skew, slope

#######################

#Description.

#Coordinate transformations between hour angle/declination and Utmost NS and MD angles.

#######################

#Define functions.

def HA(T, M):
        '''Returns hour angle (in rad) corresponding to Utmost NS and MD coordinates (called "T" and "M", both in rad), assuming the telescope is at latitude "lat", with skew "skew" and slope "slope" (all in rad).'''
	return arctan( (sin(M)*cos(skew)*cos(slope)+sin(T)*cos(M)*sin(skew)-cos(T)*cos(M)*cos(skew)*sin(slope))*1./(-sin(T)*cos(M)*sin(lat)*cos(skew)+cos(T)*cos(M)*(cos(lat)*cos(slope)-sin(lat)*sin(skew)*sin(slope))+sin(M)*(cos(lat)*sin(slope)+sin(lat)*sin(skew)*cos(slope))) )

def DEC(T, M):
        '''Returns declination (in rad) corresponding to Utmost NS and MD coordinates (called "T" and "M", both in rad), assuming the telescope is at latitude "lat", with skew "skew" and slope "slope" (all in rad).'''
	return arcsin( sin(T)*cos(M)*cos(lat)*cos(skew)+cos(T)*cos(M)*(sin(lat)*cos(slope)+cos(lat)*sin(skew)*sin(slope))+sin(M)*(sin(lat)*sin(slope)-cos(lat)*sin(skew)*cos(slope)) )

def NS(H, D):
        '''Returns Utmost NS angle (in rad) corresponding to sky position with hour angle "H" (in rad) and declination "D" (in rad), assuming the telescope is at latitude "lat", with skew "skew" and slope "slope" (all in rad).'''
	return arctan( (sin(H)*cos(D)*sin(skew)-cos(H)*cos(D)*sin(lat)*cos(skew)+sin(D)*cos(lat)*cos(skew))*1./(-sin(H)*cos(D)*cos(skew)*sin(slope)+cos(H)*cos(D)*(cos(lat)*cos(slope)-sin(lat)*sin(skew)*sin(slope))+sin(D)*(sin(lat)*cos(slope)+cos(lat)*sin(skew)*sin(slope))) )

def MD(H, D):
        '''Returns Utmost MD angle (in rad) corresponding to sky position with hour angle "H" (in rad) and declination "D" (in rad), assuming the telescope is at latitude "lat", with skew "skew" and slope "slope" (all in rad).'''
	return arcsin( sin(H)*cos(D)*cos(skew)*cos(slope)+cos(H)*cos(D)*(cos(lat)*sin(slope)+sin(lat)*sin(skew)*cos(slope))+sin(D)*(sin(lat)*sin(slope)-cos(lat)*sin(skew)*cos(slope)) )

