from __future__ import division
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from scipy.optimize import fsolve, least_squares
from scipy.io import savemat
import math
import pdb

import sys

ALPHA_ITS = 0.01
ALPHA_UPB = 56/180*np.pi
ALPHA_LWB = -np.pi/2
D_ITS = 0.01
A_ITS = 0.001
R = 0.035
F = 0.135
L = 0.145
class CSpace:
    def __init__(self,h=0.1, di = 0.3):
        #for h < f
        f = F
        b = F
        l = L
        r = R
        self.h = h
        self.d0 = di
        #self.alphai = np.arange(-np.pi/2, np.pi/2, 0.01)
        self.di = np.arange(D_ITS,di+D_ITS, D_ITS)
        self.ai = np.arange(r,h+r+A_ITS, A_ITS)

        dm, am = np.meshgrid(self.di, self.ai)
        dm = dm.reshape((-1))
        am = am.reshape((-1))
        da = np.stack([dm,am])
        da_lst = da.T.tolist()
        self.special = []
    def compute_alpha(self, d,a):
        '''
            note in this function, don't neet to compare with upper bound and lower bound because outside of this function, we will check the boundary, don't need to check a range, we will check outside
        '''
        #for h < f
        f = F
        l = L
        r = R
        h = self.h
        #condition
        alpha2_lb = -np.arcsin(h/l)
        #if (d-0.24)**2 < 0.0000001 and (a-r)**2<0.000001:
        #    pdb.set_trace()

        #print('x1',f+l+r)
        if d>=f+l+r:#R1
            return np.arange(0, ALPHA_UPB+ALPHA_ITS, ALPHA_ITS).reshape((-1,1)) if a == r else None 
        elif d >= np.sqrt(f**2-(h-r)**2)+l+r:#R2
            #print('x2',np.sqrt(f**2-(h-r)**2)+l+r)
            alpha_r2 = np.arccos((d-l-r)/f)
            return np.arange(alpha_r2, ALPHA_UPB+ALPHA_ITS, ALPHA_ITS).reshape((-1,1)) if a == r else None 
        else:
            #compute the x3
            def equations(p):
                dx, k, v = p
                return ((k+h-r)/f - (r-v)/r,\
                        (dx-l-r+v)/f-k/r,\
                        k**2+(r-v)**2-r**2)
            d_x3, k_r3, v_r3 = least_squares(equations, (0.03,0.01, 0.01), bounds=((0,0,0),(np.sqrt(f**2-(h-r)**2)+l+r,r,r))).x.tolist() 

            #print('x3',d_x3)
            if d >= d_x3:# R3
                def equations(p):
                    alpha, k, v = p
                    return (np.sin(alpha)-(k+h-r)/f,\
                            np.cos(alpha)-(d-l-r+v)/f,\
                            k**2+(r-v)**2-r**2)
                alpha_r3, k_r3, v_r3 = least_squares(equations, (0.1,0.01, 0.01), bounds=((0,0,0),(np.pi/2,r,r))).x.tolist()
                return np.arange(alpha_r3, ALPHA_UPB+ALPHA_ITS, ALPHA_ITS).reshape((-1,1)) if a == r else None  

            else:
                #compute the x4
                def equations(p):
                    theta, k, v = p
                    return (np.cos(theta) - k/r,\
                            np.sin(theta) - (k+h-r)/(f+l),\
                            k**2+(r-v)**2-r**2) #np.tan(theta) - (h - r + k)/(d-r+v))
                theta_x4, k1, v1 = least_squares(equations, (0.1,0.01, 0.01), bounds=((0,0,0),(np.pi/2, r,r))).x.tolist()
                d_x4 = (l+f)*np.cos(theta_x4)-v1+r 
                #print('x4',d_x4)
                if d >= d_x4:# R4
                    def equations(p):
                        theta, k, v, alpha, phi = p
                        return (np.pi/2-theta - (alpha+phi),\
                                l*np.sin(theta) + f*np.sin(theta+ alpha) - (k+h-r),\
                                r*np.sin(phi) - k,\
                                r*np.cos(phi) - (r - v),\
                                l*np.cos(theta)+f*np.cos(theta+alpha) - (d-r+v))
                    theta_r4, k1, v1, alpha_r4, phi1 = least_squares(equations, (0.1,0.01, 0.01, 0.1, 0.1),bounds=((0,0,0,0,0),(np.pi/2, r, r, np.pi/2, np.pi/2))).x.tolist()
                    return np.arange(alpha_r4, ALPHA_UPB+ALPHA_ITS, ALPHA_ITS).reshape((-1,1)) if a == r else None 
                else:
                    #compute the x5
                    d_x5 = np.sqrt(l**2-h**2)+f
                    #print('x5',d_x5)
                    if d >= d_x5:# R5
                        def equations(p):
                            theta, k, v, beta, phi = p
                            return (l*np.sin(theta) + f*np.sin(theta-beta)-(k+h-r),\
                                    2*np.pi - np.pi/2 - phi - theta - (np.pi-beta),\
                                    r*np.sin(phi) - k,\
                                    r*np.cos(phi) - (r-v),\
                                    l*np.cos(theta) + f*np.cos(theta - beta) - (d-r+v))
                        theta_r5, k1, v1, beta1, phi1 = least_squares(equations, (0.1,0.01, 0.01, 0.1,0.1),bounds=((0,0,0,0,0),(np.pi/2, r, r, np.pi/2, np.pi/2))).x.tolist()
                        alpha_r5 = -beta1
                        return np.arange(alpha_r5, ALPHA_UPB+ALPHA_ITS, ALPHA_ITS).reshape((-1,1)) if a == r else None
                    elif d >= r+l:# R6
                        #print('x6',r+l)
                        def equations(p):
                            alpha, k, v = p
                            return (np.tan(alpha)-(h-r+k)/(d-l-r+v),\
                                    np.cos(alpha)-(k/r),\
                                    np.sin(alpha)-(r-v)/r)
                        alpha_r6_up, k1, v1 = least_squares(equations, (0.1,0.01, 0.01), bounds=((0,0,0),(np.pi/2, r,r))).x.tolist() 
                        alpha_r6_lb = -np.arcsin(h/l)
                        if alpha_r6_up <= ALPHA_UPB:
                            alpha_r6_up = ALPHA_UPB
                        else:
                            alpha_r6_up = ALPHA_UPB

                        return np.arange(alpha_r6_lb, alpha_r6_up+ALPHA_ITS, ALPHA_ITS).reshape((-1,1)) if a == r else None

                    else:
                        # compute x7
                        d_x7 = np.sqrt(l**2-(h-r)**2)+r
                        #print('x7',d_x7)
                        if d >= d_x7:#R7
                            alpha_r7_ub = np.arcsin((d-r)/l)
                            alpha_r7_lb = alpha2_lb#-np.arcsin(h/l)
                            return np.arange(alpha_r7_lb, alpha_r7_ub+ALPHA_ITS, ALPHA_ITS).reshape((-1,1)) if a == r else None 
                        else:
                            if a <= h:
                                if d<=r:
                                    return None
                                #compute x8
                                def equations(p):
                                    theta, k = p
                                    return (np.sin(theta) - (h-a+k) / l,\
                                            np.cos(theta) - k/r)
                                theta_x8, k1 = least_squares(equations, (0.1,0.01),bounds=((0,0),(np.pi/2, r))).x.tolist()
                                d_x8 = l * np.cos(theta_x8) + r * np.sin(theta_x8) 
                                #if a == r:
                                #print('x8',d_x8)
                                self.special.append((d_x8, a, 0))
                                if d > d_x8:#R8, here is > not >=
                                    def equations(p):
                                        thetap, thetapp, v = p
                                        return (np.cos(thetapp) - (r-v)/r,\
                                                np.cos(thetap) - (d-r+v)/l,\
                                                np.sin(thetap) - (h-r+r*np.sin(thetapp))/l)
                                    theta1, theta2, v1 = least_squares(equations, (0.1,0.1,0.001), bounds = ((0, 0, 0),(np.pi/2, np.pi/2, r))).x.tolist()#fsolve(equations, (0.1, 0.1, 0.01))
                                    alpha_r8_ub = np.pi / 2 - theta1 - theta2 
                                    alpha_r8_lb = alpha2_lb#-np.arcsin(h/l)
                                    return np.arange(alpha_r8_lb,alpha_r8_ub+ALPHA_ITS, ALPHA_ITS).reshape((-1,1)) if a==r else None 
                                else:#from here we can use a>=r
                                    #compute x9
                                    theta_x9 = np.arcsin((h-a+r)/l)
                                    k1 = np.sin(np.pi/2-theta_x9)*r
                                    v1 = r - np.cos(np.pi/2-theta_x9)*r
                                    d_x9 = (h-a+k1)/np.tan(theta_x9) - v1 + r 
                                    #print('x9',d_x9)
                                    if d >= d_x9:#R9
                                            def equations(p):
                                                theta, k, v = p
                                                return (np.tan(theta) - (h-a+k)/(d-r+v),\
                                                        np.sin(np.pi/2 - theta) - k/r,\
                                                        np.cos(np.pi/2 - theta) - (r-v)/r)
                                            theta5, k1, v1 = least_squares(equations, (1.2, 0.01, 0.01),bounds=((0,0,0),(np.pi/2,r,r))).x.tolist()
                                            
                                            kp = np.sin(theta5) * l - h + a
                                            vp = l*np.cos(theta5) - (d-r)
                                            phi1 = np.arctan(kp/(r-vp))
                                            alpha_r9 = np.pi - 2*(np.pi/2 - theta5) - 2*(np.pi/2 - phi1)
                                            alpha_r9 = -alpha_r9
                                            return np.array(alpha_r9).reshape((-1,1)) 
                                    else:#R10
                                            def equations(p):
                                                theta, k, v = p
                                                return (np.tan(theta) - (h-a+k)/(d-r+v),
                                                     np.sin(np.pi/2-theta) - k/r,
                                                     np.cos(np.pi/2-theta) - (r-v)/r)
                                            theta6, k1, v1 = least_squares(equations, (1.2, 0.01,0.01), bounds=((0.,0.,0.),(np.pi/2, r,r))).x.tolist()  
                                            kpp = np.sin(theta6) * l - h + a
                                            vpp = np.cos(theta6) * l - d + r
                                            rho1 = np.arccos((kpp-r)/f)
                                            alpha_r10 = - (np.pi - rho1 - (np.pi/2 - theta6))
                                            return np.array(alpha_r10).reshape((-1,1)) 

                            else: # a > h:
                                if (d**2+(a-h)**2) < (r)**2:
                                    return None
                                #compute x8
                                mu1 = a-h
                                
                                def equations(p):
                                    mu, theta = p
                                    return (np.sin(theta) - mu/l,\
                                            np.cos(theta) - (mu1+mu)/r)                                
                                mu2, theta7 = least_squares(equations, (0.01,0.1),bounds=((0,0),(r,np.pi/2))).x.tolist()

                                if theta7 <0 or mu2<0:
                                    print('theta 7 solve invalid, skip with',d,a,theta7,mu2)
                                    return None
                                d_x8 = l*np.cos(theta7) + r * np.sin(theta7) 
                                if d > d_x8:
                                    return None
                                else:
                                    mu1 = a-h
                                    theta8 = np.arcsin((r-mu1)/l)
                                    d_x9 =  r*np.tan(theta8/2)+l*np.cos(theta8)
                                    if d >= d_x9:#R9
                                        def equations(p):
                                            theta, mu, u = p
                                            return (np.tan(theta) - mu/(d-r+u),\
                                                    np.cos(theta) - (mu1 + mu)/r,\
                                                    np.sin(theta) - (r-u)/r)
                                        theta9, mu2, u1 = least_squares(equations,(1.2, 0.01, 0.01),bounds=((0,0,0),(np.pi/2,r,r))).x.tolist()
                                        phi2 = np.arctan(r/(l-(mu2/np.sin(theta9))))
                                        alpha_r9 = -(np.pi - 2 * phi2)
                                    
                                        return np.array(alpha_r9).reshape((-1,1))
                                    else:#R10
                                        
                                        def equations(p):
                                            theta, mu, u = p
                                            return (np.tan(theta) - mu/(d-r+u),\
                                                    np.cos(theta) - (mu1 + mu)/r,\
                                                    np.sin(theta) - (r-u)/r)
                                        theta10, mu2, u1 = least_squares(equations,(1.2, 0.01, 0.01), bounds=((0,0,0),(np.pi/2, r, r))).x.tolist()
                                        if theta10 <0 or mu2<0 or u1<0:
                                            print('theta 10 solve invalid, skip with',d,a,theta10, mu2, u1)
                                            return None
                                        phi2 = np.arccos((l*np.sin(theta10)+a-h-r)/f)
                                        alpha_r10 = - (np.pi - phi2 - (np.pi/2 - theta10))
                                        return np.array(alpha_r10).reshape((-1,1))

    def compute_back_theta_t(self, d,a, alpha):
        '''
            note in this function, don't neet to compare with upper bound and lower bound because outside of this function, we will check the boundary, don't need to check a range, we will check outside
        '''
        #for h < f
        f = F
        b = F
        l = L
        r = R
        h = self.h
        #condition
        alpha2_lb = -np.arcsin(h/l)
        #if (d-0.17)**2 < 0.0000001:
        #    pdb.set_trace()

        if d>=f+l+r:#R1
            return alpha, 0, l
        elif d >= np.sqrt(f**2-(h-r)**2)+l+r:#R2
            return alpha, 0, l
        else:
            #compute the x3
            def equations(p):
                dx, k, v = p
                return ((k+h-r)/f - (r-v)/r,\
                        (dx-l-r+v)/f-k/r,\
                        k**2+(r-v)**2-r**2)
            d_x3, k_r3, v_r3 = least_squares(equations, (0.03,0.01, 0.01), bounds=((0,0,0),(np.sqrt(f**2+(h-r)**2)+l+r,r,r))).x.tolist()

            if d >= d_x3:# R3
                return 0,0,l
            else:
                #compute the x4
                def equations(p):
                    theta, k, v = p
                    return (np.cos(theta) - k/r,\
                            np.sin(theta) - (k+h-r)/(f+l),\
                            k**2+(r-v)**2-r**2)
                theta_x4, k1, v1 = least_squares(equations, (0.1,0.01, 0.01), bounds=((0,0,0),(np.pi/2, r,r))).x.tolist()
                d_x4 = (l+f)*np.cos(theta_x4)-v1+r 
                if d >= d_x4:# R4
                    def equations(p):#FP1
                        theta, k, v = p
                        return (np.sin(theta)-(r-v)/r,\
                                np.cos(theta)-k/r,\
                                np.tan(theta)-(h-r+k)/(d-l-r+v))
                    theta1p, k1, v1 = least_squares(equations, (0.1,0.01, 0.01), bounds=((0,0,0),(np.pi/2,r,r))).x.tolist() 
                    if alpha >= theta1p:
                        return 0, 0, l
                    else:
                        def equations(p):#FP2
                            theta, k, v = p
                            return(np.sin(theta+alpha)-(r-v)/r,\
                                   np.cos(theta+alpha)-k/r,\
                                   np.tan(theta+alpha) - (h-r+k-l*np.sin(theta))/(d-r+v-l*np.cos(theta)))
                        theta1p, k1, v1 = least_squares(equations, (0.01,0.01, 0.01),bounds=((0,0,0),(np.pi/2,r,r))).x.tolist()
                        return theta1p, theta1p, l 
                else:
                    #compute the x5
                    d_x5 = np.sqrt(l**2-h**2)+f
                    if d >= d_x5:# R5
                        def equations(p):#FP1
                            theta, k, v = p
                            return (np.sin(theta)-(r-v)/r,\
                                    np.cos(theta)-k/r,\
                                    np.tan(theta)-(h-r+k)/(d-l-r+v))
                        theta1p, k1, v1 = least_squares(equations, (0.1,0.01, 0.01), bounds=((0,0,0),(np.pi/2,r,r))).x.tolist() 
                        if alpha >= theta1p:
                            return 0, 0, l 
                        else:
                            def equations(p):#FP3 (actually the same as FP2)
                                theta, k, v = p
                                return(np.sin(theta+alpha)-(r-v)/r,\
                                       np.cos(theta+alpha)-k/r,\
                                       np.tan(theta+alpha) - (h-r+k-l*np.sin(theta))/(d-r+v-l*np.cos(theta)))
                            theta1p, k1, v1 = least_squares(equations, (0.1,0.01, 0.01),bounds=((0,0,0),(np.pi/2,r,r))).x.tolist()
                            return theta1p, theta1p, l 

                    elif d >= r+l:# R6
                        def equations(p):#FP1
                            theta, k, v = p
                            return (np.sin(theta)-(r-v)/r,\
                                    np.cos(theta)-k/r,\
                                    np.tan(theta)-(h-r+k)/(d-l-r+v))
                        theta1p, k1, v1 = least_squares(equations, (0.1,0.01, 0.01), bounds=((0,0,0),(np.pi/2,r,r))).x.tolist()
                        if alpha >= theta1p:
                            return 0, 0, l
                        else:
                            def equations(p):#FP2/FP3
                                theta, k, v = p
                                return(np.sin(theta+alpha)-(r-v)/r,\
                                       np.cos(theta+alpha)-k/r,\
                                       np.tan(theta+alpha) - (h-r+k-l*np.sin(theta))/(d-r+v-l*np.cos(theta)))
                            theta1p, k1, v1 = least_squares(equations, (0.01,0.01, 0.01),bounds=((0,0,0),(np.pi/2,r,r))).x.tolist()
                            return theta1p, theta1p, l  

                    else:
                        # compute x7
                        d_x7 = np.sqrt(l**2-(h-r)**2)+r
                        if d >= d_x7:#R7
                            def equations(p):#FP2/FP3
                                theta, k, v = p
                                return(np.sin(theta+alpha)-(r-v)/r,\
                                       np.cos(theta+alpha)-k/r,\
                                       np.tan(theta+alpha) - (h-r+k-l*np.sin(theta))/(d-r+v-l*np.cos(theta)))
                            theta1p, k1, v1 = least_squares(equations, (0.01,0.01, 0.01),bounds=((0,0,0),(np.pi/2,r,r))).x.tolist()
                            return theta1p, theta1p, l  
                        else:
                            if a <= h:
                                if d<=r:
                                    return None
                                #compute x8
                                def equations(p):
                                    theta, k = p
                                    return (np.sin(theta) - (h-a+k) / l,\
                                            np.cos(theta) - k/r)
                                theta_x8, k1 = least_squares(equations, (0.1,0.01),bounds=((0,0),(np.pi/2, r))).x.tolist()
                                d_x8 = l * np.cos(theta_x8) + r * np.sin(theta_x8) 
                                if d > d_x8:#R8, here is > not >=
                                    def equations(p):#FP2/FP3
                                        theta, k, v = p
                                        return(np.sin(theta+alpha)-(r-v)/r,\
                                               np.cos(theta+alpha)-k/r,\
                                               np.tan(theta+alpha) - (h-r+k-l*np.sin(theta))/(d-r+v-l*np.cos(theta)))
                                    theta1p, k1, v1 = least_squares(equations, (0.01,0.01, 0.01),bounds=((0,0,0),(np.pi/2,r,r))).x.tolist()
                                    if k1 < 1e-6 and v1 < 1e-6:#error1 that is on the upper bound
                                        def equations(p):
                                            thetap, thetapp, v = p
                                            return (np.cos(thetapp) - (r-v)/r,\
                                                    np.cos(thetap) - (d-r+v)/l,\
                                                    np.sin(thetap) - (h-r+r*np.sin(thetapp))/l)
                                        theta1, theta2, v1 = least_squares(equations, (0.1,0.1,0.001), bounds = ((0, 0, 0),(np.pi/2, np.pi/2, r))).x.tolist()#fsolve(equations, (0.1, 0.1, 0.01))
                                        return theta1, theta1, l
                                    else:
                                        return theta1p, theta1p, l   
                                else:#from here we can use a>=r
                                    #compute x9
                                    theta_x9 = np.arcsin((h-a+r)/l)
                                    k1 = np.sin(np.pi/2-theta_x9)*r
                                    v1 = r - np.cos(np.pi/2-theta_x9)*r
                                    d_x9 = (h-a+k1)/np.tan(theta_x9) - v1 + r 
                                    if d >= d_x9:#R9
                                        def equations(p):
                                            theta, k, v = p
                                            return (np.tan(theta) - (h-a+k)/(d-r+v),\
                                                    np.sin(np.pi/2 - theta) - k/r,\
                                                    np.cos(np.pi/2 - theta) - (r-v)/r)
                                        theta5, k1, v1 = least_squares(equations, (1.2, 0.01, 0.01),bounds=((0,0,0),(np.pi/2,r,r))).x.tolist()
                                        kp = np.sin(theta5) * l - h + a 
                                        return theta5+np.arccos((a-r)/b)+np.pi/2-np.pi, theta5, l-(kp-k1)/np.sin(theta5) 
                                    else:#R10
                                        def equations(p):
                                            theta, k, v = p
                                            return (np.tan(theta) - (h-a+k)/(d-r+v),
                                                 np.sin(np.pi/2-theta) - k/r,
                                                 np.cos(np.pi/2-theta) - (r-v)/r)
                                        theta6, k1, v1 = least_squares(equations, (1.2, 0.01,0.01), bounds=((0.,0.,0.),(np.pi/2, r,r))).x.tolist()  
                                        kpp = np.sin(theta6) * l - h + a 
                                        return theta6+np.arccos((a-r)/b)+np.pi/2-np.pi, theta6, l-(kpp-k1)/np.sin(theta6) 

                            else: # a > h:
                                if (d**2+(a-h)**2) <= (r)**2:
                                    return None
                                #compute x8
                                mu1 = a-h
                                
                                def equations(p):
                                    mu, theta = p
                                    return (np.sin(theta) - mu/l,\
                                            np.cos(theta) - (mu1+mu)/r)                                
                                mu2, theta7 = least_squares(equations, (0.01,0.1),bounds=((0,0),(r,np.pi/2))).x.tolist()

                                if theta7 <0 or mu2<0:
                                    print('theta 7 solve invalid, skip with',d,a,theta7,mu2)
                                    return None
                                d_x8 = l*np.cos(theta7) + r * np.sin(theta7) 
                                if d > d_x8:
                                    return None
                                else:
                                    mu1 = a-h
                                    theta8 = np.arcsin((r-mu1)/l)
                                    d_x9 =  r*np.tan(theta8/2)+l*np.cos(theta8)
                                    if d >= d_x9:#R9
                                        def equations(p):
                                            theta, mu, u = p
                                            return (np.tan(theta) - mu/(d-r+u),\
                                                    np.cos(theta) - (mu1 + mu)/r,\
                                                    np.sin(theta) - (r-u)/r)
                                        theta9, mu2, u1 = least_squares(equations,(1.2, 0.01, 0.01),bounds=((0,0,0),(np.pi/2,r,r))).x.tolist()
                                        return theta9+np.arccos((a-r)/b)+np.pi/2-np/pi, theta9, mu2/np.sin(theta9)  
                                    else:#R10
                                        
                                        def equations(p):
                                            theta, mu, u = p
                                            return (np.tan(theta) - mu/(d-r+u),\
                                                    np.cos(theta) - (mu1 + mu)/r,\
                                                    np.sin(theta) - (r-u)/r)
                                        theta10, mu2, u1 = least_squares(equations,(1.2, 0.01, 0.01), bounds=((0,0,0),(np.pi/2, r, r))).x.tolist()

                                        #check if back flipper touch curve
                                        if d < r:# in idea/back_touch_curve.png
                                            def equations(p):
                                                omega, k, v = p
                                                return(np.sin(omega)-k/r,\
                                                        np.cos(omega)-(r-v)/r,\
                                                        np.tan(omega)-(r-d-v)/(a-h-k))
                                            omega_10, k10, v10 = least_squares(equations, (0.3,0.01,0.01),bounds=((0,0,0),(np.pi/2, r,r))).x.tolist()
                                            if omega_10 > np.arccos((a-r)/b):
                                                print("!!!!!!!")
                                                return theta10+omega_10+np.pi/2-np.pi, theta10, mu2/np.sin(theta10)

                                        return theta10+np.arccos((a-r)/b)+np.pi/2-np.pi, theta10, mu2/np.sin(theta10)  
 




    def elem_func(self, x):
        if x[1]>self.h+R:
            return None
        alpham = self.compute_alpha(x[0],x[1])
        if alpham is None:
            #sys.stdout.write("None\n")
            return None
        else:
            alpha_nm = alpham.shape[0]
            if alpha_nm == 0:
                #sys.stdout.write("Empty in d=%f\n"%(x[0]))
                return None

            alphamc = [alf for alf in alpham[:,0].tolist() if alf <= ALPHA_UPB and alf >= ALPHA_LWB]
            alphamc = np.array(alphamc).reshape((-1,1))
            alpha_nm = alphamc.shape[0]
            if alpha_nm == 0:
                #sys.stdout.write("Empty in d=%f\n"%(x[0]))
                return None

            dam = np.tile(np.array(x).reshape((1,2)),(alpha_nm,1))
            return np.hstack([dam,alphamc]).reshape((-1,3)) 
                                     
                    
    @property
    def path(self):
        return self.get_path()
    def get_space(self):
        #for each d, we find its suitable a, alpha pair that is closest the current
        current_posi = np.array((self.d0,R,ALPHA_UPB)).reshape((1,3))
        self.space_lst = []
        add_lst = None
        for i in range(self.di.shape[0]):
            d = self.di[-(i+1)]
            print(d)
            dm = np.ones(self.ai.shape[0])*d
            da = np.stack([dm,self.ai])
            da_lst = da.T.tolist()        
            pts_lst = list(map(self.elem_func, da_lst))
            valid_pts_lst = []
            for i in range(len(da_lst)):
                if not pts_lst[i] is None:
                    if not np.isnan(pts_lst[i][0,2]):
                        valid_pts_lst.append(pts_lst[i])
                    else:
                        print('WARNING:find NAN')
            if (len(valid_pts_lst) == 0):
                print('No valid in d = %f'%d)
                continue
            pts = np.vstack(valid_pts_lst).reshape((-1,3))
            self.space_lst.append(pts)
        self.space = np.concatenate(self.space_lst, axis=0)
        return self.space
         

    def get_path(self):
        #for each d, we find its suitable a, alpha pair that is closest the current
        current_posi = np.array((self.d0,R,ALPHA_UPB)).reshape((1,3))
        self.path_lst = []
        self.path_lst.append(current_posi)
        for i in range(self.di.shape[0]):
            d = self.di[-(i+1)]
            dm = np.ones(self.ai.shape[0])*d
            da = np.stack([dm,self.ai])
            da_lst = da.T.tolist()        
            pts_lst = list(map(self.elem_func, da_lst))
            valid_pts_lst = []
            for i in range(len(da_lst)):
                if not pts_lst[i] is None:
                    if not np.isnan(pts_lst[i][0,2]):
                        valid_pts_lst.append(pts_lst[i])
                    else:
                        print('WARNING:find NAN')
            if (len(valid_pts_lst) == 0):
                print('No valid in d = %f'%d)
                continue
            pts = np.vstack(valid_pts_lst)
            #dist = np.sum((pts-current_posi)**2,axis=1)
            dist = (pts-current_posi)**2
            dist[:,1]*=100
            dist = np.sum(dist,axis=1)
            min_id = np.argmin(dist)
            current_posi = pts[min_id,:].reshape((1,3))
            self.path_lst.append(current_posi)
        return self.path_lst
         
    def tuple3_to_tuple6(self, p):
        '''
           p: np.array[1,3]
           return: np.array[1,6]
        '''
        print(p)
        d,a,alpha = p[0].tolist()
        btt = self.compute_back_theta_t(d,a,alpha)
        return np.array([d,a,alpha,btt[0],btt[1],btt[2]]).reshape((-1,1))

if __name__ == '__main__':
    draw_path = True
    if draw_path:
        cp = CSpace(h=0.1,di=0.34)
        path_lst = cp.path
        pts = np.concatenate(path_lst)


        fig = plt.figure()
        ax = Axes3D(fig)

        ax.set_xlabel('d')
        ax.set_ylabel('a')
        ax.set_zlabel('alpha')

        ax.scatter(pts[:,0].tolist(), pts[:,1].tolist(), pts[:,2].tolist(),'b.')
        #ax.scatter(cp.special[0,0],cp.special[0,1],cp.special[0,2],'r.')

        #plt.savefig(__file__+'.png')
        plt.savefig('path.png')
        plt.show()
        print("finish cspace")

        savemat(__file__+'.mat', {'pc':pts,'D_ITS':D_ITS,'A_ITS':A_ITS,'ALPHA_ITS':ALPHA_ITS, 'st':np.array((0.3, R, np.pi/3)).reshape((1,3)), 'ed':np.array((0.,0.05+R, np.pi/3))})

        alpha_l = []
        b_alpha_l = []
        for p in path_lst: 
            print(p)

        f = open('../../../../PairData.txt','w')
        f.write('[')

        print('generating extention path ################')
        f_first=True
        for p in path_lst:
            if f_first:
                f_first=False
            else:
                f.write(',')
            p = p[0]
            #print(p)
            d,a,alpha = p[0], p[1] ,p[2] 
            back_alpha, theta, t = cp.compute_back_theta_t(d,a,alpha)
            alpha_l.append(alpha)
            b_alpha_l.append(back_alpha)
            print(d,a,alpha,back_alpha, theta, t)
            f.write('(%f,%f,%f,%f,%f,%f)'%(d,a,alpha,back_alpha, theta, t))
        f.write(']')
        f.close()
        import scipy.io
        scipy.io.savemat('tmp.mat',{'alpha':np.stack([np.array(alpha_l), np.array(b_alpha_l)])})

            #print(d,a,alpha, btt)
    else: #draw space
        cp = CSpace(h=0.095,di=0.4)
        pts = cp.get_space()


        fig = plt.figure()
        ax = Axes3D(fig)

        ax.set_xlabel('d')
        ax.set_ylabel('a')
        ax.set_zlabel('alpha')

        ax.scatter(pts[:,0].tolist(), pts[:,1].tolist(), pts[:,2].tolist(),c=pts[:,2].tolist(),s=2,cmap='cool')
        #ax.scatter([0.13656500283928868], [0.035], [1.],c='b',s=2) 
        #ax.plot(pts[:,0].tolist(), pts[:,1].tolist(),'b.')

        cp.special = np.array(cp.special)
        pdb.set_trace()
        ax.scatter(cp.special[:,0].tolist(),cp.special[:,1].tolist(),cp.special[:,2].tolist(),c='r',s=2)

        plt.savefig('space.png')
        plt.show()
        print("finish cspace")
 
 
