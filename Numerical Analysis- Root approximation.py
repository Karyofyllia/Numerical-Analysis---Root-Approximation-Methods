#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math
import matplotlib.pyplot as plt
import numpy as np
import random
from scipy.linalg import solve
import decimal


# In[2]:


def f(x):
    result =math.e**((math.sin(x))**3) + (x**6) - (2*(x)**4) - (x**3) - 1
    return result
print(f(-2),' f(-2)')
print(f(2),' f(2)')
print(f(1),' f(1)')


# In[3]:


def f2(x):
    result= 94*(math.cos(x)**3)-24*math.cos(x)+177*(math.sin(x)**2)-108*(math.sin(x)**4)-72*(math.cos(x)**3)*(math.sin(x)**2)-65
    return result


# In[4]:


e = 1/2 * 10**-5

def bisection_method(a,b,iterat):
    
   
    if  f(a)*f(b)<0:
        m = (a + b ) / 2
        if (np.abs(f(m)) < e):
            m= m + e
            return m, iterat

        if f(a)*f(m)<0:
            iterat=iterat+1
            return bisection_method(a,m,iterat)
        else:
            iterat=iterat+1
            return bisection_method(m,b,iterat)


# In[5]:


bisection_method(-2,1,0)


# In[6]:


bisection_method(1,2,0)


# In[7]:


def derf(x):
    result = 3*math.e**((math.sin(x))**3)*math.cos(x)*(math.sin(x))**2 + 6*x**5 -8*x**3 - 3*x**2
    return result


# In[8]:



def newton_raphson_method(x):
    iterat=0
    while (np.abs(f(x)) > e):
        iterat=iterat+1
        x = x - f(x)/derf(x)
        #print(x, iterat)
    return x, iterat


# In[9]:


newton_raphson_method(2)


# In[10]:


newton_raphson_method(-2)


# In[11]:


newton_raphson_method(0)


# In[12]:


def secant_method(a,b):
    iterat=0
    while (np.abs(f(a)- f(b)) > e):
        iterat=iterat+1
        temp=a
        a = a - (a-b) * f(a)/(f(a)-f(b))
        b= temp
    return a,iterat


# In[13]:


secant_method(-2,2)


# In[14]:


secant_method(-2,0) 


# In[15]:


def dderf(x):
    result= 9*math.e**((math.sin(x))**3)* ((math.cos(x))**2)*(math.sin(x)**4)-3*math.e**((math.sin(x))**3)*(math.sin(x)**3)+6*math.e**((math.sin(x))**3)*((math.cos(x))**2)*(math.sin(x))+30*x**4-24*x**2-6*x
    return result


# In[16]:


def derf2(x):
    result= 216*(math.cos(x)**2)*(math.sin(x)**3)-432*math.cos(x)*(math.sin(x)**3)-144*(math.cos(x)**4)*math.sin(x)-282*(math.cos(x)**2)*math.sin(x)+354*math.cos(x)*math.sin(x)+24*math.sin(x)
    return result


# In[17]:


def dderf2(x):
    result = (432-432*math.cos(x))*(math.sin(x)**4)+(1224*(math.cos(x)**2)-1296*(math.cos(x)**2)+564*math.cos(x)-354)*math.sin(x)**2-144*math.cos(x)**5-282*math.cos(x)**3+354*math.cos(x)**2+24*math.cos(x)
    return result


# In[18]:


def modified_NR(f2,derf2,dderf2,x):
    iterat=0
    while (np.abs(f2(x)) > e):
        iterat=iterat+1
        x= x - 1/((derf2(x)/f2(x)) - 1/2*(dderf2(x)/derf2(x))) #modified
    return x,iterat


# In[19]:


modified_NR(f2,derf2,dderf2,1)


# In[20]:


modified_NR(f2,derf2,dderf2,3)


# In[21]:


modified_NR(f2,derf2,dderf2,2)


# In[22]:


def modified_BM(a,b,iterat):    
   
    if  f2(a)*f2(b)<0:
        m = random.uniform(a,b) #modified
        if (np.abs(f2(m)) < e):
            m= m + e
            return m, iterat

        if f2(a)*f2(m)<0:
            iterat=iterat+1
            return modified_BM(a,m,iterat)
        else:
            iterat=iterat+1
            return modified_BM(m,b,iterat)


# In[23]:


modified_BM(0.9,2,0)


# In[24]:


def modified_SM(x1,x2,x3): 
    iterat=0
    q=f2(x1)/f2(x2)
    r=f2(x3)/f2(x2)
    s=f2(x3)/f2(x1)
    while (np.abs(f2(x3)) > e):
        iterat=iterat+1
        x4 = x3 - (r*(r-q)*(x3-x2)+(1-r)*s*(x3-x1))/((q-1)*(r-1)*(s-1))
        x1=x2
        x2=x3
        x3=x4
    return x3,iterat


# In[25]:


modified_SM(0,1,3)


# In[26]:


def find_LPU(A):
    
    L = np.zeros((A.shape[0],A.shape[1]))  #initialize with 0
    P= np.diag(np.ones(A.shape[0]))        #initialize: 1 diagonial 0 elsewhere.
    
    for j in range(A.shape[0]-1): 
            if (A[j][1]==0 or A[j][1]<A[j+1][1]):
                Pnew= np.diag(np.ones(A.shape[0]))
                A[[j+1,j],:]=A[[j,j+1],:] 
                
                Pnew[[j+1,j],:]=Pnew[[j,j+1],:] #for every line shift of A , do the same to P.
                P = np.dot(Pnew,P)              #multiple older with present P

    for i in range(0, A.shape[0]): 
       
            
        L[i][i]=1

        for j in range(i + 1, A.shape[0]): #Gauss
            if (A[i][j]): 
                l=A[j][i]/A[i][i]
                A[j] = A[j] - A[i] * (A[j][i]/A[i][i])
                L[j][i] = l
        U=A
                
    return L,P,U


# In[32]:


def PALU(A,b):
    L,P,U= find_LPU(A)
    b = np.dot(P,b) 
    print("L\n",L)
    print("P\n",P)
    print("U\n",U)
    print(b)
    y = solve(L,b) #Ly=b
    print("y= ",y) 
    x=solve(U,y)   #Ux=y
    
    return x
    


# In[33]:


M= np.array([[2,1,5],[4,4,-4],[1,3,1]])
b=np.array([5,0,6])
PALU(M,b)


# In[29]:


#K=np.array([[7, 3, -1, 2], [3, 8, 1, -4], [-1, 1, 4, -1], [2, -4, -1, 6]])
#find_LPU(K)


# In[30]:


def cholesky(A):
    L = np.array([[0.0] * A.shape[0] for i in range(A.shape[0])])
   
    for i in range(A.shape[0]):
        for j in range(i+1):
            sum_ch = sum(L[i][k] * L[j][k] for k in range(j))
            
            if (i == j): 
                L[i][j] = math.sqrt(A[i][i] - sum_ch)
            else:
                L[i][j] = (1.0 / L[j][j] * (A[i][j] - sum_ch))
    return L
 


# In[31]:


A = np.array([[6, 3, 4, 8], [3, 6, 5, 1], [4, 5, 10, 7], [8, 1, 7, 25]])
L = cholesky(A)
print("A:\n",A)


print("L:\n",L)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




