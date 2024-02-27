#!/usr/bin/env python
# coding: utf-8

# In[1]:


from py_vollib.black_scholes import black_scholes as bs
from py_vollib.black_scholes.greeks.analytical import delta, gamma, vega, theta,rho


# In[ ]:





# In[2]:


import numpy as np
from scipy.stats import norm

r=0.01
S=30
K=40
T = 240/365
sigma = 0.30

def blackscholes(r,S,K,T,sigma,Type="c"):
    d1 = ((np.log(S/K))+(r+(sigma**2/2)*(T)))/(sigma*np.sqrt(T))
    d2 = d1-sigma*(np.sqrt(T))
    
    try:
        if Type == "c":
            price  = (S*norm.cdf(d1,0,1))-(K*np.exp(-r*T)*norm.cdf(d2,0,1))
        elif Type == "p":
            price = (K*np.exp(-r*T)*norm.cdf(-d2,0,1)) - (S*norm.cdf(-d1,0,1))
            
        return price ,bs(Type,S,K,T,r,sigma) 
    
    except:
        print("Please check the parameters 'p' or 'c' ")  




# In[3]:


print("option price: ", blackscholes(r,S,K,T,sigma,"c"))


# Delta

# In[4]:


def delta_calc(r,S,K,T,sigma,Type="c"):
    d1 = ((np.log(S/K))+(r+(sigma**2/2)*(T)))/(sigma*np.sqrt(T))
    d2 = d1-sigma*(np.sqrt(T))
    
    try:
        if Type == "c":
            delta_calc = norm.cdf(d1,0,1)
        elif Type == "p":
            delta_calc = -norm.cdf(-d1,0,1)
            
        return delta_calc ,delta(Type,S,K,T,r,sigma) 
    
    except:
        print("Please check the parameters 'p' or 'c' ")  


# In[5]:


print("the delta is: ",delta_calc(r,S,K,T,sigma,"c"))


# Gamma

# In[6]:


def gamma_calc(r,S,K,T,sigma,Type="c"):
    d1 = ((np.log(S/K))+(r+(sigma**2/2)*(T)))/(sigma*np.sqrt(T))
    d2 = d1-sigma*(np.sqrt(T))
    
    try:        
        gamma_calc = norm.pdf(d1,0,1)/(S*sigma*np.sqrt(T))   
        return gamma_calc ,gamma(Type,S,K,T,r,sigma) 
    
    except:
        print("Please check the parameters 'p' or 'c' ")  


# In[7]:


print("The gamma is :" ,gamma_calc(r,S,K,T,sigma,"c"))


# Vega

# In[15]:


def vega_calc(r,S,K,T,sigma,Type="c"):
    d1 = ((np.log(S/K))+(r+(sigma**2/2)*(T)))/(sigma*np.sqrt(T))
    d2 = d1-sigma*(np.sqrt(T))
    
    try:
        vega_calc  = (S*norm.pdf(d1,0,1)*np.sqrt(T))
        
        return vega_calc*0.01 ,vega(Type,S,K,T,r,sigma) 
    
    except:
        print("Please check the parameters 'p' or 'c' ")


# In[9]:


print("the vega is: ",vega_calc(r,S,K,T,sigma,Type="c"))


# Theta

# In[18]:


def theta_calc(r,S,K,T,sigma,Type="c"):
    d1 = ((np.log(S/K))+(r+(sigma**2/2)*(T)))/(sigma*np.sqrt(T))
    d2 = d1-sigma*(np.sqrt(T))
    
    try:
        if Type == "c":
            theta_calc  = (-S*norm.pdf(d1,0,1)*sigma/(2*np.sqrt(T)))-(r*K*np.exp(-r*T)*norm.cdf(d2,0,1))
        elif Type == "p":
            theta_calc =  (-S*norm.pdf(d1,0,1)*sigma/(2*np.sqrt(T)))+(r*K*np.exp(-r*T)*norm.cdf(d2,0,1))
            
        return theta_calc/365 ,theta(Type,S,K,T,r,sigma) 
    
    except:
        print("Please check the parameters 'p' or 'c' ")  


# In[11]:


print("the theta is : ",theta_calc(r,S,K,T,sigma,Type="c"))


# Rho

# In[16]:


def rho_calc(r,S,K,T,sigma,Type="c"):
    d1 = ((np.log(S/K))+(r+(sigma**2/2)*(T)))/(sigma*np.sqrt(T))
    d2 = d1-sigma*(np.sqrt(T))
    
    try:
        if Type == "c":
            rho_calc  = (K*T*np.exp(-r*T)*norm.cdf(d2,0,1))
        elif Type == "p":
            rho_calc = (-K*T*np.exp(-r*T)*norm.cdf(-d2,0,1))
            
        return rho_calc*0.01 ,rho(Type,S,K,T,r,sigma) 
    
    except:
        print("Please check the parameters 'p' or 'c' ")  


# All Together

# In[19]:


print("Option price: ",[round(x,3) for x in blackscholes(r,S,K,T,sigma,Type="c")])
print("Delta: ",[round(x,3) for x in delta_calc(r,S,K,T,sigma,Type="c")])
print("Gamma: ",[round(x,3) for x in gamma_calc(r,S,K,T,sigma,Type="c")])
print("Vega: ",[round(x,3) for x in vega_calc(r,S,K,T,sigma,Type="c")])
print("Theta: ",[round(x,3) for x in theta_calc(r,S,K,T,sigma,Type="c")])
print("Rho: ",[round(x,3) for x in rho_calc(r,S,K,T,sigma,Type="c")])


# In[ ]:





# In[ ]:





# In[ ]:




