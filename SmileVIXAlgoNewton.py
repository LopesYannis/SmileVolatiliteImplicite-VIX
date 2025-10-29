import matplotlib.pyplot as plt
import random as rd
import pandas as pd
import numpy as np
import math

link="H:\Desktop\python\Ing3\TD1_SmileVolatilite\sp-index.txt"
fichier=np.loadtxt(link, skiprows=1)
N=len(fichier)
T=np.zeros(N)
        
def d1(S_0,K,r,sigma,T,t):
    #valeur d1 formule BS
    d1=(np.log(S_0/K)+(r+sigma**2/2)*(T-t))/(sigma*np.sqrt(T-t));
    return d1

def d2(S_0,K,r,sigma,T,t):
    #valeur d2 formule BS
    d2=(np.log(S_0/K)+(r-sigma**2/2)*(T-t))/(sigma*np.sqrt(T-t));
    return d2

def N(x):
    #fonction de répartition d'une loi normale centrée réduite
    return (1/2)*(1+math.erf(x/np.sqrt(2)))

def Valeur_BScall(S_0,K,r,sigma,T,t):
    #CallEU Black&Scholes
    if (T-t)==0:
        V=np.max(S_0-K,0)
    else:
        V=S_0*N(d1(S_0,K,r,sigma,T,t))-K*np.exp(-r*(T-t))*N(d2(S_0,K,r,sigma,T,t))
    return V

def Valeur_BSput(S_0,K,r,sigma,T,t):
    #PutEU Black&Scholes
    if (T-t)==0:
        V=np.max(K-S_0,0)
    else:
        V=K*np.exp(-r*(T-t))*N(-d2(S_0,K,r,sigma,T,t))-S_0*N(-d1(S_0,K,r,sigma,T,t))
    return V

def Vega_BS(S_0,K,r,sigma,T,t):
    #Vega Black&Scholes Call&PutEU
    Vega=S_0*np.sqrt(T/(2*np.pi))*np.exp((-d1(S_0,K,r,sigma,T,t)**2)/2)
    return Vega

def graph():
    sigmaTEST = 0.5
    TTEST = 0.5
    rTEST = 0.1
    KTEST = 10
    tTEST = 0 

    #abcisse S0
    S0_vals = np.linspace(0.01, 20, 200)  # éviter S0=0 pour log

    #calcul ordonnée
    valeurs = [Valeur_BScall(i, KTEST, rTEST, sigmaTEST, TTEST, tTEST) for i in S0_vals]
    vegas = [Vega_BS(i, KTEST, rTEST, sigmaTEST, TTEST, tTEST) for i in S0_vals]

    # Tracé
    plt.figure(figsize=(10,5))

    plt.subplot(1,2,1)
    plt.plot(S0_vals, valeurs, label="Valeur BS (Call)")
    plt.xlabel("S0")
    plt.ylabel("Valeur de l'option")
    plt.title("Valeur du call en fonction de S0")
    plt.grid(False)
    plt.legend()
    
    plt.subplot(1,2,2)
    plt.plot(S0_vals, vegas, label="Vega", color="orange")
    plt.xlabel("S0")
    plt.ylabel("Vega")
    plt.title("Vega en fonction de S0")
    plt.grid(False)
    plt.legend()
    
    plt.tight_layout()
    plt.show()

def AlgoNewton(S_0,T,r,K,p,nu):
    #calcule la volatilité implicite pour chaque strike pour un CallEU
    sigma0=np.sqrt(2*np.abs((np.log(S_0/K)+r*T)/T))
    sigma0 = max(0.0001, min(sigma0, 5)) #ajout pour ne pas avoir un sigma qui explose
    
    F=Valeur_BScall(S_0,K,r,sigma0,T,0)-p
    dF=S_0*np.sqrt((T)/2*np.pi)*np.exp((-(d1(S_0,K,r,sigma0,T,0)**2))/2)
    while (np.abs(F)>nu) and (np.abs(dF) > 1e-8):
        sigma0=sigma0-(F/dF)
        F=Valeur_BScall(S_0,K,r,sigma0,T,0)-p
        dF=S_0*np.sqrt((T)/2*np.pi)*np.exp((-(d1(S_0,K,r,sigma0,T,0)**2))/2)
    return sigma0

def AlgoNewtonPut(S_0,T,r,K,p,nu):
    #calcule la volatilité implicite pour chaque strike pour un put EU
    sigma0=np.sqrt(2*np.abs((np.log(S_0/K)+r*T)/T))
    sigma0 = max(0.0001, min(sigma0, 5)) #ajout pour ne pas avoir un sigma qui explose
    
    F=Valeur_BSput(S_0,K,r,sigma0,T,0)-p
    dF=S_0*np.sqrt((T)/2*np.pi)*np.exp((-(d1(S_0,K,r,sigma0,T,0)**2))/2)
    while (np.abs(F)>nu) and (np.abs(dF) > 1e-8):
        sigma0=sigma0-(F/dF)
        F=Valeur_BSput(S_0,K,r,sigma0,T,0)-p
        dF=S_0*np.sqrt((T)/2*np.pi)*np.exp((-(d1(S_0,K,r,sigma0,T,0)**2))/2)
    return sigma0

def newtonTest():
    S_0=5430.3
    T=4/12
    r=0.05    
    K=[5125, 5225, 5325, 5425, 5525, 5625, 5725, 5825]
    p=[475, 405, 340, 280.5, 226, 179.5, 139, 105]
    nu=0.01
    V=[] 

    for i in range(len(K)):
        vol=AlgoNewton(S_0,T,r,K[i],p[i],nu)
        V.append(vol)
        
    plt.figure(figsize=(8,5))
    plt.plot(K, V, marker="o", linestyle="-", color="b", label="Vol implicite")
    plt.xlabel("Strike (K)")
    plt.ylabel("Volatilité implicite")
    plt.title("Smile de volatilité implicite")
    plt.grid(False)
    plt.legend()
    plt.show()

def GraphVolImpliciteSPX():
    q=0.0217
    eps=0.001
    t=0
    S_0=1260
    link="H:\Desktop\python\Ing3\TD1_SmileVolatilite\sp-index.txt"
    fichier=np.loadtxt(link, skiprows=1)
    N=len(fichier)
    K=np.zeros(N)
    r=np.zeros(N)
    T=np.zeros(N)
    S_0L=np.zeros(N)
    vol_imp=np.zeros(N)
    vol_imp2=np.zeros(N)
    Ca=np.zeros(N)
    Cb=np.zeros(N)
    Pb=np.zeros(N)
    Pa=np.zeros(N)
    pC=np.zeros(N)
    pP=np.zeros(N)
    Nindice=[58,115,149,171,214,237,265]
    for i in range(N):
        T[i]=fichier[i][0]
        K[i]=fichier[i][1]
        Cb[i]=fichier[i][2]
        Ca[i]=fichier[i][3]
        Pb[i]=fichier[i][4]
        Pa[i]=fichier[i][5]
        r[i]=fichier[i][6]
        S_0L[i]=S_0*np.exp(-q*T[i])
    for i in range(N):
        pC[i]=(Ca[i]+Cb[i])/2
        if pC[i]<S_0L[i] and pC[i]>=np.max(S_0L[i]-K[i]*np.exp(-(r[i]/100)*T[i]),0):
            vol_imp[i]=AlgoNewton(S_0L[i],T[i],r[i]/100,K[i],pC[i],eps)
        else :
            vol_imp[i]=0
    for i in range(N):
        pP[i]=(Pa[i]+Pb[i])/2
        if pP[i]<K[i]*np.exp(-(r[i]/100)*T[i]) and pP[i]>np.max(K[i]*np.exp(-(r[i]/100)*T[i])-S_0L[i],0):
            vol_imp2[i]=AlgoNewtonPut(S_0L[i],T[i],r[i]/100,K[i],pP[i],eps)
        else:
            vol_imp2[i]=0
    vol_impBis=[]
    vol_imp2Bis=[]
    Kbis=[]
    Tbis=[]
    Kbis2=[]
    Tbis2=[]
    for i in range(N):
        if vol_imp[i]!=0:
            vol_impBis.append(vol_imp[i])
            Kbis.append(K[i])
            Tbis.append(T[i])
        if vol_imp2[i]!=0:
            vol_imp2Bis.append(vol_imp2[i])
            Kbis2.append(K[i])
            Tbis2.append(T[i])
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(Kbis, Tbis, vol_impBis, c='m', label="Vol implicite Call")
    ax.scatter(Kbis2, Tbis2, vol_imp2Bis, c='c', label="Vol implicite Put")
    ax.set_title("Smile volatilité implicite en 3D")
    ax.set_xlabel('K (Strike)')
    ax.set_ylabel('T (Maturité)')
    ax.set_zlabel('Volatilité implicite')
    ax.legend()
    plt.tight_layout()
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(Kbis, Tbis, vol_impBis, c='m', label="Vol implicite Call")
    ax.set_title("Smile volatilité Call implicite en 3D")
    ax.set_xlabel('K (Strike)')
    ax.set_ylabel('T (Maturité)')
    ax.set_zlabel('Volatilité implicite')
    ax.legend()
    plt.tight_layout()
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(Kbis2, Tbis2, vol_imp2Bis, c='c', label="Vol implicite Put")
    ax.set_title("Smile volatilité Put implicite en 3D")
    ax.set_xlabel('K (Strike)')
    ax.set_ylabel('T (Maturité)')
    ax.set_zlabel('Volatilité implicite')
    ax.legend()
    plt.tight_layout()


    plt.show()
