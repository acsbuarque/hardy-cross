from deap import tools, algorithms, base, creator
import numpy as np
import pandas as pd
import random
import time
from pyeasyga import pyeasyga

#Dados
dados = pd.read_csv('dados.csv', encoding = "ISO-8859-1", index_col='Nó')
dados_observados = pd.read_csv('dados_observados_def.csv', index_col='Vazão (l/s)')
dados_trechos = pd.read_csv('dados_trechos_def.csv', index_col='Trecho')
Nivel_Reservatorio = 240
g = 9.8
vazoes_reservatorio = np.array([35.62,28.11,29.24,64.63,79.51,90.07,79.73,74.24,71.5,82.17,66.3,50.86])/1000

#Inserindo vazões aleatórias em cada anel
#vazoes_observadas = dados_observados['Vazão (l/s)'].as_matrix()*1000
#Separando os dados de cada anel
anel1, anel2, anel3 = dados_trechos.iloc[[1,2,15,14,16,17,18,19]], dados_trechos.iloc[[1,2,4,5,6,7,8]], dados_trechos.iloc[[3,9,10,11,12,13,14,15]] 
Q1  = dados_trechos.iloc
#E = 0.0015

def Reynolds(Q,D):
    return abs(Q*4*(10**6)/(np.pi*D/1000))

def fator_de_atrito(E,D,Rey):
    return 0.25/(np.log10(abs((E/(3.7*D))+(5.74/(Rey**0.9))))**2)

def perda_de_carga(f,L,Q,D):
    V = 4*Q/(np.pi*((D/1000)**2))
    return f*L*(V**2)/((D/1000)*2*g)

def hardy_cross(E,Q):
    """Recebe a rugosidade do material e a vazão no reservatório e retorna a pressão nos
    pontos 5, 12 e 16"""
    somaH = np.full(3,10)
    ReyR1 = Reynolds(Q,dados_trechos.loc[['R-1'],'Diametro (mm)']) #Reynolds do trecho R-1
    FR1 = fator_de_atrito(E,dados_trechos.loc[['R-1'],'Diametro (mm)'],ReyR1) #Fator de atrito do trecho R-1
    H1 = perda_de_carga(FR1,dados_trechos.loc[['R-1'],'Comprimento'],Q,dados_trechos.loc[['R-1'],'Diametro (mm)'])
    perda1 = pd.Series(H1, index=['R-1'])
    
    start = time.time()
    PERIOD_OF_TIME = 60

    while np.all(somaH>np.full(3,0.01)):
        dados_trechos['coeficientes'] = (dados_trechos['vazoes'].values)/abs((dados_trechos['vazoes'].values))
        
        #Calculando a perda de carga no Anel 1
        anel1['Reynolds']=Reynolds(anel1['vazoes'],anel1['Diametro (mm)'])
        anel1['Fator de Atrito']=fator_de_atrito(E,anel1['Diametro (mm)'],anel1['Reynolds'])
        anel1['Perda de Carga']=perda_de_carga(anel1['Fator de Atrito'],anel1['Comprimento'],anel1['vazoes'],anel1['Diametro (mm)'])*anel1['coeficientes']
        somaH1 = anel1['Perda de Carga'].sum()
            
        #Calculando a perda de carga no Anel 2
        anel2['Reynolds']=Reynolds(anel2['vazoes'],anel2['Diametro (mm)'])
        anel2['Fator de Atrito']=fator_de_atrito(E,anel2['Diametro (mm)'],anel2['Reynolds'])
        anel2['Perda de Carga']=perda_de_carga(anel2['Fator de Atrito'],anel2['Comprimento'],anel2['vazoes'],anel2['Diametro (mm)'])*anel2['coeficientes']
        somaH2 = anel2['Perda de Carga'].sum()
        
        #Calculando a perda de carga no Anel 3
        anel3['Reynolds']=Reynolds(anel3['vazoes'],anel3['Diametro (mm)'])
        anel3['Fator de Atrito']=fator_de_atrito(E,anel3['Diametro (mm)'],anel3['Reynolds'])
        anel3['Perda de Carga']=perda_de_carga(anel3['Fator de Atrito'],anel3['Comprimento'],anel3['vazoes'],anel3['Diametro (mm)'])*anel3['coeficientes']
        somaH3 = anel3['Perda de Carga'].sum()
    
        #Calculando as vazões demandadas no ponto
        Q_demandas = (dados['Fração demanda (%)'].values)*Q/100
        
        #Calculando o fator de correção    
        somaH1Q1 = np.sum(np.divide((anel1['Perda de Carga'].values),(anel1['vazoes'].values)))
        somaH2Q2 = np.sum(np.divide((anel2['Perda de Carga'].values),(anel2['vazoes'].values)))
        somaH3Q3 = np.sum(np.divide((anel3['Perda de Carga'].values),(anel3['vazoes'].values)))
        
        #Organizando o somatório das perdas por anel
        somaH = np.array([somaH1,somaH2,somaH3])
        
        Qa1 = anel1['vazoes'].sum()
        Qa2 = anel2['vazoes'].sum()
        Qa3 = anel3['vazoes'].sum()        
        
        DeltaQ1 = -somaH1/(17*somaH1Q1)
        DeltaQ2 = -somaH2/(17*somaH2Q2)
        DeltaQ3 = -somaH3/(17*somaH3Q3)
        
        #Atualizando os dados de vazão               
        dados_trechos.loc[['1-2','2-3'],['vazoes']]=(dados_trechos.loc[['1-2','2-3'],['vazoes']].values)+DeltaQ1+DeltaQ2 #delta1+2
        dados_trechos.loc[['4-3','5-4'],['vazoes']]=(dados_trechos.loc[['4-3','5-4'],['vazoes']].values)+DeltaQ2+DeltaQ3 #delta2+3
        dados_trechos.loc[['5-6','6-7','7-8','8-1'],['vazoes']] = (dados_trechos.loc[['5-6','6-7','7-8','8-1'],['vazoes']].values)+DeltaQ1 #delta1
        dados_trechos.loc[['3-13','14-15','15-16','16-17','17-5'],['vazoes']] = (dados_trechos.loc[['3-13','14-15','15-16','16-17','17-5'],['vazoes']].values)+DeltaQ2 #delta2
        dados_trechos.loc[['9-1','10-9','11-10','12-11','3-12'],['vazoes']] = (dados_trechos.loc[['9-1','10-9','11-10','12-11','3-12'],['vazoes']].values)+DeltaQ3 #delta3
        
        anel1.update(dados_trechos)
        anel2.update(dados_trechos)
        anel3.update(dados_trechos)
        
        if time.time() > start + PERIOD_OF_TIME : break
    p = pd.concat([perda1, anel1['vazoes'],anel2['vazoes'],anel3['vazoes']])
    p = p[~p.index.duplicated(keep='first')]
    no5 = p.loc[['R-1','8-1','7-8','6-7','5-6']]
    no12 = p.loc[['R-1','9-1','10-9','11-10','5-6','12-11']]
    no16 = p.loc[['R-1','8-1','7-8','6-7','5-6','17-5','16-17']]
    
    #Perda de carga em cada nó
    H_no5 = no5.sum()
    H_no12 = no12.sum()
    H_no16  = no16.sum()
    
    #Pressão em cada nó
    P_no5 = (Nivel_Reservatorio - dados.loc[[5],'Cota topografica (m)'].get(5) - H_no5)*9800
    P_no12 = (Nivel_Reservatorio - dados.loc[[12],'Cota topografica (m)'].get(12) - H_no12)*9800
    P_no16 = (Nivel_Reservatorio - dados.loc[[16],'Cota topografica (m)'].get(16) - H_no16)*9800
        
    return (P_no5,P_no12,P_no16)


individuos = np.random.uniform(0.0015, 0.1, 10)
data = []
c = 0
for n in individuos:
    for q in vazoes_reservatorio:
        dados_trechos['vazoes']=np.array([35.62,10.74,9.65,5.57,-0.91,2.57,4.2,5.25,7,2.67,0.35,-4.23,-6.66,-7.35,-0.2,-1,-9.63,-11.63,-13.59,-15])/1000
        c = c + 1
        i = {}
        i['vazao'] = q
        i['rugosidade'] = n
        i['P5'],i['P12'],i['P16'] = hardy_cross(n,q)
        data.append(i)
    
ga = pyeasyga.GeneticAlgorithm(data)

def fitness(individual, data):
    for (selected, item) in zip(individual, data):
        #solucao = {}
        erro = (item['P5'] - (dados_observados.loc[[str(item['vazao'])],'pressao5'].get(0))**2/1000) + (item['P12'] - (dados_observados.loc[[str(item['vazao'])],'pressao12'].get(0))**2/1000) + (item['P16'] - (dados_observados.loc[[str(item['vazao'])],'pressao16'].get(0))**2/1000)
        if erro > 0.001:
            erro = 0
    return erro

ga.fitness_function = fitness               # set the GA's fitness function
ga.run()                                    # run the GA
valor, vetor = ga.best_individual()            # print the GA's best solution

solucoes_binario = np.asarray(vetor)
total_solucoes = np.multiply(np.arange(120),solucoes)
solucoes_validas = total_solucoes[total_solucoes != 0]
solucoes_finais = pd.DataFrame()
for i in solucoes_validas:
    solucoes_finais.append(i)



