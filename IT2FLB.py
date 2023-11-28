from pyit2fls import Mamdani, product_t_norm, max_s_norm, crisp, IT2FS, T1FS, tri_mf,ltri_mf, rtri_mf, trapezoid_mf, IT2FS_plot, T1Mamdani, gaussian_mf
import numpy as np
import math

def define_SIF():

    # AND NODE
    domain_freq = np.linspace(0., 22., 200)
    domain = np.linspace(0., 5., 200)

    U = IT2FS(domain_freq, ltri_mf, [14, 22, 1], ltri_mf, [15, 22, 1])
    R = IT2FS(domain_freq, tri_mf, [11, 13, 17, 1], tri_mf, [12, 13, 16, 1])
    VL = IT2FS(domain_freq, tri_mf, [9, 11, 14, 1], tri_mf, [10, 11, 13, 1])
    L = IT2FS(domain_freq, tri_mf, [6, 8, 11, 1], tri_mf, [7, 8, 10, 1])
    M = IT2FS(domain_freq, tri_mf, [3.5, 6, 9, 1], tri_mf, [4.5, 6, 8, 1])
    H = IT2FS(domain_freq, tri_mf, [1.5, 3, 7, 1], tri_mf, [2.5, 3, 6, 1])
    VH = IT2FS(domain_freq, rtri_mf, [3, 0, 1], rtri_mf, [2, 0, 1])

    #plot the membership functions
    #change the labels to match the ones in the paper
    IT2FS_plot(U, R, VL, L, M, H, VH)
    
    V = IT2FS(domain, trapezoid_mf, [3.5,4.5,5,6,1], trapezoid_mf, [3.75,4.5,5,6,1])
    IV = IT2FS(domain, tri_mf, [2.5, 3.5, 4.5, 1], tri_mf, [2.75, 3.5, 4.25, 1])
    III = IT2FS(domain, tri_mf, [1.5, 2.5, 3.5, 1], tri_mf, [1.75, 2.5, 3.25, 1])
    II = IT2FS(domain, tri_mf, [0.5, 1.5, 2.5, 1], tri_mf, [0.75, 1.5, 2.25, 1])
    I = IT2FS(domain, trapezoid_mf, [-3,0,0.5,1.5,1], trapezoid_mf, [-3,0,0.5,1.25,1])

    #plot the membership functions
    IT2FS_plot(V, IV, III, II, I)

    AR = U
    BR = R
    CR = VL
    DR = L
    ER = M
    FR = H
    GR = VH

    IT2FLS_AND = Mamdani(product_t_norm, max_s_norm, method="Centroid")
    IT2FLS_AND.add_input_variable("input1")
    IT2FLS_AND.add_input_variable("input2")
    IT2FLS_AND.add_output_variable("result")

    # RULES

    #G: Very High
    IT2FLS_AND.add_rule([("input1", VH), ("input2", I)], [("result", GR)])
    IT2FLS_AND.add_rule([("input1", VH), ("input2", II)], [("result", FR)])
    IT2FLS_AND.add_rule([("input1", VH), ("input2", III)], [("result", FR)])
    IT2FLS_AND.add_rule([("input1", VH), ("input2", IV)], [("result", FR)])
    IT2FLS_AND.add_rule([("input1", VH), ("input2", V)], [("result", FR)])

    #F: High
    IT2FLS_AND.add_rule([("input1", H), ("input2", I)], [("result", FR)])
    IT2FLS_AND.add_rule([("input1", H), ("input2", II)], [("result", FR)])
    IT2FLS_AND.add_rule([("input1", H), ("input2", III)], [("result", ER)])
    IT2FLS_AND.add_rule([("input1", H), ("input2", IV)], [("result", ER)])
    IT2FLS_AND.add_rule([("input1", H), ("input2", V)], [("result", DR)])

    #E: Moderate
    IT2FLS_AND.add_rule([("input1", M), ("input2", I)], [("result", ER)])
    IT2FLS_AND.add_rule([("input1", M), ("input2", II)], [("result", DR)])
    IT2FLS_AND.add_rule([("input1", M), ("input2", III)], [("result", DR)])
    IT2FLS_AND.add_rule([("input1", M), ("input2", IV)], [("result", CR)])
    IT2FLS_AND.add_rule([("input1", M), ("input2", V)], [("result", CR)])

    #D: Low
    IT2FLS_AND.add_rule([("input1", L), ("input2", I)], [("result", DR)])
    IT2FLS_AND.add_rule([("input1", L), ("input2", II)], [("result", DR)])
    IT2FLS_AND.add_rule([("input1", L), ("input2", III)], [("result", CR)])
    IT2FLS_AND.add_rule([("input1", L), ("input2", IV)], [("result", CR)])
    IT2FLS_AND.add_rule([("input1", L), ("input2", V)], [("result", BR)])

    #C: Very Low
    IT2FLS_AND.add_rule([("input1", VL), ("input2", I)], [("result", CR)])
    IT2FLS_AND.add_rule([("input1", VL), ("input2", II)], [("result", BR)])
    IT2FLS_AND.add_rule([("input1", VL), ("input2", III)], [("result", BR)])
    IT2FLS_AND.add_rule([("input1", VL), ("input2", IV)], [("result", BR)])
    IT2FLS_AND.add_rule([("input1", VL), ("input2", V)], [("result", AR)])

    #B: Unlikely
    IT2FLS_AND.add_rule([("input1", R), ("input2", I)], [("result", BR)])
    IT2FLS_AND.add_rule([("input1", R), ("input2", II)], [("result", AR)])
    IT2FLS_AND.add_rule([("input1", R), ("input2", III)], [("result", AR)])
    IT2FLS_AND.add_rule([("input1", R), ("input2", IV)], [("result", AR)])
    IT2FLS_AND.add_rule([("input1", R), ("input2", V)], [("result", AR)])

    #A: Remote
    IT2FLS_AND.add_rule([("input1", U), ("input2", I)], [("result", AR)])
    IT2FLS_AND.add_rule([("input1", U), ("input2", II)], [("result", AR)])
    IT2FLS_AND.add_rule([("input1", U), ("input2", III)], [("result", AR)])
    IT2FLS_AND.add_rule([("input1", U), ("input2", IV)], [("result", AR)])
    IT2FLS_AND.add_rule([("input1", U), ("input2", V)], [("result", AR)])

    # OR NODE

    IT2FLS_OR = Mamdani(product_t_norm, max_s_norm, method="Centroid")
    IT2FLS_OR.add_input_variable("input1")
    IT2FLS_OR.add_input_variable("input2")
    IT2FLS_OR.add_output_variable("result")

    # RULES

    #G: Very High
    IT2FLS_OR.add_rule([("input1", VH), ("input2", VH)], [("result", GR)])
    IT2FLS_OR.add_rule([("input1", VH), ("input2", H)], [("result", GR)])
    IT2FLS_OR.add_rule([("input1", VH), ("input2", M)], [("result", GR)])
    IT2FLS_OR.add_rule([("input1", VH), ("input2", L)], [("result", GR)])
    IT2FLS_OR.add_rule([("input1", VH), ("input2", VL)], [("result", GR)])
    IT2FLS_OR.add_rule([("input1", VH), ("input2", R)], [("result", GR)])
    IT2FLS_OR.add_rule([("input1", VH), ("input2", U)], [("result", GR)])

    #F: High
    IT2FLS_OR.add_rule([("input1", H), ("input2", VH)], [("result", GR)])
    IT2FLS_OR.add_rule([("input1", H), ("input2", H)], [("result", FR)])
    IT2FLS_OR.add_rule([("input1", H), ("input2", M)], [("result", FR)])
    IT2FLS_OR.add_rule([("input1", H), ("input2", L)], [("result", FR)])
    IT2FLS_OR.add_rule([("input1", H), ("input2", VL)], [("result", FR)])
    IT2FLS_OR.add_rule([("input1", H), ("input2", R)], [("result", FR)])
    IT2FLS_OR.add_rule([("input1", H), ("input2", U)], [("result", FR)])

    #E: Moderate
    IT2FLS_OR.add_rule([("input1", M), ("input2", VH)], [("result", GR)])
    IT2FLS_OR.add_rule([("input1", M), ("input2", H)], [("result", FR)])
    IT2FLS_OR.add_rule([("input1", M), ("input2", M)], [("result", ER)])
    IT2FLS_OR.add_rule([("input1", M), ("input2", L)], [("result", ER)])
    IT2FLS_OR.add_rule([("input1", M), ("input2", VL)], [("result", ER)])
    IT2FLS_OR.add_rule([("input1", M), ("input2", R)], [("result", ER)])
    IT2FLS_OR.add_rule([("input1", M), ("input2", U)], [("result", ER)])

    #D: Low
    IT2FLS_OR.add_rule([("input1", L), ("input2", VH)], [("result", GR)])
    IT2FLS_OR.add_rule([("input1", L), ("input2", H)], [("result", FR)])
    IT2FLS_OR.add_rule([("input1", L), ("input2", M)], [("result", ER)])
    IT2FLS_OR.add_rule([("input1", L), ("input2", L)], [("result", DR)])
    IT2FLS_OR.add_rule([("input1", L), ("input2", VL)], [("result", DR)])
    IT2FLS_OR.add_rule([("input1", L), ("input2", R)], [("result", DR)])
    IT2FLS_OR.add_rule([("input1", L), ("input2", U)], [("result", DR)])

    #C: Very Low
    IT2FLS_OR.add_rule([("input1", VL), ("input2", VH)], [("result", GR)])
    IT2FLS_OR.add_rule([("input1", VL), ("input2", H)], [("result", FR)])
    IT2FLS_OR.add_rule([("input1", VL), ("input2", M)], [("result", ER)])
    IT2FLS_OR.add_rule([("input1", VL), ("input2", L)], [("result", DR)])
    IT2FLS_OR.add_rule([("input1", VL), ("input2", VL)], [("result", CR)])
    IT2FLS_OR.add_rule([("input1", VL), ("input2", R)], [("result", CR)])
    IT2FLS_OR.add_rule([("input1", VL), ("input2", U)], [("result", CR)])

    #B: Unlikely
    IT2FLS_OR.add_rule([("input1", R), ("input2", VH)], [("result", GR)])
    IT2FLS_OR.add_rule([("input1", R), ("input2", H)], [("result", FR)])
    IT2FLS_OR.add_rule([("input1", R), ("input2", M)], [("result", ER)])
    IT2FLS_OR.add_rule([("input1", R), ("input2", L)], [("result", DR)])
    IT2FLS_OR.add_rule([("input1", R), ("input2", VL)], [("result", CR)])
    IT2FLS_OR.add_rule([("input1", R), ("input2", R)], [("result", BR)])
    IT2FLS_OR.add_rule([("input1", R), ("input2", U)], [("result", BR)])

    #A: Remote
    IT2FLS_OR.add_rule([("input1", U), ("input2", VH)], [("result", GR)])
    IT2FLS_OR.add_rule([("input1", U), ("input2", H)], [("result", FR)])
    IT2FLS_OR.add_rule([("input1", U), ("input2", M)], [("result", ER)])
    IT2FLS_OR.add_rule([("input1", U), ("input2", L)], [("result", DR)])
    IT2FLS_OR.add_rule([("input1", U), ("input2", VL)], [("result", CR)])
    IT2FLS_OR.add_rule([("input1", U), ("input2", R)], [("result", BR)])
    IT2FLS_OR.add_rule([("input1", U), ("input2", U)], [("result", AR)])

    return IT2FLS_AND, IT2FLS_OR

def find_out_word(x, xmf):  # Adjetivo a partir de x(valor) e xmf(funcao de pert)
        aux = {}
        for word in xmf.terms.keys():
            aux[word] = fuzz.interp_membership(
                xmf.universe, xmf.terms[word].mf, x)
        word = list(aux.keys())[np.argmax(list(aux.values()))]
        return word

def fuzzy_bowtie(f_treats, pfd):

    #Define SIF
    FSand,FSor = define_SIF()
    pathway = []
    
    # ### Organização da matriz com os valores de frequência

    #criar uma lista auxiliar à f_list com 23 elementos
    pfd_aux = [0]*15

    pfd_aux[0] = pfd[0] #Chemical treatment
    pfd_aux[1] = pfd[1] #Corrosion monitoring
    pfd_aux[2] = pfd[2] #Steel containment envelope
    pfd_aux[3] = pfd[3] #Corrosion protection system
    pfd_aux[4] = pfd[4] #Operation conditions monitoring
    pfd_aux[5] = pfd[5] #Filtering
    pfd_aux[6] = pfd[6] #Particle detection
    pfd_aux[7] = pfd[7] #Vibration and tension monitoring system
    pfd_aux[8] = pfd[8] #PSH, LSH
    pfd_aux[9] = pfd[9] #SIS
    pfd_aux[10] = pfd[10] #Relief devices
    pfd_aux[11] = pfd[11] #Pipeline barrier
    pfd_aux[12] = pfd[2] #Steel containment envelope
    pfd_aux[13] = pfd[12] #Checklist
    pfd_aux[14] = pfd[7] #Vibration and tension monitoring system

    # print(pfd_aux)

    ### Separação barreiras lado esquerdo e lado direito
    path_len_1 = [2,2,3,1,3,2,2]
    ini = 0

    #Definição do Bow-tie

    for i in range(len(f_treats)):
        result = -math.loVH0(f_treats[i])
        fim = ini+path_len_1[i]
        for j in pfd_aux[ini:fim]:
            input = -math.loVH0(j)
            output = FSand.evaluate({'input1': result, 'input2': input})
            result = crisp(output[1]['result'])
        pathway.append(result)
        ini = fim
    
    Ftop = pathway[0]
    for i in pathway[1:]:
        output = FSor.evaluate({'input1': Ftop, 'input2': i})
        Ftop = crisp(output[1]['result'])
    
    return Ftop

def model_bowtie(sample):

    treats=[0.01,0.01,0.01,0.01,0.03,0.01,0.02]

    sample = np.array(sample,ndmin=2)
    Y = np.empty((sample.shape[0],1))

    #Caso todas as amostras sejam passadas em um único array
    sample = sample.tolist()

    for k in range(len(sample)):
        print(f'Amostra {k+1} de {len(sample)}')
        Y[k] = 10**-fuzzy_bowtie(treats, sample[k])

    return Y
