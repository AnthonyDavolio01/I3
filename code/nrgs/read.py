import numpy as np
import matplotlib.pyplot as plt

cat=['0', '1a', '1b', '1c', '1d', '1e', '2ab', '2ac', '2ad', '2ae', '2bc', '2bd', '2be', '2cd', '2ce', '2de', '3abc', '3abd', '3abe', '3acd', '3ace', '3ade', '3bcd', '3bce', '3cde', '4abcd', '4abce', '4abde', '4acde', '4bcde', '5']

for tag in cat:
    print(tag)
    file = 'output/observables/nrgs/nrgs'+tag+'.txt'
    with open(file) as fi:
        f=fi.readlines()



    na=488
    ns=5
    nf=4000

    nrgs = np.zeros(shape=(na*ns,nf))


    r=[]
    fr=0
    for l in f:
        if l=='___\n':
            nrgs[:,fr]=r
            r=[]
            print('Frame '+str(fr))
            fr+=1
        else:
            r.append(l)

    np.save('output/observables/nrgs/nrgs'+tag+'.npy',nrgs)



#    for l in range(len(f)):
#        frame=int(np.floor(l/(na*ns)))
#        aa=int(l%(na*ns))
#        if aa!=488:
#            nrgs[aa,frame]=f[l]
#        else:
#            print('Frame '+str(frame))
