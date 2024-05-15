import numpy as np

import matplotlib.pyplot as plt

import argparse


parser = argparse.ArgumentParser()

parser.add_argument('-t', '--tag', default='profile_check', help='txt file name generated from previous charmm command')
parser.add_argument('-n','--numPoints',default=4, help='number of points in ion channel')
parser.add_argument('-f','--frames',default=3900, help='number of frames after equilibration')

parser.add_argument('-c','--check_profile', action='store_true')
parser.add_argument('--no-profile', dest='feature', action='store_false')
parser.set_defaults(check_profile=False)

args = parser.parse_args()

s=int(args.numPoints)
fr=int(args.frames)

if not args.check_profile:
    file = 'output/observables/barrier/channel/ionchan'+args.tag+'.txt'
    with open(file) as fi:
        f=fi.readlines()

    nrgs = np.zeros(shape=(fr,s))

    for l in range(len(f)):
        frame=int(np.floor(l/s))
        aa=int(l%(s+1))
        if frame<fr:
            #print(aa,f[l])
            if aa!=s:
                nrgs[frame,aa]=float(f[l])
            #else:
                #print('Frame '+str(frame))        

    # Average over the points in the barrier domain
    barrier = np.zeros(fr)
    for i in range(fr):
        barrier[i] = np.average(nrgs[i])

    np.save('output/observables/barrier/channel/barrier'+args.tag+'.npy',barrier)

elif args.check_profile:
    file = 'output/observables/barrier/channel/'+args.tag+'.txt'
    with open(file) as fi:
        f=fi.readlines()

    nrgs = np.zeros(shape=(fr,s))
    for l in range(len(f)):
        frame=int(np.floor(l/s))
        aa=int(l%s)
        if frame<fr:
            #print(aa,f[l])
            if aa!=s-1:
                #print(aa,s-1)
                nrgs[frame,aa]=float(f[l])

    plt.plot(nrgs[0])
    plt.show()