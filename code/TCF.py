import numpy as np
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-p", "--plot", default=3,type=int, help="0 for no plots, 1 for individual plots, 2 for combined plots")
parser.add_argument("-s", "--start", default=100,type=int, help="starting snapshot")
parser.add_argument("-e", "--end", default=4000,type=int, help="ending snapshot")
parser.add_argument("-f", "--pdf", default='gmm',type=str, help="probability density function estimation method")
args = vars(parser.parse_args())
plotwhat=args["plot"]
snapStart=args["start"]
snapEnd=args["end"]
pdfmode=args["pdf"]

nf=snapEnd-snapStart

cat=['0', '1a', '1b', '1c', '1d', '1e', '2ab', '2ac', '2ad', '2ae', '2bc', '2bd', '2be', '2cd', '2ce', '2de', '3abc', '3abd', '3abe', '3acd', '3ace', '3ade', '3bcd', '3bce', '3cde', '4abcd', '4abce', '4abde', '4acde', '4bcde', '5']

actcv=np.zeros(nf*len(cat))
ligcv=np.zeros(nf*len(cat))
nrg=np.zeros(shape=(488*5,nf*len(cat)))

count=0
for tag in cat:
    barrier=np.load('output/observables/barrier/channel/barrier'+tag+'.npy')
    actcv[count*nf:(count+1)*nf]=barrier#[snapStart:snapEnd]
    nLig=int(tag[0])
    ligcv[count*nf:(count+1)*nf]=nLig
    
    nrgs_temp=np.load('output/observables/nrgs/nrgs'+tag+'.npy')
    nrg[:,count*nf:(count+1)*nf]=nrgs_temp[:,snapStart:snapEnd]

    count+=1

aCV=[]
bCV=[]
nrg_E=[]
for i in range(len(actcv)):
    if actcv[i]!=0:
        aCV.append(actcv[i])
        bCV.append(ligcv[i])
        nrg_E.append(nrg[:,i])
#print(nrg)

nrg_E=np.array(nrg_E).T
bCV=np.array(bCV)
aCV=np.array(aCV)
nrgs_E=(nrg_E.T)#-nrg_E.T.mean(axis=0)).T
nrgs = np.exp(-nrgs_E/(0.00198721319*310))

if False: # Make heatmap
    heatmap, xedges, yedges = np.histogram2d(bCV, aCV, bins=20)

    xmesh, ymesh = np.meshgrid(xedges[:-1], yedges[:-1])

    plt.pcolormesh(xmesh, ymesh, heatmap.T, cmap='viridis')
    plt.colorbar(label='# of snapshots')

    plt.xlabel('# of Ligands Bound')
    plt.ylabel('Channel Energy Barrier (KCal / Mol)')
    plt.title('Energy Barrier against # of Ligands Bound Heat Density Map')

    plt.show()
    exit()

numAminos=488*5
numPoints=len(aCV)
if pdfmode=='gmm':
    num_components=2 ### Need to update other numbers if this is changed
    # Fit a Gaussian Mixture Model with two components
    b_gmm = GaussianMixture(n_components=6).fit(bCV.reshape(-1, 1))
    a_gmm = GaussianMixture(n_components=2).fit(aCV.reshape(-1, 1))
    a_p=np.rint(a_gmm.predict_proba(aCV.reshape(-1, 1)))
    b_p=np.rint(b_gmm.predict_proba(bCV.reshape(-1, 1)))
    anb_p=[]
    for sn in range(numPoints):
        anb_p.append(np.outer(a_p[sn],b_p[sn]))
    anb_p=np.array(anb_p)
    
    # 2 bins for activation, 6 bins for ligand
    pa=np.zeros(2)
    pb=np.zeros(6)
    pab=np.zeros(shape=(2,6))
    
    for i in range(2):
        pa[i]=np.sum(a_p[:,i]/numPoints)
        for j in range(6):
            pb[j]=np.sum(b_p[:,j]/numPoints)
            pab[i,j]=np.sum(anb_p[:,i,j]/numPoints)
    p=1/numPoints

    tot_anb=np.zeros(shape=(2,6))
    for sn in range(numPoints):
        tot_anb=tot_anb+anb_p[sn]
    #print(tot_anb)
    #exit()
        
    sfa=np.zeros(numAminos)
    sfb=np.zeros(numAminos)
    sfab=np.zeros(numAminos)
    sf=np.zeros(numAminos)
    for sn in range(numPoints):
        sfa=sfa+nrgs[sn]/np.sum(np.dot(a_p[sn],tot_anb))#*np.dot(pa,a_p[sn])
        sfb=sfb+nrgs[sn]/np.sum(np.dot(tot_anb,b_p[sn]))#*np.dot(pb,b_p[sn])
        sfab=sfab+nrgs[sn]/np.dot(anb_p[sn].flatten(),tot_anb.flatten())#*np.dot(pab.flatten(),anb_p[sn].flatten())
        sf=sf+nrgs[sn]/np.sum(tot_anb.flatten())#*p
        #print(sfa, sfb, sfab, sf, nrgs[sn])
    #i3=-(1/0.00198721319)*np.log(sfab*sf/(sfa*sfb))
    #i3=sfab*sf/(sfa*sfb)
    
    i3=np.log(sfab*sf/(sfa*sfb))


tcf=[]
for i in range(len(i3)):
    if i%5==0:
        tcf.append(i3[i])

i3=np.array(tcf)
i3=np.nan_to_num(i3,nan=0)

np.save("output/i3.npy",i3)
print(i3)
min_val = i3.min()
max_val = i3.max()
print(min_val)
print(max_val)
normalized_arr = (i3 - min_val) / (max_val - min_val)
print(normalized_arr)

#min2=0.2
#max2=0.4
mean = np.mean(normalized_arr)
std_dev = np.std(normalized_arr)

# Scale the array to have mean 0.5 and standard deviation 0 and 1
normalized_arr = (normalized_arr - mean) / std_dev * 0.5 + 0.5
normalized_arr = np.clip(normalized_arr, 0, 1)
print(normalized_arr)

#print(normalized_arr)
with open("output/tcf.txt", "w") as file:
    for aa in range(488):
        file.write("set_color me"+str(aa)+", ["+str(normalized_arr[aa])+","+str(1-normalized_arr[aa])+",0]\n")

with open("output/temp.txt", "w") as file:
    for aa in range(488):
        file.write("cmd.color(\'me"+str(aa)+"\', \'resi "+str(aa+1)+"\')\n")

exit()
if False:
    if plotwhat==3:
        X=bCV
        gmm=a_gmm
        plt.hist(X, bins=30, density=True, alpha=0.5, color='skyblue', label='Histogram')

        # Generate a range of values for X
        x_vals = np.linspace(-0, 50, 1000) #b CV
        #x_vals = np.linspace(35, 45, 1000) #a CV

        # Plot the Gaussian Mixture Model
        for i in range(2):
            # Extract the parameters of the i-th Gaussian component
            mean = gmm.means_[i][0]
            std_dev = np.sqrt(gmm.covariances_[i][0][0])
            weight = gmm.weights_[i]
            
            # Evaluate the PDF of the i-th Gaussian component
            pdf_vals = weight * np.exp(-0.5 * ((x_vals - mean) / std_dev) ** 2) / (std_dev * np.sqrt(2 * np.pi))
            
            # Plot the component PDF
            plt.plot(x_vals, pdf_vals, label=f'Component {i+1}')

        plt.show()

exit()
if False:
    a_p=a_gmm.predict_proba(aCV.reshape(-1, 1))
    b_p=b_gmm.predict_proba(aCV.reshape(-1, 1))
    #anb_p=anb_gmm.predict_proba(np.column_stack((aCV,bCV)))
    anb_p=np.zeros(shape=(snapEnd-snapStart,num_components,num_components))
    for s in range(snapEnd-snapStart):
        anb_p[s]=np.outer(a_p[s],b_p[s])

    numAminos=len(nrgs)
    AA_DDAb=np.zeros(shape=(numAminos,num_components,num_components))
    S_AA_DDAb=np.zeros(shape=(numAminos))
    for aa in range(numAminos):
        DDAb=np.zeros(shape=(num_components,num_components))
        for n1 in range(num_components):
            for n2 in range(num_components):
                f_anb=np.sum(nrgs[aa]*anb_p.T[n1][n2])
                f_a=np.sum(nrgs[aa]*a_p.T[n1])
                f_b=np.sum(nrgs[aa]*b_p.T[n2])
                f=np.sum(nrgs[aa])
                DDAb[n1][n2]=(f_anb/(f_a))*(f/(f_b))
                #if(not DDAb[n1][n2]):
                #   print(f_anb,f_a,f_b,f,DDAb[n1][n2])
        AA_DDAb[aa]=-np.log(DDAb)
        thing=AA_DDAb[aa]*anb_p
        S_AA_DDAb[aa]=np.average(AA_DDAb[aa]*anb_p/np.average(anb_p,axis=0))

        
# Min-max normalization
S_AA_DDAb = np.exp(S_AA_DDAb)
min_val = S_AA_DDAb.min()
max_val = S_AA_DDAb.max()
normalized_arr = (S_AA_DDAb - min_val) / (max_val - min_val)

numAminos=488*5

with open("check.txt", "w") as file:
    for aa in range(numAminos):
        file.write(str(aa/5)+','+str(normalized_arr[aa])+"\n")
exit()
with open("tcf.txt", "w") as file:
    for aa in range(numAminos):
        file.write("set_color me"+str(aa)+", ["+str(normalized_arr[aa])+","+str(1-normalized_arr[aa])+",0]\n")

with open("temp.txt", "w") as file:
    for aa in range(numAminos):
        file.write("cmd.color(\'me"+str(aa)+"\', \'resi "+str(aa+1)+"\')\n")
exit()
if pdfmode=='bins':
    # Set the number of bins
    num_bins = 20
    aHist, aBins = np.histogram(aCV, bins=num_bins)
    bHist, bBins = np.histogram(bCV, bins=num_bins)
    abHist, a_edges, b_edges = np.histogram2d(aCV, bCV, bins=(num_bins, num_bins))
    Ifxy=np.zeros(shape=(numAminos,num_bins,num_bins))
    Ifx=np.zeros(shape=(numAminos,num_bins))
    Ify=np.zeros(shape=(numAminos,num_bins))
    If=np.zeros(shape=(numAminos))
    tcfArg=np.zeros(shape=(numAminos,num_bins,num_bins))
    for s in range(len(aCV)):
        if s%1000==0:
            print(str(100*s/len(aCV))+"%")
        aaE=nrgs[:,s]#np.load("snapshots/e"+str(s)+".npy")
        aVal=aCV[s]
        bVal=bCV[s]
        aBin=np.digitize(aVal,aBins,right=True)-1
        bBin=np.digitize(bVal,bBins,right=True)-1
        for amino in range(numAminos):
            ebu=aaE[amino]
#            print(-aaE[amino])
#            print(ebu)
            Ifxy[amino][aBin][bBin]+=ebu/abHist[aBin][bBin] #Getting divide by zero here...
            Ifx[amino][aBin]+=ebu/aHist[aBin]
            Ify[amino][bBin]+=ebu/bHist[bBin]
            If[amino]+=ebu/(snapEnd-snapStart)
    i3=np.zeros(shape=(numAminos))
    for am in range(numAminos):
        for x in range(num_bins):
            for y in range(num_bins):
                if(Ifx[am][x]*Ify[am][y]==0):
                    tcfArg[am][x][y]=1
                tcfArg[am][x][y]=If[am]*Ifxy[am][x][y]/(Ifx[am][x]*Ify[am][y])
                i3[am]+=np.log(tcfArg[am][x][y])*abHist[x][y]/(snapEnd-snapStart) #filled with nan's
    print(tcfArg)
    np.save("output/tcfArg.npy",tcfArg)
    print(i3)
    np.save("output/i3.npy",i3)
    exit()

exit()
if pdfmode=='gmm':
    num_components=2 ### Need to update other numbers if this is changed
    # Fit a Gaussian Mixture Model with two components
    b_gmm = GaussianMixture(n_components=2).fit(bCV.reshape(-1, 1))
    a_gmm = GaussianMixture(n_components=2).fit(aCV.reshape(-1, 1))

    if plotwhat==3:
        X=bCV
        gmm=b_gmm
        plt.hist(X, bins=30, density=True, alpha=0.5, color='skyblue', label='Histogram')

        # Generate a range of values for X
        x_vals = np.linspace(-400, 400, 1000) #b CV
        #x_vals = np.linspace(35, 45, 1000) #a CV

        # Plot the Gaussian Mixture Model
        for i in range(2):
            # Extract the parameters of the i-th Gaussian component
            mean = gmm.means_[i][0]
            std_dev = np.sqrt(gmm.covariances_[i][0][0])
            weight = gmm.weights_[i]
            
            # Evaluate the PDF of the i-th Gaussian component
            pdf_vals = weight * np.exp(-0.5 * ((x_vals - mean) / std_dev) ** 2) / (std_dev * np.sqrt(2 * np.pi))
            
            # Plot the component PDF
            plt.plot(x_vals, pdf_vals, label=f'Component {i+1}')

        plt.show()



    a_p=a_gmm.predict_proba(aCV.reshape(-1, 1))
    b_p=b_gmm.predict_proba(aCV.reshape(-1, 1))
    #anb_p=anb_gmm.predict_proba(np.column_stack((aCV,bCV)))
    anb_p=np.zeros(shape=(snapEnd-snapStart,num_components,num_components))
    for s in range(snapEnd-snapStart):
        anb_p[s]=np.outer(a_p[s],b_p[s])

    numAminos=len(nrgs)
    AA_DDAb=np.zeros(shape=(numAminos,num_components,num_components))
    S_AA_DDAb=np.zeros(shape=(numAminos))
    for aa in range(numAminos):
        DDAb=np.zeros(shape=(num_components,num_components))
        for n1 in range(num_components):
            for n2 in range(num_components):
                f_anb=np.sum(nrgs[aa]*anb_p.T[n1][n2])
                f_a=np.sum(nrgs[aa]*a_p.T[n1])
                f_b=np.sum(nrgs[aa]*b_p.T[n2])
                f=np.sum(nrgs[aa])
                DDAb[n1][n2]=(f_anb/(f_a))*(f/(f_b))
                #if(not DDAb[n1][n2]):
                #   print(f_anb,f_a,f_b,f,DDAb[n1][n2])
        AA_DDAb[aa]=-np.log(DDAb)
        thing=AA_DDAb[aa]*anb_p
        S_AA_DDAb[aa]=np.average(AA_DDAb[aa]*anb_p/np.average(anb_p,axis=0))

        
# Min-max normalization
S_AA_DDAb = np.exp(S_AA_DDAb)
min_val = S_AA_DDAb.min()
max_val = S_AA_DDAb.max()
normalized_arr = (S_AA_DDAb - min_val) / (max_val - min_val)

numAminos=488*5

with open("check.txt", "w") as file:
    for aa in range(numAminos):
        file.write(str(aa/5)+','+str(normalized_arr[aa])+"\n")
exit()
with open("tcf.txt", "w") as file:
    for aa in range(numAminos):
        file.write("set_color me"+str(aa)+", ["+str(normalized_arr[aa])+","+str(1-normalized_arr[aa])+",0]\n")

with open("temp.txt", "w") as file:
    for aa in range(numAminos):
        file.write("cmd.color(\'me"+str(aa)+"\', \'resi "+str(aa+1)+"\')\n")






exit()

exit()

# show heatmap

heatmap, xedges, yedges = np.histogram2d(bCV, aCV, bins=20)

xmesh, ymesh = np.meshgrid(xedges[:-1], yedges[:-1])

plt.pcolormesh(xmesh, ymesh, heatmap.T, cmap='viridis')
plt.colorbar(label='Density')

plt.xlabel('# of Ligands Bound')
plt.ylabel('Channel Energy Barrier')
plt.title('Energy Barrier vs # of Ligands Bound Heat Density Map')

plt.show()

exit()
import numpy as np
from numpy import array, linspace
import sys
from scipy.signal import argrelextrema
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
import math
from scipy.stats import gaussian_kde
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from sklearn.mixture import GaussianMixture


parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-p", "--plot", default=0,type=int, help="0 for no plots, 1 for individual plots, 2 for combined plots")
parser.add_argument("-f", "--pdf", default='gmm',type=str, help="probability density function estimation method")
parser.add_argument("-s", "--start", default=400,type=int, help="starting snapshot")
parser.add_argument("-e", "--end", default=4400,type=int, help="ending snapshot")
parser.add_argument("-a", "--activation", default='namdBarrier/barrier.npy',type=str, help="activation CV")
parser.add_argument("-b", "--binding", default='ligCV/bindingPCA1.npy',type=str, help="binding CV")
parser.add_argument("-c", "--bias", default='nrgs/nrgs.npy',type=str, help="amino acid interaction energies")
args = vars(parser.parse_args())
plotwhat=args["plot"]
pdfmode=args["pdf"]
b=args["binding"]
a=args["activation"]
c=args["bias"]
snapStart=args["start"]
snapEnd=args["end"]

bCV=np.load(b)[0:snapEnd-snapStart]
aCV_E=np.load(a)[0:snapEnd-snapStart]
#aCV = np.exp(-aCV_E/(0.00198721319*310)) # Using just energy for now...
aCV=aCV_E
nrgs_E=np.load(c).T[0:snapEnd-snapStart].T
nrgs_E=(nrgs_E.T-nrgs_E.T.mean(axis=0)).T
nrgs = np.exp(-nrgs_E/(0.00198721319*310))


###THIS PROGRAM CALCULATES THE INTEGRATED TCF FOR EACH AMINO ACID ACCORDING TO THE INPUT PROJECTIONS

if plotwhat==2:
    # Create a 2D histogram to represent the density
    heatmap, xedges, yedges = np.histogram2d(bCV, aCV, bins=400)

    # Create a meshgrid for plotting
    xmesh, ymesh = np.meshgrid(xedges[:-1], yedges[:-1])

    # Plot the heat density map
    plt.pcolormesh(xmesh, ymesh, heatmap.T, cmap='viridis')
    plt.colorbar(label='Density')

    # Add labels and title
    plt.xlabel(b)
    plt.ylabel(a)
    plt.title('Heat Density Map')

    # Show the plot
    plt.show()

elif plotwhat==1:
    #I need to plot the two CVs, to see what kind of kernel I must pick.
    # Create subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Plot histogram for the first array
    ax1.hist(bCV, bins=400, edgecolor='black')
    ax1.set_title(b)
    ax1.set_xlabel('Values')
    ax1.set_ylabel('Frequency')

    # Plot histogram for the second array
    ax2.hist(aCV, bins=400, edgecolor='black')
    ax2.set_title(a)
    ax2.set_xlabel('Values')
    ax2.set_ylabel('Frequency')

    plt.tight_layout()

    # Show the plot
    plt.show()




#print(bCV.shape)
#print(bCV.reshape(-1,1).shape)
#exit()
if pdfmode=='gmm':
    num_components=2 ### Need to update other numbers if this is changed
    # Fit a Gaussian Mixture Model with two components
    b_gmm = GaussianMixture(n_components=2).fit(bCV.reshape(-1, 1))
    #print(b_gmm.means_)
    #print(b_gmm.covariances_)
    #print(b_gmm.weights_)
    #print('-')
    a_gmm = GaussianMixture(n_components=2).fit(aCV.reshape(-1, 1))
    #print(a_gmm.means_)
    #print(a_gmm.covariances_)
    #print(a_gmm.weights_)
    #print('-')
    # anb_gmm = GaussianMixture(n_components=2).fit(np.column_stack((aCV,bCV)))
    #print(anb_gmm.means_)
    #print(anb_gmm.covariances_)
    #print(anb_gmm.weights_)
    #print('-')

    if plotwhat==3:
        X=bCV
        gmm=b_gmm
        plt.hist(X, bins=30, density=True, alpha=0.5, color='skyblue', label='Histogram')

        # Generate a range of values for X
        x_vals = np.linspace(-400, 400, 1000) #b CV
        #x_vals = np.linspace(35, 45, 1000) #a CV

        # Plot the Gaussian Mixture Model
        for i in range(2):
            # Extract the parameters of the i-th Gaussian component
            mean = gmm.means_[i][0]
            std_dev = np.sqrt(gmm.covariances_[i][0][0])
            weight = gmm.weights_[i]
            
            # Evaluate the PDF of the i-th Gaussian component
            pdf_vals = weight * np.exp(-0.5 * ((x_vals - mean) / std_dev) ** 2) / (std_dev * np.sqrt(2 * np.pi))
            
            # Plot the component PDF
            plt.plot(x_vals, pdf_vals, label=f'Component {i+1}')

        plt.show()


    # if plotwhat==4:
    #     X=np.column_stack((aCV,bCV))
    #     gmm=anb_gmm
    #     # Create a scatter plot of the data points
    #     plt.scatter(X[:, 0], X[:, 1], alpha=0.5, label='Data Points')

    #     # Generate a range of values for X and Y to plot the GMM components
    #     x_vals = np.linspace(35, 45, 100)
    #     y_vals = np.linspace(-400, 400, 100)
    #     X_grid, Y_grid = np.meshgrid(x_vals, y_vals)
    #     XY = np.column_stack([X_grid.ravel(), Y_grid.ravel()])
    #     Z = -gmm.score_samples(XY)
    #     Z = Z.reshape(X_grid.shape)

    #     # Plot contour plot of the GMM
    #     plt.contour(X_grid, Y_grid, Z, levels=np.logspace(0, 3, 10), cmap='viridis', alpha=0.5)

    #     # Plot each component of the GMM as a separate cluster
    #     for i in range(2):
    #         mean = gmm.means_[i]
    #         cov = gmm.covariances_[i]
    #         plt.scatter(mean[0], mean[1], s=100, marker='x', color='red', label=f'Cluster {i+1} Mean')
    #     plt.show()

    a_p=a_gmm.predict_proba(aCV.reshape(-1, 1))
    b_p=b_gmm.predict_proba(aCV.reshape(-1, 1))
    #anb_p=anb_gmm.predict_proba(np.column_stack((aCV,bCV)))
    anb_p=np.zeros(shape=(snapEnd-snapStart,num_components,num_components))
    for s in range(snapEnd-snapStart):
        anb_p[s]=np.outer(a_p[s],b_p[s])

    numAminos=len(nrgs)
    AA_DDAb=np.zeros(shape=(numAminos,num_components,num_components))
    S_AA_DDAb=np.zeros(shape=(numAminos))
    for aa in range(numAminos):
        DDAb=np.zeros(shape=(num_components,num_components))
        for n1 in range(num_components):
            for n2 in range(num_components):
                f_anb=np.sum(nrgs[aa]*anb_p.T[n1][n2])
                f_a=np.sum(nrgs[aa]*a_p.T[n1])
                f_b=np.sum(nrgs[aa]*b_p.T[n2])
                f=np.sum(nrgs[aa])
                DDAb[n1][n2]=(f_anb/(f_a))*(f/(f_b))
                #if(not DDAb[n1][n2]):
                #   print(f_anb,f_a,f_b,f,DDAb[n1][n2])
        AA_DDAb[aa]=-np.log(DDAb)
        thing=AA_DDAb[aa]*anb_p
        S_AA_DDAb[aa]=np.average(AA_DDAb[aa]*anb_p/np.average(anb_p,axis=0))

#print(np.average(anb_p,axis=0).shape)
#exit()
#print(S_AA_DDAb)

#top10=np.argpartition(S_AA_DDAb, -10)[-10:]
#print(top10)
#print(S_AA_DDAb[top10])
        
# Min-max normalization
S_AA_DDAb = np.exp(S_AA_DDAb)
min_val = S_AA_DDAb.min()
max_val = S_AA_DDAb.max()
normalized_arr = (S_AA_DDAb - min_val) / (max_val - min_val)

with open("check.txt", "w") as file:
    for aa in range(numAminos):
        file.write(str(aa)+','+str(normalized_arr[aa])+"\n")




exit()
with open("tcf.txt", "w") as file:
    for aa in range(numAminos):
        file.write("set_color me"+str(aa)+", ["+str(normalized_arr[aa])+","+str(1-normalized_arr[aa])+",0]\n")

with open("temp.txt", "w") as file:
    for aa in range(numAminos):
#        file.write("set_color "+str(aa)+", [0,0,"+str(normalized_arr[aa])+"]\n")
        file.write("cmd.color(\'me"+str(aa)+"\', \'resi "+str(aa+1)+"\')\n")



exit()



exit()
aaE=np.load("snapshots/e1.npy")
numAminos=len(aaE)-1
# Begin TCF calculation
if pdfmode=='bins':
    # Set the number of bins
    num_bins = 20
    aHist, aBins = np.histogram(aCV, bins=num_bins)
    bHist, bBins = np.histogram(bCV, bins=num_bins)
    abHist, a_edges, b_edges = np.histogram2d(aCV, bCV, bins=(num_bins, num_bins))
    Ifxy=np.zeros(shape=(numAminos,num_bins,num_bins))
    Ifx=np.zeros(shape=(numAminos,num_bins))
    Ify=np.zeros(shape=(numAminos,num_bins))
    If=np.zeros(shape=(numAminos))
    tcfArg=np.zeros(shape=(numAminos,num_bins,num_bins))
    for s in range(snapStart,snapEnd):
        if s%1000==0:
            print(str(100*(s-snapStart)/(snapEnd-snapStart))+"%")
        aaE=np.load("snapshots/e"+str(s)+".npy")
        aVal=aCV[s-snapStart]
        bVal=bCV[s-snapStart]
        aBin=np.digitize(aVal,aBins,right=True)-1
        bBin=np.digitize(bVal,bBins,right=True)-1
        for amino in range(numAminos):
            ebu=np.exp(-aaE[amino])
#            print(-aaE[amino])
#            print(ebu)
            Ifxy[amino][aBin][bBin]+=ebu/abHist[aBin][bBin] #Getting divide by zero here...
            Ifx[amino][aBin]+=ebu/aHist[aBin]
            Ify[amino][bBin]+=ebu/bHist[bBin]
            If[amino]+=ebu/(snapEnd-snapStart)
    i3=np.zeros(shape=(numAminos))
    for am in range(numAminos):
        for x in range(num_bins):
            for y in range(num_bins):
                if(Ifx[am][x]*Ify[am][y]==0):
                    tcfArg[am][x][y]=1
                tcfArg[am][x][y]=If[am]*Ifxy[am][x][y]/(Ifx[am][x]*Ify[am][y])
                i3[am]+=np.log(tcfArg[am][x][y])*abHist[x][y]/(snapEnd-snapStart) #filled with nan's
    print(tcfArg)
    np.save("tcfArg.npy",tcfArg)
    print(i3)
    np.save("i3.npy",i3)
    exit()

elif pdfmode=='kde':
    
    exit()

elif pdfmode=='clusters':
    exit()

exit()
# Fit a kernel to the CVs to obtain a probability density function
combined=np.column_stack((cv1,cv2)).T
kde_combined = gaussian_kde(combined)
kde_1 = gaussian_kde(cv1)
kde_2 = gaussian_kde(cv2)

exit()