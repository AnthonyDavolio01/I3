import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate

i3=np.load('output/I3.npy')

idd=[]
pseg=[]
for i in range(7,341):
    idd.append(i)
    pseg.append("PROA")
for i in range(10,348):
    idd.append(i)
    pseg.append("PROB")
for i in range(7,341):
    idd.append(i)
    pseg.append("PROC")
for i in range(10,348):
    idd.append(i)
    pseg.append("PROD")
for i in range(10,358):
    idd.append(i)
    pseg.append("PROE")


shi={}
for i in range(len(idd)):
    shi[str(idd[i])+pseg[i][3]]=i


plip={'6X3S':['65B','130B','67B','200A','157A'],
      '6X3T':['289A','228B','236B','265A','282A'],
      '6X3U':['100B','142E','205B','77E','210B'],
      '6X3V':['265A','285A','290A','289A','228B'],
      '6X3W':['228A','291B','294B','223A','305E','231A','280E'],
      '6X3X':['233B','262A','265A','265B','269B','285A','289A','58E','77E','100B','160B','203B','210B','205B','206B','102B','223A','227A','231A','284E','280E']
      }


rows=[]
apndd=[]
for a in plip.keys():
    for b  in plip[a]:
        if b[-1]=='A':
            s='C'
            if b not in apndd:
                apndd.append(b)
                rows.append(b+'/'+s)
        elif b[-1]=='B':
            s='D'
            if b not in apndd:
                apndd.append(b)
                rows.append(b+'/'+s)
        else:
            if b not in apndd:
                apndd.append(b)
                rows.append(b+'  ')
rows = sorted(rows, key=lambda x: (x[-1], int(x[0:-3])))


for a in plip.keys():
    for b  in plip[a]:
        if b[-1]=='A':
            s='C'
            plip[a].append(b[0:-1]+s)
        if b[-1]=='B':
            s='D'
            plip[a].append(b[0:-1]+s)

cmi={}
for a in plip.keys():
    cmi[a]=[]
    for b in plip[a]:
        cmi[a].append('%.3f'%(i3[shi[b]]))

r=[]
t=['res','6X3S','6X3T','6X3U','6X3V','6X3W','6X3X','CMI']

for row in rows:
    rr=[row]
    for a in plip.keys():
        if row[0:-2] in plip[a]:
            rr.append('X')
            cc=str(cmi[a][plip[a].index(row[0:-2])])
            if row[-2]=='/':
                cc=cc+' / '+str(cmi[a][plip[a].index(row[0:-3]+row[-1])])
        else:
            rr.append('-')
    rr.append(cc)
    r.append(rr)



hist, bins = np.histogram(i3, bins=30)
# Calculate bin centers
bin_centers = (bins[1:] + bins[:-1]) / 2

# Plot PDF
plt.bar(bin_centers, hist, width=bins[1]-bins[0], label='PDF', align='center', edgecolor='black')
plt.title('PDF of CMI')
plt.xlabel('Conditional Mutual Information')
plt.ylabel('# of amino acids')
plt.legend()
plt.grid(True)
plt.show()

print(tabulate(r, headers=t, tablefmt='grid'))
