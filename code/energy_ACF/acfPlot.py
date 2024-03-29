import numpy as np
import matplotlib.pyplot as plt

# Generate some sample data (you can replace this with your own array)
with open('output/observables/acf/acf.txt', 'r') as file:
    d2 = file.readlines()
    data=np.zeros(len(d2))
    for i in range(len(d2)):
        data[i] = float(d2[i])

# Compute autocorrelation function
def autocorr(x):
    result = np.correlate(x, x, mode='full')
    return result

autocorr_result = np.log(autocorr(data))



# Plot autocorrelation function
plt.figure(figsize=(10, 5))
plt.plot(autocorr_result)
plt.title('Autocorrelation Function')
plt.xlabel('Lag')
plt.ylabel('Autocorrelation')
plt.grid(True)
plt.show()

### 100 steps (x5000) equilibrated under the following settings

#MC IACCept 0 NSTEp 2e6 TEMP 310.15 -
#INBFrq 400 IECHeck 800 IMGFrq 800 IDOMcfrq 10 -
#IUNC 69 NSAVc 5000


#MOVE ADD MVTP TORS WEIGht 1 DMAX 180.0 FEWEr 1 -
#         SELE ALL END SELE ALL END

#MOVE ADD MVTP CROT WEIGht 10 DMAX 90.0 NLIMit 1 -
#         SELE ((TYPE N).OR.(TYPE CA).OR.(TYPE C)) END -
#         SELE (TYPE C) END  SELE (TYPE N) END -
#         SELE (RESNAME PRO .AND. TYPE CA) END -
#         SELE (RESNAME PRO .AND. TYPE  N) END - 
#         SELE (RESNAME CYS .AND. TYPE CA) END -
#         SELE (RESNAME CYS .AND. TYPE N) END

#MOVE ADD MVTP RTRN BYATom WEIGht 1
