#!/usr/bin/python  -tt
import sys
import commands
import numpy as np
from scipy.optimize import curve_fit
from pylab import *
'''
A Python code to Construct a modified Morse function fit for a hybrid semiconductor-metal-semiconductor MoS2 monolayer (2D material)
DOI: 10.1021/acs.jpcc.6b07917
'''
def intro(x):
    if x==1:
      print "Program to create a Morse Fitting along a reaction coordinate (Y axis)"
      print "'Warning: Enter all string inputs in quotes'"

########## Get Input ############################################
def extract_inputs():
'''
A function to process inputs
'''
    xyz = input("Input the xyz file [quotes] : ")
    ref = input("Input the reference atoms number (no quotes): " )
    coord = raw_input("Input the coordination number [seperated by tab space and no quotes] : ")
    Bond_length = raw_input("Input the "+str(len(coord.split()))+" coordinated bond length in the same order as above [seperated by tab space and no quotes] : ")
    coord_array=np.zeros((len(coord.split()),3))
    ref_array=np.zeros((1,3))
    xyz_data = np.loadtxt(xyz,skiprows=2,usecols=(1,2,3))
    ref_array=xyz_data[ref-1,:]
    i = 0; a = 0
    for i in coord.split():
        coord_array[a,:]=xyz_data[int(i)-1,:]
        a +=1
  
    dim = len(np.loadtxt('YvsE.dat',usecols=(0,)))
    return ref, coord.split(), ref_array, coord_array, Bond_length.split(), dim

########## Extracts distance for atomic reference ###################
def compute_distance(m,n):
'''
A function to convert distance along a reference axis to distance along the bond axis of three coordinated atoms
'''
    Y=np.loadtxt('YvsE.dat',usecols=(0,))
    E=np.loadtxt('YvsE.dat',usecols=(1,))
    unload_para=np.zeros((len(Y),3))
    iter=0;Y_val=0
    for Y_val in  Y:
        unload_para[iter,0]=np.sqrt((m[0]-n[0])**2+(Y_val-n[1])**2+(m[2]-n[2])**2)  ## Norm
        unload_para[iter,1]=float(E[iter])                                          ## Restore energy data
        print float(E[iter])-float(E[18])
        unload_para[iter,2]=float(E[iter])-float(E[18])                             ## Shift all data with reference to minimum energy!
        iter+=1
    return unload_para

############### Morse function of choice ###########################

def morse_func(x,p0,p1,p2,p3):
'''
The modified Morse function
'''
    return p0 + p1*((1 - np.exp(-p2*(x[:,0]-p3)))**2 + (1 - np.exp(-p2*(x[:,1]-p3)))**2 + (1 - np.exp(-p2*(x[:,2]-p3)))**2)


########################## Main ####################################

def main():
    (status, output) = commands.getstatusoutput('rm AvsE* fit* popt*')
    print output
    (status, output)= commands.getstatusoutput('ls -ltr')
    print output
    intro(int(1))
    (ref_num,coor_num,ref_val,rec_coord,Bond_length,Y_dim)=extract_inputs() #Received all inputs
    #Compute distances 
    k=0;p=0
    x_data = np.zeros((Y_dim,3))
    for k in coor_num:
        save_para=compute_distance(ref_val,rec_coord[p,:])
        np.savetxt('AvsE'+k+'.dat',save_para)
        x_data[:,p]=save_para[:,0]
        y_data=save_para[:,2]
        print "The bond length is :",Bond_length[p],"\n"
        p+=1
    bond_para = np.array(Bond_length).astype(np.float)
    print 'The first Bond length is ', bond_para[0]
    popt, pcov = curve_fit(morse_func, x_data, y_data,p0=(1.,1.,0.5,2.0)) #Non-linear curve fitting
    print popt
    l = 0 ;r=0
    for r in coor_num:
        print "The morse fit parameter is ", popt
        plot(x_data[:,l],morse_func(x_data,popt[0],popt[1],popt[2],popt[3]),'r-',lw=4.5,label='morse fit')
        plot(x_data[:,l],y_data,'gs',lw=4.5,label='dft')
        autoscale(tight=True)
        xlabel('Bond streching (Ang)',fontsize=18)
        ylabel('Energy (eV)',fontsize=17)
        ticklabel_format(useOffset=False)
        title('Morse fit with for Mo-'+r+' with reference to S-84',fontsize=17)
        grid(True)
        legend()
        savefig("fit"+r+".png")
        show()
        l+=1
    np.savetxt('popt_refined.dat',popt)

if __name__ == '__main__':
  main()
