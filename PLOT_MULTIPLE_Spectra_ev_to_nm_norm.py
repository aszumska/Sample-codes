import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys, getopt
import csv
import os
import math
import matplotlib.colors as clr
from matplotlib.font_manager import FontProperties

def main(argv):
#default parametres
        file = ''
	sigma = 0.1	#broadening -s
        dL = 1		#mesh grid -m
        title='' 	#title of a plot -t
        ymax = 10	#zoom plot absorption max --ymax
        ymin = 0	#zoom plot absorption min --ymin
        Lmax = 800	#zoom plot wavelength max --lmax
        Lmin = 100	#zoom plot wavelength min --lmin
        zoom = 0
#        Colors=["blue","orange","green", "yellow","black"]
        labels=0	#wavelengths labels for peaks on the plot -l=1
        hartree=0	#energy given in Hartree --hartree
        log=0		#Log scale of the plot y --log
	vlines=0
	steps=0
	legend ='upper right'
        try:
                opts, args = getopt.getopt(argv,"h:f:s:m:t:vlines:ymax:ymin:lmax:lmin:labels:hartree:log:intermediate:legend",["ymax=","ymin=","lmax=","lmin=","hartree=","labels=","log=","vlines=","intermediate=","legend="])
        except getopt.GetoptError:
                print 'test.py -f inputfile -s sigma -m mesh grid dE -t title of the plot default(file name) --ymax --ymin zoom plot, --lmax --lmin wavelength range  --labels labels with wavelengths (1 on, 0 off(default) --hartree energy given in hartree 0 off 1 on --log log scale on y (1 on 0 off) --vlines 0 off 1 on --intermediate giva a number of steps to include between plots --legend position of legend (default center left), can be 1-4'
                sys.exit(2)
        for opt, arg in opts:
                if opt == '-h':
                        print 'test.py -f inputfile -s sigma -m mesh grid dE -t title of the plot default(file name) --ymax --ymin zoom plot, --lmax --lmin wavelength range -labels labels with wavelengths (1 on, 0 off(default) --hartree energy given in hartree 0 off 1 on --log log scale on y (1 on 0 off) --vlines 1 on 0 off --intermediate give a number of steps to include between plots'
                        sys.exit()
		elif opt in ("-f"):
			file = arg
                elif opt in ("-s"):
                        sigma = float(arg)
                elif opt in ("-m"):
                        dL = float(arg)
                elif opt in ("-t"):
                        title = arg
                elif opt in ("--ymax"):
                        ymax = float(arg)
                        zoom = 1
                elif opt in ("--ymin"):
                        ymin = float(arg)
                        zoom = 1
                elif opt in ("--lmax"):
                        Lmax = float(arg)
                elif opt in ("--lmin"):
                        Lmin = float(arg)
                elif opt in ("--labels"):
                        if int(arg)==1:
                                labels = 1
                elif opt in ("--hartree"):
                        if int(arg)==1:
                                hartree = 1
                elif opt in ("--log"):
                        if int(arg)==1:
                                log = 1
		elif opt in ("--vlines"):
			if int(arg)==1:
				vlines =1
		elif opt in ("--intermediate"):
			intermediate = float(arg)
			steps =int(intermediate)
		elif opt in ("--legend"):
			legend = arg
        print 'Plot in the wavelngth range: ', Lmin, Lmax
        print 'Plot zoomed to absorbtion:', ymin, ymax

	filename = os.path.splitext(file)[0]
	outfile = filename + "sigma" + str(sigma)+".png"
	outtxt = filename + "sigma" + str(sigma)+ ".txt"
	print 'Input file: ', file 
	print 'Plot is saved to : ', outfile 
	print 'Txt data is saved to: ', outtxt
	print 'Step:', steps

	headers=pd.read_csv(file, header=None,nrows=1, sep='\s+')
	header= headers.values
	n_plots = len(header[0])
	data=pd.read_csv(file, delimiter ='\s+')
	data=data.values
	data2=np.transpose(data)
	data = np.nan_to_num(data2)
	k_plots=int(math.ceil(n_plots/2.))
	Ei=[float(0) for i in range(k_plots)]	#array with energies for peaks 
	fi=[float(0) for i in range(k_plots)]	#array with oscilator strengths for peaks
	Legend = ["" for i in range((steps+1)*(k_plots-1)+1)]
	for i in range(0,n_plots):		#sort columns to energies and oscilator strengths
		k=int(math.floor(i/2.))
		if i%2 == 0:
			Ei[k]=data[i]
			#print i, k, steps
			Legend[k*(steps+1)]=header[0][i]
			if steps!=0 and i!=(n_plots-2):
				for s in range(steps):
					Legend[k*(steps+1)+s+1]=str(k)+str(k+1)+'_'+str(s+1)
		else:
			fi[k]=data[i]
	print Legend
	#COLOR scheme
	Colors = [0 for x in range((k_plots-1)*(steps+1)+1)]
 	#cmap = plt.get_cmap('jet_r')
	grormute = ['#4DAF4A', '#377EB8', '#984EA3', '#E41A1C']
	cmap = clr.LinearSegmentedColormap.from_list('custom blue', grormute , N=256)
	for i, n in enumerate(np.linspace(1.0,0.0,((k_plots-1)*(steps+1)+1))):
		Colors[i]=cmap(float(n))
	
	#plot linespectrum
        fig = plt.figure()
        ax = fig.add_subplot(111)
	#minL=1239.8/max(Ei)
	print Ei
	print Ei.min()
	Enm=np.linspace(Lmin, Lmax, num=(Lmax-Lmin)/dL+1)
	fval=[float(0) for i in range((k_plots-1)*(steps+1)+1)]
	Einm=[float(0) for i in range(k_plots)]
        for k in range(0,int(k_plots)):
		k_steps=k*(steps+1)	#only main plots
	        for i in range(0,len(Ei[k])):
        	 	a=fi[k][i]/(2*sigma*math.sqrt(math.pi))
                	b=Ei[k][i]
	                c=2*sigma
        	        fval[k_steps] += a*np.exp(-(((1239.8/Enm)-b)**2)/(c**2))
                #plot linespectrum
		max=np.max(fval[k_steps])
		fval[k_steps] /= max
		Einm[k]=1239.8/Ei[k]
		if vlines==1:
		       plt.vlines(Einm[k],[0],fi[k])
	
	#Calculate intermediate steps 
        for k in range(0, int(k_plots)):
		fval_prev=fval[k*(steps+1)]
		if k!=(k_plots-1):
	        	fval_next=fval[(k+1)*(steps+1)]
			for s in range(steps):
				fval[k*(steps+1)+s+1]=float(steps+1-s-1)/float(steps+1)*fval_prev+float(s+1)/float(steps+1)*fval_next
	#Save txt file	
	for k in range((k_plots-1)*(steps+1)+1):
		print k 
                print Legend[k]
		uvvis = np.array([Enm,fval[k]])
                uvvis = uvvis.T
		if int(k) == 0:
			df=pd.DataFrame(Enm, columns = ['Wavelength'])
			df.to_csv(outtxt, index=False )
		df=pd.read_csv(outtxt)
		df[Legend[k]]=fval[k]
		df.to_csv(outtxt, index=False)
	#plot gaussian broaden
		LW=0
		if k%(steps+1)==0:	#linewidth
			LW=2		#for main lines
		else:
			LW=1		#for addon lines
	        if log!=1:
      			plt.plot(Enm,fval[k],label=Legend[k], color=Colors[k], lw=LW)
	        else:
        	        plt.semilogy(Enm,fval[k],label=Legend[k])
        fontP=FontProperties()
        fontP.set_size('small')
        lgd=plt.legend(loc=legend, borderaxespad=0., prop=fontP)#, bbox_to_anchor=(1,0.5))

	#adding nm labels 
	if int(labels) == 1:
		for k in range(0,int(k_plots)):
	        	for X, Y in zip(Einm[k], fi[k]):
	                	if X>=Lmin and X<=Lmax:
              	                	d=int((X-Lmin)/float(dL))
                                       # Annotate the points 
	     		                if X!=0:
                                               	Wave=int(X)
	                                        if fval[k][d]<ymax:
                                                        ax.annotate('{}'.format(Wave), xy=(Wave,fval[k][d]), xytext=(15,5), ha='right', textcoords='offset points',color=Colors[k])
               	                                else:
                       	                                ax.annotate('{}'.format(Wave), xy=(Wave,0.9*ymax), xytext=(15,5), ha='right', textcoords='offset points',color=Colors[k])
                                        Wave=0

	#formatting labelling etc
        plt.title(filename + title)
        plt.xlabel("Wavelength [nm]")
        plt.ylabel("Absorption")
        plt.xlim(Lmin,Lmax)
        if int(zoom) == 1:
                plt.ylim(ymin,ymax)
        #save the file as png or eps or anything else
        plt.savefig(outfile, bbox_extra_artists=(lgd,), bbox_inches='tight')
        #this makes the plot pop up in a new window after running the script
	plt.show()

if __name__ == "__main__":
        main(sys.argv[1:])
