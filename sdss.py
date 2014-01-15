#!/Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python


import pyfits as py
import matplotlib.pyplot as plt
import pylab
import os, sys
import numpy as np


plt.rcParams["text.usetex"] = "True"

POSITION_DELTA = 0.01
VARIANCE_DELTA = 1


def read_sdss(name):
    """    
    Read SDSS spectra fits. Plotting capability.

    :INPUTS:
       name: Filename of fits spectra
       
    :OUTPUTS:
       wave: wavelength array - Angstroms

       flux: flux array  - Flux lambda erg/s/cm^2/Angstrom
    """
    flux=py.getdata(name,0)
    wdel=py.getval(name,'CD1_1',0)
    w0=py.getval(name,'CRVAL1',0)
    wave= 10.0**(w0+wdel*np.arange(len(flux[0])))
    
    return(wave,flux[0]*1e-17)




######################################################################
#
#		PLOTTING FUNCTIONS
#
#
######################################################################


def plot(name, lmin, lmax, save=True , norm = True):



	fig = plt.figure()

	ax = fig.add_subplot(111)

	wave, flux = read_sdss(name)

	if norm:
		flux = normalise(wave, flux, lmin)

	ax.plot(wave, flux)

	plt.xlim(lmin, lmax)

	savename = "%s_%i_%i.png" %( name[:-3], lmin, lmax)

	ax.set_xlabel("$\lambda (\AA)$")
	ax.set_ylabel("Flux")

	if save: 	
		plt.savefig( savename )

	#os.system("open -a preview %s" % savename)

def plotall(lmin, lmax):

	fig = plt.figure()

	ax = fig.add_subplot(111)

	f = open("lsfile", "r")

	for line in f:
		filename = line.split()[0]

		wave, flux = read_sdss(filename)

		flux = normalise(wave, flux, lmin)

		ax.plot(wave, flux)
		plt.xlim(lmin, lmax)

	ax.set_xlabel("$\lambda (\AA)$")
	ax.set_ylabel("Flux")

	savename = "all_%i_%i.png" % (lmin, lmax)

	plt.savefig( savename )


#def all_files():

def normalise(wave, flux, l):

	i = 0

	while wave[i] <= l:

		ipt = i

		i+=1

	comp = np.array(flux)

	delta = 10
	n = len(comp[ipt-delta:ipt+delta])


	norm_factor = np.sum(comp[ipt-delta:ipt+delta])	 / n

	flux = flux / norm_factor

	return flux



def get_variance(flux):
	'''get variance of normalised flux'''

	n = len(flux)
	mean = np.sum(flux) / n

	sum_1 = 0

	for i in range(n):

		f = flux[i]

		df = abs(f - mean)

		sum_1 += df

	return sum_1 / n



######################################################################
#
#		FUNCTIONS FOR CREATING INCLINATION DEPENDENT PLOTS
#
#
######################################################################


def get_all_incl(filename = "lsfile"):

	all_files = np.loadtxt(filename, dtype="string")

	print "Total files: %i" % len(all_files)

	incs = []
	names = []
	types = []

	for i in range(len(all_files)):

		inclin, cv_type = get_incl(all_files[i])

		if inclin < 90.0 and inclin > 0.0:
			incs.append(inclin)
			names.append(all_files[i])
			types.append(cv_type)

	print "Files with inclinations: %i" % len(incs)

	incs = np.array(incs)
	names = np.array(names)

	return incs, names, types



def get_incl(name, ritter_data = "ritter_v7.19.fits"):

	from math import fabs

	ra = py.getval(name,'RAOBJ',0)
	de = py.getval(name,'DECOBJ',0)

	#print "Finding", ra, de

	array = py.getdata(ritter_data, 0)

	ra2000 = array["_RAJ2000"]
	de2000 = array["_DEJ2000"]

	bool_ra = ( np.fabs(ra2000 -ra) < POSITION_DELTA )
	bool_dec = ( np.fabs(de2000 - de) < POSITION_DELTA)

	#print ra2000

	#print np.sum(bool_ra*bool_dec)

	Found = False

	#Found = (np.sum(bool_ra*bool_dec) != 1)

	for i in range(len(ra2000)):
		if fabs(ra2000[i] - ra) < POSITION_DELTA:
			if fabs(de2000[i] - de) < POSITION_DELTA:
				Found = True
				ipt = i

	if Found:

		inc = array["Incl"]

		types = array["Type1"]

		inclin = inc[ipt]

		cv_type = types[ipt]

		#if inclin == "":
		#	print "Found match,", name, ra, de, "None"
		if inclin < 90.0 and inclin > 0.0:
			print "Found match, and inclination!", name, ra, de, inclin 
		#else:
		#	print "Found match,", name, ra, de, "--" 

		return inclin, cv_type


	else:
		print "no RA/dec match", name, ra, de
		return 0


def plot_inclinations(incs, names, types, lmin, lmax, fname = "lsfile", norm = False, save = "incs.png", ymax = 4, obj = "NL"):

	#incs, names = get_all_incl(filename = fname)

	fig = plt.figure(figsize=(8.3,11.7),dpi=80)


	ax1 = fig.add_subplot(411)
	ax2 = fig.add_subplot(412)
	ax3 = fig.add_subplot(413)
	ax4 = fig.add_subplot(414)

	ax1.set_title("80+")
	ax2.set_title("70+")
	ax3.set_title("60+")
	ax4.set_title("60-")


	for i in range(len(names)):

		wave, flux = read_sdss(names[i])

		if norm:
			flux = normalise(wave, flux, lmin)

		lab = str(incs[i]) + " " + str(types[i])

		print types[i]

		if types[i] == obj:

			if incs[i] > 80.0:
				ax1.plot(wave, flux, label = lab)

			elif incs[i] > 70.0:
				ax2.plot(wave, flux, label = lab)

			elif incs[i] > 60.0:
				ax3.plot(wave, flux, label = lab)

			else:
				ax4.plot(wave, flux, label = lab)

		#plt.xlim(lmin, lmax)

	#plt.xlim(lmin, lmax)

	ax4.set_xlabel("$\lambda (\AA)$")

	ax1.set_ylabel("Flux")
	ax2.set_ylabel("Flux")
	ax3.set_ylabel("Flux")
	ax4.set_ylabel("Flux")

	ax1.legend()
	ax2.legend()
	ax3.legend()
	ax4.legend()

	ax1.set_xlim(lmin, lmax)
	ax2.set_xlim(lmin, lmax)
	ax3.set_xlim(lmin, lmax)
	ax4.set_xlim(lmin, lmax)

	ax1.set_ylim(0.5, ymax)
	ax2.set_ylim(0.5, ymax)
	ax3.set_ylim(0.5, ymax)
	ax4.set_ylim(0.5, ymax)

	plt.savefig(save)


	import os

	os.system("open -a preview %s" % save)

	return 0


import urllib2







######################################################################
#
#		FUNCTIONS FOR PLOTS BY TYPE
#
#
######################################################################


def get_type(name, ritter_data = "ritter_v7.20.fits"):
	from math import fabs

	ra = py.getval(name,'RAOBJ',0)
	de = py.getval(name,'DECOBJ',0)

	array = py.getdata(ritter_data, 0)

	ra2000 = array["_RAJ2000"]
	de2000 = array["_DEJ2000"]

	bool_ra = ( np.fabs(ra2000 -ra) < 0.1 )
	bool_dec = ( np.fabs(de2000 - de) < 0.1)

	Found = False

	#Found = (np.sum(bool_ra*bool_dec) != 1)

	for i in range(len(ra2000)):
		if fabs(ra2000[i] - ra) < 0.01:
			if fabs(de2000[i] - de) < 0.01:
				Found = True
				ipt = i
	if Found:

		cv_type = array["Type1"][ipt]
		type2 = array["Type2"][ipt]
		obj_name = array["Name"][ipt]


		ritter_class = Ritter ( name, obj_name, ra, de, cv_type, type2 )

		return ritter_class


	else:
		print "no RA/dec match", name, ra, de
		return 0


def get_all_types(filename = "lsfile"):	

	all_files = np.loadtxt(filename, dtype="string")

	print "Total files: %i" % len(all_files)

	all_objects = []

	for i in range(len(all_files)):

		CV_object = get_type(all_files[i])

		all_objects.append(CV_object)

	# return an array of ritter class instances
	return np.array(all_objects)





def collect_stats(ritter):

	types = np.array(["NL", "DN", "Na", "Nr", "Nb", "N", "CV", "SS"])

	types2 = np.array([ np.array(["UX", "SW", "IP", "AM", "ZC", "NS", "VY", "AC","",]),
	               np.array(["SU", "WZ", "UG", "ZC"]) ])


	count1 = np.zeros(len(types))
	count2 = np.array( [ np.zeros(len(types2[0])), np.zeros(len(types2[1])) ] )


	for i in range(len(ritter)):

		CVtype = ritter[i].type

		CVtype2 = ritter[i].type2

		bool_type = ( types == CVtype)

		bool_type2_nl = ( types2[0] == CVtype2)

		bool_type2_dn = ( types2[1] == CVtype2)
		#bool_type2_dn = ( types == CVtype2[1])


		print bool_type, CVtype
		print bool_type2_nl, CVtype2

		count1 += bool_type
		count2[0] += bool_type2_nl * bool_type[0]
		count2[1] += bool_type2_dn * bool_type[1]

	for i in range(len(types)):
		print types[i], int(count1[i])

	for i in range(len(types2)):

		for j in range(len(types2[i])):
			print types[i], types2[i][j], int(count2[i][j])

	return count1, count2, types, types2


def plot_stats(pop, titles = ["Ritter", "Ritter matched w/ SDSS"] ):

	fig = plt.figure(figsize=(8.3,11.7),dpi=80)

	col = ["r", "b"]


	for iplot in range(len(pop)):

		count1, count2, t1, t2 = collect_stats(pop[iplot])
		width = 1.0 



		ax = fig.add_subplot (2*len(pop), 1, iplot+1)

		ticklocs = []
		ticklabels = []

		for i in range(len(count2)):

			imin = i*10
			imax = imin + len(count2[i])

			ind = np.arange(imin, imax) 


			ax.bar ( ind , count2[i], width, label = t1[i], color=col[i] )

			#label_tuple = 
			for j in range(len(ind)):
				ticklabels.append(t2[i][j])
				ticklocs.append( imin + j + (0.5*width) )

		ax.legend()
		ax.set_xticks( ticklocs)
		ax.set_xticklabels( ticklabels )
		ax.set_ylabel('Freq')
		ax.set_title(titles[iplot])

	#plt.savefig("subtypes_hist.png")
	#plt.clf()


	#fig = plt.figure()
	for iplot in range(len(pop)):

		count1, count2, t1, t2 = collect_stats(pop[iplot])
		width = 1.0 

		ax = fig.add_subplot (len(pop) * 2, 1, iplot+1 + 2)

		ticklocs = []
		ticklabels = []

		imin = 0
		imax = imin + len(count1)

		ind = np.arange(imin, imax) 


		ax.bar ( ind , count1, width, label = t1[i], color=col[i] )

		#label_tuple = 
		for j in range(len(ind)):
			ticklabels.append(t1[j])
			ticklocs.append( imin + j + (0.5*width) )

		#ax.legend()
		ax.set_xticks( ticklocs)
		ax.set_xticklabels( ticklabels )
		ax.set_ylabel('Freq')
		ax.set_title(titles[iplot])



	plt.savefig("types_hist.png")









class Ritter:
	'''
	Stores a number of variables from the Ritter CV catalog
	'''
	def __init__(self, _filename, _object_name, _ra, _dec, _type, _type2):
		self.filename = _filename
		self.object_name = _object_name
		self.ra = _ra 
		self.dec = _dec 
		self.type = _type 
		self.type2 = _type2



def plot_types(ritter_classes, lmin, lmax, obj="NL", fname = "lsfile", norm = False, ymin= 0.5, ymax = None, obj2="UX"):

	from math import fabs

	fig = plt.figure()


	ax1 = fig.add_subplot(111)


	for i in range(len(ritter_classes)):

		# check we want this type of object
		if ritter_classes[i].type == obj and ritter_classes[i].type2 == obj2:

			print "correct object name"

			wave, flux = read_sdss(ritter_classes[i].filename)	# read from sdss fits file


			if norm:
				flux = normalise(wave, flux, lmin)	# normalise flux

			sigma = get_variance(flux)	# get variance

			print sigma, VARIANCE_DELTA

			lab = ritter_classes[i].object_name

			if sigma < VARIANCE_DELTA:	# check not too noisy
				ax1.plot(wave, flux, label = lab)



	ax1.set_xlim(lmin, lmax)

	if ymax!= None:
		ax1.set_ylim(ymin, ymax)

	ax1.set_xlabel("$\lambda (\AA)$")
	ax1.axvline(x=6564, c="m")
	ax1.axvline(x=4861, c="m")
	ax1.axvline(x=4542, c="c")
	ax1.axvline(x=4686, c="c")
	ax1.axvline(x=4686, c="c")
	ax1.axvline(x=4471, c="y")
	ax1.axvline(x=4922, c="y")
	ax1.axvline(x=5876, c="y")
	ax1.axvline(x=6678, c="y")

	ax1.set_ylabel("Flux")



	ax1.legend()


	save = "type_%s_%s_%i_%i.png" % (obj, obj2, lmin, lmax)

	plt.savefig(save)

	#plt.clf()

	import os

	os.system("open -a preview %s" % save)

	return 0










######################################################################
#
#		RETRIEVING FUNCTIONS
#
#
######################################################################


def get_page(id):
	import urllib

	url = "http://cas.sdss.org/dr7/en/tools/explore/obj.asp?id=%s" % id

	urllib.urlretrieve (url, "id.html")

	return 0


def get_online_file_sdss(id):

	url = "http://cas.sdss.org/dr7/en/tools/explore/obj.asp?id=%s" % id

	get_online_file(url)

	return 0


def get_online_file(url):

	import urllib2

	file_name = url.split('/')[-1]

	u = urllib2.urlopen(url)
	f = open(file_name, 'wb')
	meta = u.info()
	file_size = int(meta.getheaders("Content-Length")[0])

	print "Downloading: %s Bytes: %s" % (file_name, file_size)

	file_size_dl = 0
	block_sz = 8192
	while True:
	    buffer = u.read(block_sz)
	    if not buffer:
	        break

	    file_size_dl += len(buffer)
	    f.write(buffer)
	    status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
	    status = status + chr(8)*(len(status)+1)
	    print status,

	f.close()

	return 0



def get_link(id):
	'''
	'''

	# pip install selenium, or get source from https://pypi.python.org/pypi/selenium
	from selenium import webdriver						
	from selenium.webdriver.common.keys import Keys


	# we'll open with firefox - can do this with other browsers
	browser = webdriver.Firefox()

	print "getting link for CV %s..." % id

	try: 
		browser.get('http://cas.sdss.org/dr7/en/tools/explore/obj.asp?id=%s' % id)


		# switch to the TOC fram on side of page
		browser.switch_to_frame("OETOC")

		xpath = '//*[@id="toc"]/table/tbody/tr[32]/td/a'



		allelements = browser.find_elements_by_xpath( xpath )

		if len(allelements) > 1:
			"Warning, multiple elements!"

		for i in allelements:
			new_url = i.get_attribute('href')

	except:
		browser.quit()
		print "error"
		return 1


	browser.quit()


	return new_url





def get_fits (url):

	# pip install selenium, or get source from https://pypi.python.org/pypi/selenium
	from selenium import webdriver						
	from selenium.webdriver.common.keys import Keys

	# we'll open with firefox - can do this with other browsers
	browser = webdriver.Firefox()

	print "getting fits from URL %s..." % url

	try: 
		browser.get(url)

		xpath = "/html/body/div[2]/table/tbody/tr[1]/td[2]/a"

		allelements = browser.find_elements_by_xpath( xpath )


		for i in allelements:

			fits_url = i.get_attribute('href')

	except:
		browser.quit()
		print "error"
		return 1

	if len(allelements) == 0:
		print "No fits for this CV, moving on."
		return 1


	browser.quit()

	return fits_url



def get_all_fits(filename = "catalog.dat"):


	all_files = np.loadtxt(filename, unpack = True, dtype="string")


	files = all_files[1]

	for i in range(len(files)):
		idname = files[i]

		url = get_link(idname)

		fits_url = get_fits (url)

		if fits_url != 1:
			get_online_file(fits_url)


	return 0




def signin(usr, pword, \
	url = "https://www.outlook.soton.ac.uk/owa/auth/\
	logon.aspx?replaceCurrent=1&url=https%3a%2f%2fwww\
	.outlook.soton.ac.uk%2fowa%2f"):

	'''uses selenium webdriver to'''

	from selenium import webdriver
	from selenium.webdriver.common.keys import Keys


	browser = webdriver.Firefox()

	browser.get(url)

	elem = browser.find_element_by_name('username')

	elem.send_keys(usr)

	elem = browser.find_element_by_name('password')

	elem.send_keys(pword + Keys.RETURN)

	return 0




#get_online_file(sys.argv[1])

#wave, flux = read_sdss(sys.argv[1])

#lmin = float(sys.argv[2])
#lmax = float(sys.argv[3])

