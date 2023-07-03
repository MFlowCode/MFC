import numpy as np 
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, BSpline
from scipy.fft import fft, fftfreq


str1 = "prim.3.00."
str2 = ".dat"
str3 = "cons.12.00."
t_step = 0
plist_nopv = []
plist_nopv1 = []
plist_nopv2 = []
plist_nopv3 = []
plist_nopv_poly = []
tlist_nopv = []
dt = 0.0273 * 9.93440117e-7
tfinal = 94*40
delt = 94
index = 129
while(t_step < tfinal):
	if(t_step < 10):
		rstr = str1 + "00000" + str(t_step) + str2
		nrstr = str3 + "00000" + str(t_step) + str2
	elif(t_step < 100):
		rstr = str1 + "0000" + str(t_step) + str2
		nrstr = str3 + "0000" + str(t_step) + str2
	elif(t_step < 1000):
		rstr = str1  + "000" + str(t_step) + str2
		nrstr = str3  + "000" + str(t_step) + str2
	elif(t_step < 10000):
		rstr = str1  + "00" + str(t_step) + str2
		nrstr = str3  + "00" + str(t_step) + str2
	elif(t_step < 100000):
		rstr = str1  + "0" + str(t_step) + str2
		nrstr = str3  + "0" + str(t_step) + str2
	else:
		rstr = str1  + "" + str(t_step) + str2
		nrstr = str3  + "" + str(t_step) + str2

	pres = np.loadtxt("D/" + rstr)
	plist_nopv.append(pres[index][1])

	#plt.plot(pres[:, 0], pres[:, 1], label = r'$\sigma = 0.0$')
	
	plt.plot(pres[:, 0], pres[:, 1])
	plt.xlabel(r'$ x / R_0$')
	plt.ylabel(r'$ p / p_0$')
	#plt.ylim(0.95, 1.12)
	plt.grid(linestyle = '--', linewidth = 0.5)
	plt.savefig("pres_plots/" + "/sig_0.0/" +  rstr[:-4] + ".png", dpi = 200)	
	plt.pause(1)
	plt.show(block = False)
	plt.close()
	
	


	"""
	plt.plot(pres[:, 0], pres[:, 1], label = r'$\sigma = 0.01$')
	plt.xlabel(r'$ x / R_0$')
	plt.ylabel(r'$ p / p_0$')
	plt.ylim(0.95, 1.12)
	
	plt.legend()
	plt.grid(linestyle = '--', linewidth = 0.5)
	plt.savefig("pres_plots/" + "/all/" +  rstr[:-4] + ".png", dpi = 200)
	#plt.pause(1)
	#plt.show(block = False)
	plt.close()
	"""

	tlist_nopv.append(dt*t_step)
	t_step += delt

	#pres[:, 1] = (pres[:, 1] - 1 )/ (0.01)


	#print(pres[49, 1])


plist_nopv = np.array(plist_nopv)
plist_nopv1 = np.array(plist_nopv1)
plist_nopv2 = np.array(plist_nopv2)
plist_nopv3 = np.array(plist_nopv3)
tlist_nopv = np.array(tlist_nopv)

tlist_nopv = tlist_nopv 


xnew = np.linspace(tlist_nopv.min(), tlist_nopv.max(), 200) 
spl = make_interp_spline(tlist_nopv, plist_nopv, k=3)  # type: BSpline
power_smooth = spl(xnew)
plist_nopv = power_smooth.copy()

spl = make_interp_spline(tlist_nopv, plist_nopv1, k=3)  # type: BSpline
power_smooth = spl(xnew)
plist_nopv1 = power_smooth.copy()

spl = make_interp_spline(tlist_nopv, plist_nopv2, k=3)  # type: BSpline
power_smooth = spl(xnew)
plist_nopv2 = power_smooth.copy()

spl = make_interp_spline(tlist_nopv, plist_nopv3, k=3)  # type: BSpline
power_smooth = spl(xnew)
plist_nopv3 = power_smooth.copy()

tlist_nopv = xnew.copy()

"""
plt.plot(tlist_nopv * 10 ** 6, plist_nopv, label = r'$\sigma = 0.0$')
plt.plot(tlist_nopv * 10 ** 6, plist_nopv1, label = r'$\sigma = 0.1$')
plt.plot(tlist_nopv * 10 ** 6, plist_nopv2, label = r'$\sigma = 0.05$')
plt.plot(tlist_nopv * 10 ** 6, plist_nopv3, label = r'$\sigma = 0.01$')
plt.xlabel(r'$t$' + ' ' + r'$[\mu s]$')
plt.ylabel(r'$ p / p_0$')
plt.legend()
plt.grid(linestyle = '--', linewidth = 0.5)
plt.savefig("pres_r.png", dpi = 200)
plt.show()
plt.close()

delt = delt * 40 / 200
xf = fftfreq(200, d=delt*dt)[:200//2]
yf = fft(plist_nopv)
yf1 = fft(plist_nopv1)
yf2 = fft(plist_nopv2)
yf3 = fft(plist_nopv3)
xf = xf / 10 ** 6
plt.loglog(xf,  np.abs(yf[0:200//2])** 2, label = r'$\sigma = 0.0$')
plt.loglog(xf,  np.abs(yf2[0:200//2])** 2, label = r'$\sigma = 0.1$')
plt.loglog(xf,  np.abs(yf2[0:200//2])** 2, label = r'$\sigma = 0.05$')
plt.loglog(xf,  np.abs(yf3[0:200//2])** 2, label = r'$\sigma = 0.01$')
plt.xlabel(r'$f$' + ' ' + r'$[MHz]$')
plt.ylabel(r'$|A|^2$')
plt.legend()
plt.grid(linestyle = '--', linewidth = 0.5)
plt.savefig("amp_r.png", dpi = 200)
plt.show()
"""







