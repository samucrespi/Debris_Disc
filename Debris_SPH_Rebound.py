#----------------------------------------------------------------------
#----------------------------------------------------------------------
#
# Debris_SPH_Rebound.py
# By S. Crespi, Jan 2022
# Version 3.5
#
# This algorithm is a modfied version of SPH_Rebound.py (v3.3).
#  Its porpouse is to include the debris from a specific collision in 
#   the simulation and from there to evolve the system.
#  In order to keep the code fast, the debris interact with the other
#   bodies (and vice versa) but not with other debris. This is achieved
#   with:
#	 sim.N_active = sim.N-N_debris
#	 sim.testparticle_type = 1
#
# NB: all the units are in "rebound" units (G=1) except where otherwise
#     specified
#
# Version note:
#  - v3.4: include information about the nature of the collding particle
#          (1: planet, 0: debris) in the output files
#  - v3.5: bug fixed - crashing when removing the last particle
#
#----------------------------------------------------------------------
#----------------------------------------------------------------------

#************************   SIMULATION SETUP   ************************

# starting time	( 'start' or 't1000', for example, t0=1000yr)
t0='t100000'

# integration time in years
Dt = 3000
dt = 8.e-3

# checkpoints
Nsaves=3

# scenario
scenario='2_rnBurger20_eJS1'

# number of Gas Giants
NGG=2

#************************   BOOLEAN OPTIONS   *************************

save_checkpoints=True		# generates a snapshot file at each checkpoint
new_sim_warning=False		# requires user input to delete old files

#***************************   LIBRARIES   ****************************

import rebound
import reboundx
import numpy as np
import glob

#***************************   CONSTANTS   ****************************

from CONSTANTS import *

#----------------------------------------------------------------------
#----------------------------------------------------------------------

#****************************   CLASSES   *****************************
		
class SPHcol:
	"""SPH collision"""
	def __init__(self,line):
		val=line.split()

		# collision parameters
		self.id=int(val[0])		# collision index (related to code)
		self.code=val[1]		# code with base 544322 for (v0,alpha,mtot,gamma,wt,wp)
		self.params=np.asarray(val[2:8],dtype=float) #(v0[vesc],alpha[deg],mtot[MSUN],gamma,wt,wp)
		
		# fragmented mass [mtot]
		self.mfr=float(val[9])
		
		# surviving bodies
		self.Nbig=int(val[8])	# number of surviving bodies
		
		# largest bodies
		if self.Nbig==-1: self.crashed=True
		else: self.crashed=False
		self.largest=[]

		for i in [10,18]:
			r=np.asarray(val[i:i+3],dtype=float)	# location wrt CoM in sph.coor. ([AU],[rad],[rad])
			v=np.asarray(val[i+3:i+6],dtype=float)	# velocity wrt CoM in sph.coor. ([AU/yr/2pi],[rad],[rad])
			m=float(val[i+6])						# mass [mtot]
			w=float(val[i+7])						# water [wtot]
			self.largest.append([r,v,m,w])
		
		# Perfect Merging
		if  self.largest[1][2]==-1. and not self.crashed: self.PM=True
		else: self.PM=False


#***************************   FUNCTIONS   ****************************

#-------------------
# Collisions solver
#-------------------

def collision_solver(sim_pointer, collision):
	sim = sim_pointer.contents
	
	# get collision parameters
	p1,p2=ps[collision.p1],ps[collision.p2]
	
	active=True
	# collision with Sun
	if collision.p1==0:
		ps[0].m+=p2.m
		if collision.p2>=sim.N_active: active=False
		write_rmv_file(sim.t,'SUN_col',p2.m,p2.params['wf'],active)
		return 2
	if collision.p2==0:
		ps[0].m+=p1.m
		if collision.p1>=sim.N_active: active=False
		write_rmv_file(sim.t,'SUN_col',p1.m,p1.params['wf'],active)
		return 1

	# collision with Gas Giants
	if collision.p1<=NGG:
		merge(p1,p2)
		if collision.p2>=sim.N_active: active=False
		write_rmv_file(sim.t,'GG{}_col'.format(collision.p1),p2.m,p2.params['wf'],active)
		return 2
	if collision.p2<=NGG:
		merge(p2,p1)
		if collision.p1>=sim.N_active: active=False
		write_rmv_file(sim.t,'GG{}_col'.format(collision.p2),p1.m,p1.params['wf'],active)
		return 1
	
	# collisions with debris
	if collision.p1>=sim.N_active:
		if collision.p2>=sim.N_active: active=False
		write_debris_coll_file(sim.t,p2.m,p2.params['wf'],p1.m,p1.params['wf'],active,False)
		merge(p2,p1)
		return 1
	if collision.p2>=sim.N_active:
		write_debris_coll_file(sim.t,p1.m,p1.params['wf'],p2.m,p2.params['wf'],active,False)
		merge(p1,p2)
		return 2
	
	#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	---------------	    	/\
	xrel=np.asarray([p1.x-p2.x,p1.y-p2.y,p1.z-p2.z])		# >	    BUG      <	   /  \
	if np.sqrt(xrel.dot(xrel))>2.*(p1.r+p2.r): return 0		# >	  DETECTED   <	  /	!! \
	#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	----------------	 /______\

	indeces=[collision.p1,collision.p2]	
	coll_p=get_coll_params(p1,p2)
	
	# save snapshot at the collision time
	collision_snapshot(ps,indeces)
	
	#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
	print('\n{}'.format('-'*80))
	print('Collision parameters: ')
	print('  - v0 = {:.2f}'.format(coll_p[0]))
	print('  - alpha = {:.1f}'.format(coll_p[1]))
	print('  - mtot = {:.2f} [MEAR]'.format(coll_p[2]/MEAR))
	print('  - m1/m2 = {:.2f}'.format(coll_p[3]))
	print('  - wf1 = {:.1f} [%]'.format(coll_p[4]))
	print('  - wf2 = {:.2f} [%]'.format(coll_p[5]))
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^	
	
	# collision CoM position and velocity
	x1,v1,m1,th1,R1=np.asarray([p1.x,p1.y,p1.z]),np.asarray([p1.vx,p1.vy,p1.vz]),p1.m,p1.theta,p1.r
	x2,v2,m2,th2,R2=np.asarray([p2.x,p2.y,p2.z]),np.asarray([p2.vx,p2.vy,p2.vz]),p2.m,p2.theta,p2.r
	
	mtot=m1+m2
	xCoM=(x1*m1+x2*m2)/mtot
	vCoM=(v1*m1+v2*m2)/mtot
	rcol=np.sqrt(xCoM.dot(xCoM))
	thcol=np.arctan2(xCoM[1],xCoM[0])		# coll point projected angle on the x-y plane
	inccol=np.pi/2.-np.arccos(xCoM[2]/rcol)	# coll point "inclination"

	# interpolate SPH table and find the 2 largest bodies
	largest = interpolate_SPHtable(coll_p)

	# get the surviors
	Nbig = get_Nbig(largest[0][2],largest[1][2],coll_p[3])
	survivors = []
	for i in range(Nbig): survivors.append(largest[i])
	
	# in case of a single survivor put it at the collision location
	if Nbig==1: survivors[0][0][0]=0.
	
	# check water consevation
	wtot_mtot=(coll_p[4]+coll_p[3]*coll_p[5])/(1.+coll_p[3])/100.
	if Nbig>0: survivors=water_conservation(survivors,wtot_mtot)
	
	# get mass and water of the fragments
	mfr,mwfr = 1.,1.
	for surv in survivors:
		mfr-=surv[2]
		mwfr-=surv[3]
	
	#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
	mwtot=mtot*wtot_mtot
	print('\n Mass and Water check:')
	print('--- before collision ---')
	print(' TOT: m={:.2f} MEAR'.format(mtot/MEAR))
	print('      w={:.2f} MEAR (wf:{:.2f}%)'.format(mwtot/MEAR,100.*wtot_mtot))
	print(' - T: m={:.2f} MEAR ({:.2f}%)'.format(mtot/(coll_p[3]+1.)/MEAR,1./(coll_p[3]+1.)))
	print('      w={:.3f} MEAR ({:.2f}% from this body with wf:{:.2f}%)'.format(coll_p[4]*mtot*1./(coll_p[3]+1.)/MEAR/100,coll_p[4]*1./(coll_p[3]+1.)/wtot_mtot/100.,coll_p[4]))
	print(' - P: m={:.2f} MEAR ({:.2f}%)'.format(mtot*coll_p[3]/(coll_p[3]+1.)/MEAR,coll_p[3]/(coll_p[3]+1.)))
	print('      w={:.3f} MEAR ({:.2f}% from this body with wf:{:.2f}%)'.format(coll_p[5]*mtot*coll_p[3]/(coll_p[3]+1.)/MEAR/100,coll_p[5]*coll_p[3]/(coll_p[3]+1.)/wtot_mtot/100.,coll_p[5]))
	print('--- after collision ---')
	for i in range(Nbig):
		print(' - S{}: m={:.4f} MEAR ({:.2f}%)'.format(i+1,mtot*survivors[i][2]/MEAR,survivors[i][2]))
		print('       w={:.3f} MEAR ({:.2f}% - wf:{:.2f}%)'.format(survivors[i][3]*mwtot/MEAR,survivors[i][3],survivors[i][3]*wtot_mtot*100/survivors[i][2]))
	print(' - fr: m={:.4f} MEAR ({:.4f}%)'.format(mtot*mfr/MEAR,mfr))
	if mfr>0: print('       w={:.4f} MEAR ({:.4f}% - wf:{:.3f}%)'.format(mwfr*mwtot/MEAR,mwfr,mwfr*wtot_mtot*100/mfr))
	print(' ')
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^	

	# change soc from spherical to cartisian
	for i in range(Nbig):
		survivors[i][0]=from_sph_to_cart(survivors[i][0])
		survivors[i][1]=from_sph_to_cart(survivors[i][1])

	# convert from water mass to water fraction for the surviving bodies and fragments
	for i in range(Nbig): survivors[i][3]=100.*survivors[i][3]*wtot_mtot/survivors[i][2]
	if mfr==0.: wffr=0.
	else: wffr=100.*mwfr*wtot_mtot/mfr
			
	# converting mass fraction to real mass for the surviving bodies and fragments
	for i in range(Nbig): survivors[i][2]=survivors[i][2]*mtot
	mfr=mfr*mtot

	# check for possible reaccretion
	if Nbig==2:
		drel=survivors[0][0]-survivors[1][0]
		vrel=survivors[0][1]-survivors[1][1]
		vrel2=vrel.dot(vrel)
		vesc2=2.*(survivors[0][2]+survivors[1][2])/np.sqrt(drel.dot(drel))
		if vrel2<vesc2:
			ms1=survivors[0][2]
			ms2=survivors[1][2]
			m=ms1+ms2
			x=np.zeros(3)
			v=(survivors[0][1]*ms1+survivors[1][1]*ms2)/m
			wf=(survivors[0][3]*ms1+survivors[1][3]*ms2)/m
			Nbig=1
			survivors=[[x,v,m,wf]]
	
	
	# back-tracing the surviving bodies to t=t_coll
	if Nbig==2: back_tracing_remn(survivors,d=rcol,btd=min_btd)


	# from SPH SoC to Rebound SoC and update particles
	if Nbig>0:
		chi,psi,h=angle_SPH_Rebound([x1,v1,m1],[x2,v2,m2],R1,R2)

		for i in range(Nbig):
			# rotate SoC and move to CoM
			if h[2]>0:
				survivors[i][0]=Ry(survivors[i][0],np.pi)
				survivors[i][1]=Ry(survivors[i][1],np.pi)
			survivors[i][0]=Rz(Rx(survivors[i][0],np.pi/2.-chi),psi)+xCoM
			survivors[i][1]=Rz(Rx(survivors[i][1],np.pi/2.-chi),psi)+vCoM
			
			# update parameters
			sim.particles[indeces[i]].x=survivors[i][0][0]
			sim.particles[indeces[i]].y=survivors[i][0][1]
			sim.particles[indeces[i]].z=survivors[i][0][2]
			sim.particles[indeces[i]].vx=survivors[i][1][0]
			sim.particles[indeces[i]].vy=survivors[i][1][1]
			sim.particles[indeces[i]].vz=survivors[i][1][2]
			sim.particles[indeces[i]].m=survivors[i][2]
			sim.particles[indeces[i]].params['wf']=survivors[i][3]
			sim.particles[indeces[i]].r=get_radius(survivors[i][2],survivors[i][3])
			
		#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv	
			print('\n Survivor ',i+1)
			print('   m: {:.3f} [MEAR]'.format(survivors[i][2]/MEAR))
			print('   wf: {:.2f} [%]'.format(survivors[i][3]))
			print('   r: ',survivors[i][0],' [AU]')
			print('   v: ',survivors[i][1],' [AU/yr/2pi]')
			print('   a: {:.3f} [AU]'.format(sim.particles[indeces[i]].a))
			print('   e: {:.4f} '.format(sim.particles[indeces[i]].e))
			print('   inc: {:.4f} [pi]'.format(sim.particles[indeces[i]].inc/np.pi))
		print()
		#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^	
	
	# save collisions on file
	save_collision(coll_p,xCoM,vCoM,chi,psi,h,survivors,mfr,wffr)

	if Nbig==0: return 3
	if Nbig==1: return 2	
	if Nbig==2: return 0
	

def merge(pa,pb):
	# merge the second body into the first one
	mab=pa.m+pb.m
	pa.x=(pa.m*pa.x+pb.m*pb.x)/mab
	pa.y=(pa.m*pa.y+pb.m*pb.y)/mab
	pa.z=(pa.m*pa.z+pb.m*pb.z)/mab
	pa.vx=(pa.m*pa.vx+pb.m*pb.vx)/mab
	pa.vy=(pa.m*pa.vy+pb.m*pb.vy)/mab
	pa.vz=(pa.m*pa.vz+pb.m*pb.vz)/mab
	pa.m=mab
	pa.params['wf']=(pa.m*pa.params['wf']+pb.m*pb.params['wf'])/mab

def water_conservation(survs,wtot_mtot):
	# m_water<=mass
	mfr,wfr=1.,1.
	for i in range(len(survs)):
		if survs[i][3]>survs[i][2]/wtot_mtot: survs[i][3]=survs[i][2]/wtot_mtot
		mfr=mfr-survs[i][2]
		wfr=wfr-survs[i][3]
	# m_water_fr<=mass_fr
	if wfr>mfr/wtot_mtot:
		C=(1.-mfr/wtot_mtot)/(1.-wfr)
		for i in range(len(survs)): survs[i][3]=C*survs[i][3]
	return survs

def get_coll_params(p1,p2):
	m1,m2=p1.m,p2.m
	xrel=np.asarray([p1.x-p2.x,p1.y-p2.y,p1.z-p2.z])
	vrel=np.asarray([p1.vx-p2.vx,p1.vy-p2.vy,p1.vz-p2.vz])
	xmod=np.sqrt(xrel.dot(xrel))
	vmod=np.sqrt(vrel.dot(vrel))
	
	mtot=m1+m2
	gamma=min(m1/m2,m2/m1)
	vesc=np.sqrt(2.*mtot/xmod)
	v0=vmod/vesc
	alpha=np.arccos(abs(np.dot(xrel/xmod,vrel/vmod)))*180/np.pi #[deg]
	if m1>m2: wft,wfp=p1.params['wf'],p2.params['wf']
	else: wft,wfp=p2.params['wf'],p1.params['wf']

	return v0,alpha,mtot,gamma,wft,wfp

def get_Nbig(m1,m2,g):
	tar=1/(1.+g)
	pro=g*tar
	# CC
	if m1<tar*0.1: return 0		# Catastrophic Collision
	elif m2<pro*0.1: return 1	# Projectile accretion or destruction
	else: return 2				# Hit-and-run

def back_tracing_remn(ss,d=1.,btd=3.):
	# evolves backward in time the two surviving bodies untill they are at distance btd*R_Hill_mutual

	m1,m2=ss[0][2],ss[1][2]
	rmin=btd*Hill_R_mutual(m1,m2,d)
	xs1,xs2=ss[0][0],ss[1][0]
	r0=xs2-xs1
	dist=np.sqrt(r0.dot(r0))
	
	#print('\n Distance at frm 300: {}\n  required: {}'.format(dist,rmin))
	
	if dist>rmin:
		vs1,vs2=ss[0][1],ss[1][1]
		v0=vs2-vs1
		vrel=np.sqrt(v0.dot(v0))
		
		# check the direction of the angular momentum
		if np.cross(r0,v0)[2]<0.: #solve the respective problem with hz>0
			xs1=Ry(xs1,np.pi)
			vs1=Ry(vs1,np.pi)
			xs2=Ry(xs2,np.pi)
			vs2=Ry(vs2,np.pi)
			r0=xs2-xs1
			v0=vs2-vs1
			h_neg=True
		else: h_neg=False

		# orbital elements
		mtot=m1+m2
		a=1./(vrel*vrel/mtot-2./dist)
		h=np.cross(r0,v0)
		h2=h.dot(h)
		e=np.sqrt(1.+h2/mtot/a)
		f0=np.arccos(((a*(e*e-1)/dist)-1.)/e)
		
		print('\na = ',a,' [AU]')
		print('e = ',e)
		print('f_0 = ',f0/np.pi,' [pi]')
		
		# rotating so that the orbit is in the x-y plane
		thz=np.arctan(-h[1]/h[0])	# rotation angle around z-axis
		hy0=Rz(h,thz)
		thy=np.arctan(-hy0[0]/hy0[2])	# rotation angle around y-axis
		
		# rotating so that f=0 for y=0 (omega=0)
		th0=np.arctan2(Ry(Rz(r0,thz),thy)[1],Ry(Rz(r0,thz),thy)[0])
		omega=th0-f0
		R0=Rz(Ry(Rz(r0,thz),thy),-omega)
		V0=Rz(Ry(Rz(v0,thz),thy),-omega)
		
		# check the minimum distance possible between the two bodies
		min_dist=a*(e-1.)
		if min_dist>rmin:	# back-tracing down to f_d=0
			fd=0.
			Rd=np.asarray([min_dist,0.,0.])
			Vd=np.asarray([0.,1.,0.])*np.sqrt(mtot*((2./min_dist)+(1./a)))
			
			#back-tracing time
			F0=np.arccosh((e+np.cos(f0))/(1.+e*np.cos(f0)))
			Dt=(e*np.sinh(F0)-F0)/np.sqrt(mtot/np.power(a,3.))

		else:		# back-tracing down to d=rmin
			fd=np.arccos(((a/rmin)*(e*e-1.)-1.)/e)
			cfd,sfd=np.cos(fd),np.sin(fd)
			phid=np.arctan2(1.+e*cfd,e*sfd)
			Rd=np.asarray([cfd,sfd,0.])*rmin
			Vd=np.asarray([np.cos(fd+phid),np.sin(fd+phid),0.])*np.sqrt(mtot*((2./rmin)+(1./a)))

			#back-tracing time
			F0=np.arccosh((e+np.cos(f0))/(1.+e*np.cos(f0)))
			Fd=np.arccosh((e+np.cos(fd))/(1.+e*np.cos(fd)))
			Dt=((e*np.sinh(F0)-F0)-(e*np.sinh(Fd)-Fd))/np.sqrt(mtot/np.power(a,3.))
		
		print('th_0 = ',th0/np.pi,' [pi]')
		print('omega = ',omega/np.pi,' [pi]')
		print('f_d = ',fd/np.pi,' [pi]')
		print('Dt:',Dt,' [AU/yr/2pi]')
		print()
		
		# get the CoM and move to barycenteric coordinates
		xCoM=m2*Rd/mtot
		vCoM=m2*Vd/mtot
		Rd1=-xCoM
		Vd1=-vCoM
		Rd2=Rd-xCoM
		Vd2=Vd-vCoM	
		
		# rotate back to the SPH simualtion coordinate system
		rd1=Rz(Ry(Rz(Rd1,omega),-thy),-thz)
		vd1=Rz(Ry(Rz(Vd1,omega),-thy),-thz)
		rd2=Rz(Ry(Rz(Rd2,omega),-thy),-thz)
		vd2=Rz(Ry(Rz(Vd2,omega),-thy),-thz)
		
		# correct for the 2-survivors-CoM motion
		xCoM=(m1*xs1+m2*xs2)/mtot
		vCoM=(m1*vs1+m2*vs2)/mtot
		xCoMd=xCoM-vCoM*Dt
		rd1=rd1+xCoMd
		rd2=rd2+xCoMd
		vd1=vd1+vCoM
		vd2=vd2+vCoM
		
		if h_neg:
			rd1=Ry(rd1,-np.pi)
			vd1=Ry(vd1,-np.pi)
			rd2=Ry(rd2,-np.pi)
			vd2=Ry(vd2,-np.pi)
		
		ss[0][0]=rd1
		ss[0][1]=vd1
		ss[1][0]=rd2
		ss[1][1]=vd2	

			
def Hill_R(m,d,M=1.): return d*np.power(m/M/3.,1./3.)

def Hill_R_mutual(m1,m2,d,M=1.): return d*np.power((m1+m2)/M/3.,1./3.)
	
def from_sph_to_cart(vc):
	return np.asarray([vc[0]*np.cos(vc[1])*np.sin(vc[2]),vc[0]*np.sin(vc[1])*np.sin(vc[2]),vc[0]*np.cos(vc[2])])

def Rx(v,t):	#rotates v around the y-axis through the angle t
	ct,st=np.cos(t),np.sin(t)
	return np.asarray([v[0],ct*v[1]-st*v[2],st*v[1]-ct*v[2]])

def Ry(v,t):	#rotates v around the y-axis through the angle t
	ct,st=np.cos(t),np.sin(t)
	return np.asarray([ct*v[0]+st*v[2],v[1],-st*v[0]+ct*v[2]])
	
def Rz(v,t):	#rotates v around the z-axis through the angle t
	ct,st=np.cos(t),np.sin(t)
	return np.asarray([ct*v[0]-st*v[1],st*v[0]+ct*v[1],v[2]])

def angle_SPH_Rebound(l1,l2,R1,R2):

	# move to SoC of smaller body
	if l1[2]>l2[2]:
		r0=l1[0]-l2[0]
		v0=l1[1]-l2[1]
		c=['b','r']
		R=[R2,R1]
	else:
		r0=l2[0]-l1[0]
		v0=l2[1]-l1[1]	
	
	# rotate the system so that h=[0,0,+h]
	h=np.cross(r0,v0)
	phi=np.arccos(h[0]/np.sqrt(h[0]*h[0]+h[1]*h[1]))
	th=np.arccos(h[2]/np.sqrt(h.dot(h)))
	
	r0=Ry(Rz(r0,-phi),-th)
	v0=Ry(Rz(v0,-phi),-th)

	# orbital elements
	r0_mod=np.sqrt(r0.dot(r0))
	v0_mod=np.sqrt(v0.dot(v0))
	mtot=l1[2]+l2[2]
	a=1./(v0_mod*v0_mod/mtot-2./r0_mod)
	e=np.sqrt(1.+h.dot(h)/mtot/a)
	f0=-np.arccos(((a*(e*e-1)/r0_mod)-1.)/e)
	
	# rotating so that f=0 for y=0 (omega=0)
	th0=np.arctan2(r0[1],r0[0])
	omega=th0-f0
	R0=Rz(r0,-omega)
	V0=Rz(v0,-omega)
	
	# back-tracing the colliding bodies to d=5(R1+R2) and get the velocity angle wrt x-axis (delta)
	rmin=5.*(R1+R2)
	fd=-np.arccos(((a/rmin)*(e*e-1.)-1.)/e)
	cfd,sfd=np.cos(fd),np.sin(fd)
	phid=np.arctan2(1.+e*cfd,e*sfd)

	# getting the y_versor of the SPH frame (velocity versor of the bigger body)
	SPHy=np.asarray([np.cos(fd+phid),np.sin(fd+phid),0.])
	SPHy=Rz(Ry(Rz(SPHy,omega),th),phi)
	chi=np.arccos(SPHy[2])
	psi=np.arctan2(-SPHy[0],SPHy[1])
	
	return chi,psi,h


#-------------------
# SPH table interpolator
#-------------------

def interpolate_SPHtable(params):
	
	#params=[5.2,9.,8.e-9,1.,15.,50.]	#~~~~~~~~~~~~~~~~~ NICE TEST
	
	# get the two closest indeces for each parameters of the SPH collisions and the weight (x-xa)/(xb-xa)
	ind_w=[get_interppoint_indices_and_weight(i,x) for i,x in enumerate(params)]
	
	#for i,x in enumerate(ind_w): print(params[i],'  in  ',allpar[i],'  ->  ',x)	#~~~~~~~~~~~~ prints the interpolation points
	
	# get the SPH simulation codes - interpolation points
	points=[]
	absolute_w=True
	for p1 in ind_w[0][0]:
		for p2 in ind_w[1][0]:
			for p3 in ind_w[2][0]:
				for p4 in ind_w[3][0]:
					for p5 in ind_w[4][0]:
						for p6 in ind_w[5][0]:
							code='{}{}{}{}{}{}'.format(p1,p2,p3,p4,p5,p6)
							nr=get_index_from_code(code)-1
							if SPHcat[nr].code!=code:
							   print('\n>>> WRONG SPH collision selected! <<<\n')
							   sys.exit()
							points.append([SPHcat[nr].largest,SPHcat[nr].code])
	
	for j in range(6): points=[interpolate(points[i*2],points[i*2+1],ind_w[5-j][1]) for i in range(int(len(points)/2))]
	return points[0][0]
	
def get_interppoint_indices_and_weight(i,x):
	# log for mass (i=2)
	if i==2: log=True
	else: log=False
	
	for j,xj in enumerate(allpar[i]):
		if x<xj:
			if j==0: ind=[0,1]
			else: ind=[j-1,j]
			break
		ind=[len(allpar[i])-2,len(allpar[i])-1]
	dx=weight(x,allpar[i][ind[0]],allpar[i][ind[1]],log=log)
	return [ind,dx]

def weight(x,xa,xb,log=False):
	if log: return np.log10(x/xa)/np.log10(xb/xa)
	else: return (x-xa)/(xb-xa)
	
def interpolate(pa,pb,dx):
	#SPH collision code check:
	#print('  -- {} + {}  ->  {}'.format(pa[1],pb[1],pa[1][:-1]))	#~~~~~~~ prints SPH code while interpolating
	if pa[1][:-1]!=pb[1][:-1]:
		print('\n>>> WRONG SPH collision coupling! <<<\n')
		sys.exit()
	
	# --- 1st largest ---
	ya,yb=pa[0][0],pb[0][0]
	
	#solving crashed
	if ya[2]==-1.: return pb[0],pa[1][:-1]
	if yb[2]==-1.: return pa[0],pa[1][:-1]
	
	#solving not crashed
	r=[lin_interpol(ya[0][0],yb[0][0],dx),interpol_angle(ya[0][1],yb[0][1],dx),interpol_angle(ya[0][2],yb[0][2],dx)]
	v=[lin_interpol(ya[1][0],yb[1][0],dx),interpol_angle(ya[1][1],yb[1][1],dx),interpol_angle(ya[1][2],yb[1][2],dx)]
	m=lin_interpol(ya[2],yb[2],dx)
	gwf=lin_interpol(ya[3],yb[3],dx)
	largest=[[r,v,m,gwf]]
	
	# --- 2nd largest ---
	ya,yb=pa[0][1],pb[0][1]
	
	#solving pa:PM
	if ya[2]==-1.:								# pa: PM
		if yb[2]==-1.: largest.append(ya)		#pa & pb PM	
		else:
			if dx<=0.: largest.append(ya)		#extrapolation->PM
			elif dx>=1.: largest.append(yb)		#extrapolation->yb
			else:								#interpolate with a PM
				r=[lin_interpol(0.,yb[0][0],dx),yb[0][1],yb[0][2]]
				v=[lin_interpol(0.,yb[1][0],dx),yb[1][1],yb[1][2]]
				m=lin_interpol(0.,yb[2],dx)
				gwf=yb[3]
				largest.append([r,v,m,gwf])
		adjust_result(largest)
		return largest,pa[1][:-1]
	
	#solving pb:PM
	if yb[2]==-1.:								# pb: PM
		if dx<=0.: largest.append(ya)			#extrapolation->ya
		elif dx>=1.: largest.append(yb)			#extrapolation->PM
		else:									#interpolate with a PM
			r=[lin_interpol(ya[0][0],0.,dx),ya[0][1],ya[0][2]]
			v=[lin_interpol(ya[1][0],0.,dx),ya[1][1],ya[1][2]]
			m=lin_interpol(ya[2],0.,dx)
			gwf=ya[3]
			largest.append([r,v,m,gwf])
		adjust_result(largest)
		return largest,pa[1][:-1]
	
	#solving not PM
	r=[lin_interpol(ya[0][0],yb[0][0],dx),interpol_angle(ya[0][1],yb[0][1],dx),interpol_angle(ya[0][2],yb[0][2],dx)]
	v=[lin_interpol(ya[1][0],yb[1][0],dx),interpol_angle(ya[1][1],yb[1][1],dx),interpol_angle(ya[1][2],yb[1][2],dx)]
	m=lin_interpol(ya[2],yb[2],dx)
	gwf=lin_interpol(ya[3],yb[3],dx)
	largest.append([r,v,m,gwf])
	adjust_result(largest)
	return largest,pa[1][:-1]
	
def lin_interpol(a,b,dx,log=False):
	if log: return np.power(10.,lin_interpol(np.log10(a),np.log10(b),dx))
	return a*(1.-dx)+b*dx

def interpol_angle(a,b,dx):
	return np.arctan2(lin_interpol(np.sin(a),np.sin(b),dx),lin_interpol(np.cos(a),np.cos(b),dx))

def adjust_result(largest):
	# when extrapolating, it can happen that: m1<m2		
	if largest[0][2]<largest[1][2]: largest=[largest[1],largest[0]]

	# ... or m1<0
	if largest[0][2]<0:		# and so m2<0
		largest[0][2]=1./SPHRES		#smallest SPH particle
		largest[1][2]=1./SPHRES		#smallest SPH particle
		
	# ... or m2<0
	if largest[1][2]<0:
		if largest[0][2]>=1. or largest[1][2]==-1.: 	#PM
			largest[0][2]=1.
			largest=[largest[0],[[-1.,-1.,-1.],[-1.,-1.,-1.],-1.,-1.]]
		else: largest[1][2]=1./SPHRES

	# ...or m1+m2>1
	if largest[0][2]+largest[1][2]>1.:
		corr=largest[0][2]+largest[1][2]
		largest[0][2],largest[1][2]=largest[0][2]/corr,largest[1][2]/corr
	
	# ...or PM with m1>1
	if largest[1][2]==-1. and largest[0][2]!=1.: largest[0][2]=1.

	# ... and r<0 and/or v<0
	for i in range(2):
		for j in range(2):
			if largest[i][j][0]<0.: largest[i][j][0]=0.
				
	# .... and gwf<0 and/or gwf>0
	if largest[0][2]==1.: largest[0][3]=1.	#PM
	else:
		for i in range(2):
			if largest[i][3]<0: largest[i][3]=0.
		if largest[0][3]+largest[1][3]>1.:
			corr=largest[0][3]+largest[1][3]
			largest[0][3],largest[1][3]=largest[0][3]/corr,largest[1][3]/corr

def get_radius(m,wf):	# m in Solar masses
	# independent from wf
	# from Chen & Kipping 2017
	C=1.008*REAR
	S=0.279
	R=C*np.power(m/MEAR,S)
	# depending on wf
	# ... still to be done		<<------------
	return R	# radius [AU]


#-------------------
# SPH catalogue
#-------------------

def load_SPHcat(fname):
	f=open(fname,'r')
	lines=f.readlines()
	f.close()
	head=lines.pop(0)
	return [SPHcol(line) for line in lines]

def get_index_from_code(code):
	# returns SPHcol.id given SPH.code
	val=0
	for i in range(len(bases)-1): val=(val+int(code[i]))*(int(bases[i+1]))
	return val+int(code[5])+1


#-------------------
# Debries
#-------------------

def get_collision_outcome(coll_details):
	
	
	print('\n ...getting debris from ',get_SPHdir(coll_details[1:7]))
	debris=np.loadtxt(get_SPHdir(coll_details[1:7])+'/all_aggregates.txt',skiprows=1)
	
	# find number of actual debris
	Ndebris=len(debris)-get_Nbig(debris[0,7],debris[1,7],coll_details[4])

	# adjust units
	move_to_Rebound_units(debris)
	
	# adjust mass
	SPH_totmass=sum(debris[:,6])		# [MSUN]
	debris[:,6]=debris[:,6]*(coll_details[3]*MEAR)/SPH_totmass
	
	# adjust water
	SPH_totwater=sum(debris[:,6]*debris[:,8])		# [MSUN]
	col_totwater=(coll_details[3]*MEAR)*(coll_details[5]+coll_details[4]*coll_details[6])/(1.+coll_details[4])/100	# [MSUN]
	if SPH_totwater>col_totwater: debris[:,8]=debris[:,8]*col_totwater/SPH_totwater
	else: increase_water_where_possible(debris,col_totwater)

	# collision CoM
	xCoM_col=np.asarray(coll_details[7:10])
	vCoM_col=np.asarray(coll_details[10:13])
	chi,psi,hz=coll_details[13:16]

	res=[]
	for i in range(len(debris)):
		x=np.asarray(debris[i,:3])
		v=np.asarray(debris[i,3:6])
	
		# rotate SoC and move to CoM
		if hz>0:
			x=Ry(x,np.pi)
			v=Ry(v,np.pi)
		x=Rz(Rx(x,np.pi/2.-chi),psi)+xCoM_col
		v=Rz(Rx(v,np.pi/2.-chi),psi)+vCoM_col
		res.append([x,v,debris[i,6],debris[i,8]*100.])	
	
	return res,Ndebris

def get_SPHdir(vals):
	
	v0,alpha,mtot,gamma,wt,wp=vals
	mtot=mtot*MEAR	# mtot in MSUN
	# v0
	fname='v'
	mid_v0s=[(v0s[i]+v0s[i+1])/2. for i in range(len(v0s)-1)]
	if v0<mid_v0s[0]: fname+='1.0'
	elif v0<mid_v0s[1]: fname+='1.5'
	elif v0<mid_v0s[2]: fname+='2.0'
	elif v0<mid_v0s[3]: fname+='3.0'
	else: fname+='5.0'
	# alpha
	fname+='_a'
	mid_alphas=[(alphas[i]+alphas[i+1])/2. for i in range(len(alphas)-1)]
	if alpha<mid_alphas[0]: fname+='0'
	elif alpha<mid_alphas[1]: fname+='20'
	elif alpha<mid_alphas[2]: fname+='40'
	else: fname+='60'
	# mtot
	fname+='_m'
	mid_mtots=[2.*np.sqrt(mtots[i]*mtots[i+1]) for i in range(len(mtots)-1)]
	if mtot<mid_mtots[0]: fname+='21'
	elif mtot<mid_mtots[1]: fname+='23'
	elif mtot<mid_mtots[2]: fname+='24'
	else: fname+='25'
	# gamma
	fname+='_g'
	mid_gammas=[(gammas[i]+gammas[i+1])/2. for i in range(len(gammas)-1)]
	if gamma<mid_gammas[0]: fname+='0.1'
	elif gamma<mid_gammas[1]: fname+='0.5'
	else: fname+='1.0'
	# wfts
	fname+='_wt'
	mid_wfts=[(wfts[i]+wfts[i+1])/2. for i in range(len(wfts)-1)]
	if wt<mid_wfts[0]: fname+='10.0'
	else: fname+='20.0'
	# wfps
	fname+='_wp'
	mid_wfps=[(wfps[i]+wfps[i+1])/2. for i in range(len(wfps)-1)]
	if wp<mid_wfps[0]: fname+='10.0'
	else: fname+='20.0'
	
	return glob.glob('SPH_catalogue/*'+fname)[0]

def move_to_Rebound_units(debris):
	for i in range(3): debris[:,i]=debris[:,i]/AU
	for i in range(3): debris[:,3+i]=debris[:,3+i]*YEAR/AU/2./np.pi
	debris[:,6]=debris[:,6]/MSUN

def increase_water_where_possible(debris,wtarget):
	mw_wnot1,mw_w1=0.,0.
	for i in range(len(debris)):
		if debris[i,8]<1.: mw_wnot1+=debris[i,6]*debris[i,8]
		else:
			debris[i,8]=1.
			mw_w1+=debris[i,6]
	R=(wtarget-mw_w1)/mw_wnot1
	for i in range(len(debris)):
		if debris[i,8]<1.: debris[i,8]=debris[i,8]*R
	if max(debris[:,8])>1.: increase_water_where_possible(debris,wtarget)
		

#-------------------
# remove particles
#-------------------

def remove_ps():
	sim.move_to_hel()
	for i in reversed(range(1,sim.N)):
		active=True
		# unbound particles
		if ps[i].e>=1.:
			if i>=sim.N_active: active=False
			write_rmv_file(sim.t,'HYPER',ps[i].m,ps[i].params['wf'],active)
			sim.remove(i)
			continue
		# too far
		if ps[i].a>100.:
			if i>=sim.N_active: active=False
			write_rmv_file(sim.t,'Agt100',ps[i].m,ps[i].params['wf'],active)
			sim.remove(i)
			continue
		if ps[i].a*(1.-ps[i].e)>12.:
			if i>=sim.N_active: active=False
			write_rmv_file(sim.t,'PERIgt12',ps[i].m,ps[i].params['wf'],active)
			sim.remove(i)
			continue
		# too close to the Sun
		if ps[i].a<0.2:
			if i>=sim.N_active: active=False
			write_rmv_file(sim.t,'Alt02',ps[i].m,ps[i].params['wf'],active)
			ps[0].m+=ps[i].m
			sim.remove(i)
			continue
		if ps[i].a*(1.-ps[i].e)<0.03:
			if i>=sim.N_active: active=False
			write_rmv_file(sim.t,'PERIlt003',ps[i].m,ps[i].params['wf'],active)
			ps[0].m+=ps[i].m
			sim.remove(i)
			continue
	sim.move_to_com()


#-------------------
# output files
#-------------------

def initialize_collisions_file():
	f=open(coll_file,'w+')
	f.write('# time [yr]\tv0 [v_esc]\talpha [degrees]\tmtot [MEAR]\tgamma [1]\twft [%]\twfp [%]')
	f.write('\tx_CoM [AU]\ty_CoM [AU]\tz_CoM [AU]\tvx_CoM [AU/yr/2pi]\tvy_CoM [AU/yr/2pi]\tvz_CoM [AU/yr/2pi]')
	f.write('\tchi [rad]\tpsi [rad]\thz [AU^2/yr/2pi]')
	f.write('\tNbig [1]\tm1 [MEAR]\twf1 [%]\tm2 [MEAR]\twf2 [%]\tmfr [MEAR]\twffr [%]\n')
	f.close() 

def initialize_removed_particles_file():
	f=open(rmv_file,'w+')
	f.write('# time [yr]\tTYPE\tm [SUN]\twf [%]\tactive\n')
	f.close()

def initialize_debris_coll_file():
	f=open(debris_coll_file,'w+')
	f.write('# time [yr]\tm1 [SUN]\twf1 [%]\tm2 [SUN]\twf2 [%]\tactive1\tactive2\n')
	f.close()

def save_collision(params,xCoM,vCoM,chi,psi,h,survs,mfr,wffr):
	# position and velocity of the collision's CoM are given in the CoM system of coordinate of the simulation
	v0,alpha,mtot,gamma,wft,wfp=params
	f=open(coll_file,'a')
	f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(sim.t/2./np.pi,v0,alpha,mtot/MEAR,gamma,wft,wfp))
	for val in xCoM: f.write('\t{}'.format(val))
	for val in vCoM: f.write('\t{}'.format(val))
	f.write('\t{}\t{}\t{}'.format(chi,psi,h[2]))
	Nbig=len(survs)
	if Nbig==0: f.write('\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(Nbig,0,0,0,0,mfr/MEAR,wffr))
	if Nbig==1: f.write('\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(Nbig,survs[0][2]/MEAR,survs[0][3],0,0,mfr/MEAR,wffr))
	if Nbig==2: f.write('\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(Nbig,survs[0][2]/MEAR,survs[0][3],survs[1][2]/MEAR,survs[1][3],mfr/MEAR,wffr))
	f.close()

def save_data(t,ps,path):
    sim.move_to_hel()
    f=open('{}/{}_debris.t{}'.format(path,path,int(np.around(t))),'w+')
    f.write('m [MSUN]\tR [AU]\ta [AU]\te [1]\tinc [rad]\tOmega [rad]\tomega [rad]\tM [rad]\twf [%]\tactive\t')
    f.write('x [AU]\ty [AU]\tz [AU]\tvx [AU/yr/2pi]\tvy [AU/yr/2pi]\tvz [AU/yr/2pi]\n')
    line1='{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t'
    line2='{}\t{}\t{}\t{}\t{}\t{}\n'
    f.write(line1.format(ps[0].m,ps[0].r,0.,0.,0.,0.,0.,0.,0.,1.))
    f.write(line2.format(0.,0.,0.,0.,0.,0.))
    for i,p in enumerate(ps[1:]):
        if i<sim.N_active-1: f.write(line1.format(p.m,p.r,p.a,p.e,p.inc,p.Omega,p.omega,p.M,p.params['wf'],1))
        else: f.write(line1.format(p.m,p.r,p.a,p.e,p.inc,p.Omega,p.omega,p.M,p.params['wf'],0))
        f.write(line2.format(p.x,p.y,p.z,p.vx,p.vy,p.vz))
    sim.move_to_com()
    return f
	
def write_rmv_file(t,ty,m,wf,active):
	f=open(rmv_file,'a')
	f.write('{}\t{}\t{}\t{}\t{}\n'.format(t/2./np.pi,ty,m,wf,int(active)))
	f.close()
	
def write_debris_coll_file(t,m1,wf1,m2,wf2,active1,active2):
	f=open(debris_coll_file,'a')
	f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(t/2./np.pi,m1,wf1,m2,wf2,int(active1),int(active2)))
	f.close()
	
def collision_snapshot(ps,indeces):
	previous_files=glob.glob(scenario+'/*.snapshot')
	collN=len(previous_files)+1
	new_snap='{}/collision_{}.snapshot'.format(scenario,collN)
	sim.move_to_hel()
	
	f=open(new_snap,'w+')
	f.write('# m [MSUN]\tR [AU]\ta [AU]\te [1]\tinc [rad]\tOmega [rad]\tomega [rad]\tM [rad]\twf [%]')
	f.write(' x [AU]\ty [AU]\tz [AU]\tvx [AU/yr/2pi]\tvy [AU/yr/2pi]\tvz [AU/yr/2pi]\n')
	f.write('# Colliding bodies\n')
	linea='{}'+'\t{}'*8
	lineb='\t{}'*6+'\n'
	for i in indeces:
		f.write(linea.format(ps[i].m,ps[i].r,ps[i].a,ps[i].e,ps[i].inc,ps[i].Omega,ps[i].omega,ps[i].M,ps[i].params['wf']))
		f.write(lineb.format(ps[i].x,ps[i].y,ps[i].z,ps[i].vx,ps[i].vy,ps[i].vz))
	f.write('# Other bodies\n')
	f.write(linea.format(ps[0].m,ps[0].r,0,0,0,0,0,0,0))
	f.write(lineb.format(0,0,0,0,0,0))
	for i in range(1,len(ps)-1):
		if i==0 or i==indeces[0] or i==indeces[1]: continue
		f.write(linea.format(ps[i].m,ps[i].r,ps[i].a,ps[i].e,ps[i].inc,ps[i].Omega,ps[i].omega,ps[i].M,ps[i].params['wf']))
		f.write(lineb.format(ps[i].x,ps[i].y,ps[i].z,ps[i].vx,ps[i].vy,ps[i].vz))
	f.close()	
	sim.move_to_com()
	
	

#----------------------------------------------------------------------
#----------------------------------------------------------------------


#***************   INITIALIZE SIMULATION OUTPUT FILES   ***************

# old file warning
if new_sim_warning:
	import os
	import glob
	old_files=glob.glob(scenario+'/*.dat')+glob.glob(scenario+'/*.snapshot')
	if old_files!=[]:
		print('\n\n {}\n WARNING!!!\n \n {}\n'.format('~'*50,'~'*50))
		print(' Previos simulation file detected!\n')
		for f in old_files: print(f)
		user=str(input('\nRemove old files? (y,n) '))
		if user=='y':
			for f in old_files: os.remove(f)

# collisions list file
coll_file='{}/coll.dat'.format(scenario)
if glob.glob(coll_file)==[]: initialize_collisions_file()

# removed particles file
rmv_file='{}/rmv.dat'.format(scenario)
if glob.glob(rmv_file)==[]: initialize_removed_particles_file()

# debris collisions file
debris_coll_file='{}/debris_coll.dat'.format(scenario)
if glob.glob(debris_coll_file)==[]: initialize_debris_coll_file()


#***********************   SIMULATION SETUP   *************************
											     
sim = rebound.Simulation()

#integrator options
sim.integrator = "mercurius"
sim.ri_mercurius.hillfac = 3.
sim.dt=dt

#collision and boundary options
sim.collision = "direct"				# The way it looks for collisions
sim.collision_resolve_keep_sorted = 1
sim.boundary = "none"
sim.track_energy_offset = 1

#collision solver
sim.collision_resolve = collision_solver	# Custom collision solver
min_btd = 3.							# minimum back-tracking distance in mutual Hill-radii

#reboundx
rebx = reboundx.Extras(sim)		# add the extra parameter water fraction 'wf'

#SPH catalogue
SPHcat = load_SPHcat('SPH.table')		# load the SPH coll. catalogue from SPH.table

#removing particles
rem_freq=1000				# checks for particles to be romeved every rem_freq time steps

#----------------------------------------------------------------------
#----------------------------------------------------------------------

#*****************************   INPUT   ******************************

start=np.loadtxt('{}/{}_debris.{}'.format(scenario,scenario,t0),skiprows=1)
Ndebries,i=0,0

# --- Sun
sim.add(m=start[0,0],r=start[0,1],hash='Sun')
sim.particles[0].params['wf']=start[0,8]
# --- Active Particles
for val in start[1:]:
	if bool(val[9]):
		i+=1
		sim.add(m=val[0],r=val[1],a=val[2],e=val[3],inc=val[4],Omega=val[5],omega=val[6],M=val[7])
		sim.particles[i].params['wf']=val[8]
		
# --- Debris
for val in start[1:]:
	if not bool(val[9]):
		Ndebries+=1
		i+=1
		sim.add(m=val[0],r=val[1],a=val[2],e=val[3],inc=val[4],Omega=val[5],omega=val[6],M=val[7])
		sim.particles[i].params['wf']=val[8]

# --- "deactivating" the debries
sim.N_active = sim.N-Ndebries
sim.testparticle_type = 1

ps = sim.particles
#sim.move_to_com()
E0 = sim.calculate_energy() # Calculate initial energy 

#----------------------------------------------------------------------
#----------------------------------------------------------------------		

#*****************************   OUTPUT   *****************************

# starting time	( 'start' or 't1000', for example, t0=1000yr)

tmin=float(t0[1:])*np.pi*2.
sim.t=tmin
tend=tmin+Dt*2.*np.pi
times = np.linspace(tmin,tend,Nsaves+1)[1:]

import time
clock0=time.time()
print(' - t={:.0f} yr   N={}  Ndebris={}   partial: {:.2f} s   running: {:.2f} s'.format(np.round(sim.t/2./np.pi),sim.N_active,sim.N-sim.N_active,0,0))
for t in times:
	clock1=time.time()
	
	# removing particles stepping 
	tnext=sim.t+rem_freq*sim.dt
	while tnext<t:
		sim.integrate(tnext)
		remove_ps() 	# remove particles
		tnext+=rem_freq*sim.dt
	
	# saving data stepping 
	sim.integrate(t)
	remove_ps() 	# remove particles

	if save_checkpoints: save_data(sim.t/np.pi/2.,ps,scenario)
	
	clock2=time.time()
	E1 = sim.calculate_energy()
	print(' - t={:.0f} yr   N={}  Ndebris={}   partial: {:.2f} s   running: {:.2f} s   err_rel={:.2e}'.format(np.round(sim.t/2./np.pi),sim.N_active,sim.N-sim.N_active,clock2-clock1,clock2-clock0,(E1-E0)/E0))
	
	

















