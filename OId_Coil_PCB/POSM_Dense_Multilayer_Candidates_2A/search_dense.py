from pathlib import Path
import numpy as np,json
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parent
MU=4*np.pi*1e-7;rho=1.724e-8
W=.3;P=.6008;T=.07

def stack(k):
 L=next(l for l in [2,4,6,8] if l>=k+1)
 thickness={2:.25,4:1.2,6:1.6,8:2.0}[L]
 z=-np.arange(k)*(thickness-T)/(L-1)
 return L,thickness,z

def coeff(a,n,h,zs):
 a=a*.001;b=a+n*P*.001;h=(h-zs)*.001
 def f(r):return np.arcsinh(r/h)-r/np.sqrt(r*r+h*h)
 def d(r):return r**3/(r*r+h*h)**1.5
 K=MU/(2*P*.001)*np.sum(f(b)-f(a))
 D=-MU/(2*P*.001)*np.sum((d(b)-d(a))/h)
 return K,D

results=[];options=[]
for h in [6.,8.,10.,12.]:
 choices=[]
 for body in [35.,45.,55.,65.,75.,85.]:
  for k in [2,3,4,5,6,7]:
   L,th,zs=stack(k);maxr=body/2-.8-W/2
   for a in [2.6,3.5,5.]:
    for n1 in range(10,int((maxr-a-1.2)/P)-9):
     k1,d1=coeff(a,n1,h,zs);lo=a+n1*P+W+1.0
     for n2 in range(10,int((maxr-lo)/P)+1):
      hi=maxr-n2*P
      if hi<=lo:continue
      # All outer radii here lie on or largely beyond the axial-field maximum;
      # scan both branches rather than assume monotonicity.
      xx=np.linspace(lo,hi,7);ff=[coeff(q,n2,h,zs)[0]-k1 for q in xx]
      for j in range(6):
       if ff[j]*ff[j+1]>0:continue
       a2=brentq(lambda q:coeff(q,n2,h,zs)[0]-k1,xx[j],xx[j+1],xtol=1e-10)
       _,d2=coeff(a2,n2,h,zs);g=abs(d1-d2)*100
       if g<=0:continue
       i20=20/g
       if i20>1.8:continue
       length=k*2*np.pi*(n1*(a+n1*P/2)+n2*(a2+n2*P/2))*.001
       r=rho*length/(W*.001*T*.001);power=i20*i20*r
       d=dict(target_height_mm=h,body_size_mm=body,coil_layers=k,copper_layers=L,board_thickness_mm=th,layer_z_mm=zs.tolist(),trace_width_mm=W,clear_gap_mm=.3,pitch_mm=P,copper_thickness_mm=T,inner_start_radius_mm=a,outer_start_radius_mm=a2,turns_inner=n1,turns_outer=n2,current20_winding_A=i20,current10_winding_A=i20/2,gradient_per_A_G_per_cm=g,winding_R20_ohm=r,winding_P_at20G_W=power,intercoil_gap_mm=a2-(a+n1*P)-W,outer_end_radius_mm=a2+n2*P)
       choices.append(d)
  # Stop growing footprint once at least one feasible body has <=10 W at 20G.
  if choices and min(d['winding_P_at20G_W'] for d in choices)<=10:break
 if not choices:
  print(h,'NO FEASIBLE',flush=True);continue
 # Prefer low loss, with a modest explicit penalty on footprint and layers.
 # Report lowest-loss design, not an asserted global optimum.
 best=min(choices,key=lambda d:d['winding_P_at20G_W'])
 results.append(best)
 options.append({'height':h,'feasible_sampled_roots':len(choices),'smallest_feasible_body':min(d['body_size_mm'] for d in choices),'lowest_loss':best})
 print(h,json.dumps(best),flush=True)
(ROOT/'dense_parameters.json').write_text(json.dumps(results,indent=2)+'\n')
(ROOT/'search_summary.json').write_text(json.dumps(options,indent=2)+'\n')
