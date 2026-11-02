"""Export actual KiCad copper to ordered filaments; calculate free-space DC fields.
Run from this directory's parent, or set --source to the dense-board directory.
Dependencies: numpy scipy matplotlib sexpdata numba; optional magpylib>=5.
"""
from pathlib import Path
import json, math, argparse, collections
import numpy as np
import sexpdata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numba import njit, prange

ROOT=Path(__file__).resolve().parent

def child(n,key):
    return next(x[1:] for x in n if isinstance(x,list) and str(x[0])==key)
def point(n,key,z):
    x,y=map(float,child(n,key));return np.array([x-100,100-y,z])
def extract(board,meta,step=2.5):
    nodes=sexpdata.loads(board.read_text());edges=[];adj=collections.defaultdict(list)
    zm=dict(zip(meta['layer_names'],meta['physical_layer_z_mm'])); used=collections.defaultdict(set)
    for n in nodes:
        if not isinstance(n,list) or str(n[0]) not in ['segment','arc']:continue
        z=zm[str(child(n,'layer')[0])];a=point(n,'start',z);b=point(n,'end',z)
        if str(n[0])=='arc':
            m=point(n,'mid',z);A=a[:2];M=m[:2];B=b[:2]
            ctr=np.linalg.solve(2*np.vstack([M-A,B-A]),[M@M-A@A,B@B-A@A]);r=np.linalg.norm(A-ctr)
            t0,tm,t1=[np.arctan2(*(v-ctr)[::-1]) for v in [A,M,B]]
            sweep=(t1-t0)%(2*np.pi)
            if (tm-t0)%(2*np.pi)>sweep:sweep-=2*np.pi
            ts=t0+np.linspace(0,sweep,max(2,int(np.ceil(abs(sweep)/np.deg2rad(step)))+1))
            pts=np.c_[ctr[0]+r*np.cos(ts),ctr[1]+r*np.sin(ts),np.full(len(ts),z)];pts[0]=a;pts[-1]=b
        else:pts=np.array([a,b])
        edges.append(pts)
        for q in [a,b]:used[tuple(np.round(q[:2],6))].add(z)
    # Only barrel portions carrying the series current; dangling plated stubs omitted.
    for n in nodes:
        if not isinstance(n,list) or str(n[0])!='via':continue
        q=point(n,'at',0);zz=sorted(used[tuple(np.round(q[:2],6))]);assert len(zz)==2,(q,zz)
        edges.append(np.array([[*q[:2],zz[0]],[*q[:2],zz[1]]]))
    for i,p in enumerate(edges):
        adj[tuple(np.round(p[0],6))].append(i);adj[tuple(np.round(p[-1],6))].append(i)
    ends=[q for q,v in adj.items() if len(v)==1];assert len(ends)==2 and all(len(v) in [1,2] for v in adj.values())
    q=max(ends,key=lambda p:p[1]);seen=set();parts=[]
    while True:
        todo=[i for i in adj[q] if i not in seen]
        if not todo:break
        assert len(todo)==1;i=todo[0];seen.add(i);p=edges[i]
        if tuple(np.round(p[0],6))!=q:p=p[::-1]
        parts.append(p if not parts else p[1:]);q=tuple(np.round(p[-1],6))
    assert len(seen)==len(edges)
    v=np.vstack(parts)*1e-3
    assert np.all(np.linalg.norm(np.diff(v,axis=0),axis=1)>1e-12)
    return v

@njit(parallel=True,cache=True)
def field(v,obs,current):
    out=np.zeros((len(obs),3))
    for j in prange(len(obs)):
        for i in range(len(v)-1):
            ax=obs[j,0]-v[i,0];ay=obs[j,1]-v[i,1];az=obs[j,2]-v[i,2]
            bx=obs[j,0]-v[i+1,0];by=obs[j,1]-v[i+1,1];bz=obs[j,2]-v[i+1,2]
            aa=math.sqrt(ax*ax+ay*ay+az*az);bb=math.sqrt(bx*bx+by*by+bz*bz)
            den=aa*bb*(aa*bb+ax*bx+ay*by+az*bz)
            if den>1e-40:
                f=1e-7*current*(aa+bb)/den
                out[j,0]+=(ay*bz-az*by)*f;out[j,1]+=(az*bx-ax*bz)*f;out[j,2]+=(ax*by-ay*bx)*f
    return out

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--source',type=Path,default=ROOT/'source_boards');args=parser.parse_args()
    checks=[];figall,allax=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for idx,folder in enumerate(sorted(args.source.glob('DENSE_*'))):
        meta=json.loads((folder/'parameters_and_checks.json').read_text());name=folder.name;out=ROOT/name;out.mkdir(exist_ok=True)
        board=folder/(name+'.kicad_pcb');v=extract(board,meta);closed=np.vstack([v,v[0]])
        np.savetxt(out/'vertices_open_m.txt',v,fmt='%.12g',header='x_m y_m z_m; IN to OUT, physical PCB only')
        np.savetxt(out/'vertices_closed_m.txt',closed,fmt='%.12g',header='x_m y_m z_m; IN to OUT plus illustrative local external return')
        for label,I in [('1A',1.),('10G',meta['current_10G_A']),('20G',meta['current_20G_A'])]:
            # vuthalab biot-savart 0.3 format: centimeters, amperes, comma delimiter, no header.
            np.savetxt(out/f'biot_savart_{label}_closed_cm.txt',np.c_[closed*100,np.full(len(closed),I)],delimiter=',',fmt='%.12g')
        np.savetxt(out/'biot_savart_1A_open_cm.txt',np.c_[v*100,np.ones(len(v))],delimiter=',',fmt='%.12g')
        h=meta['target_height_mm'];I=meta['current_20G_A'];q=np.array([0.,0.,h/1000]);d=1e-6
        obs=np.vstack([q,q+np.eye(3)*d,q-np.eye(3)*d]);b=field(closed,obs,I);J=((b[1:4]-b[4:7])/(2*d)).T
        fine=extract(board,meta,step=.5);fine=np.vstack([fine,fine[0]]);bf=field(fine,obs,I);Jf=((bf[1:4]-bf[4:7])/(2*d)).T
        err=np.max(abs(J-Jf))/np.max(abs(Jf));assert err<.001
        check={'name':name,'open_vertices':len(v),'current_A':I,'B_at_target_G':(b[0]*1e4).tolist(),'gradient_tensor_G_per_cm':(J*100).tolist(),'gradient_discretization_relative_error':float(err),'arc_step_degrees':2.5}
        try:
            import magpylib as magpy
            src=magpy.current.Polyline(current=I,vertices=closed)
            test=src.getB(obs[:2]);check['magpylib_max_abs_difference_T']=float(np.max(abs(test-b[:2])));assert np.allclose(test,b[:2],rtol=2e-6,atol=1e-12)
        except ImportError:pass
        checks.append(check)
        # Local x-z section, above PCB; uncorrected coil field, no mirror or bias.
        x=np.linspace(-4,4,101);z=np.linspace(max(.5,h-4),h+4,101);X,Z=np.meshgrid(x,z)
        p=np.c_[X.ravel(),np.zeros(X.size),Z.ravel()]*1e-3;B=field(closed,p,I).reshape(*X.shape,3)*1e4
        zz=np.linspace(max(.5,h-5),h+7,241);az=field(closed,np.c_[np.zeros(len(zz)),np.zeros(len(zz)),zz]*1e-3,I)*1e4
        xx=np.linspace(-4,4,161);rx=field(closed,np.c_[xx,np.zeros(len(xx)),np.full(len(xx),h)]*1e-3,I)*1e4
        ry=field(closed,np.c_[np.zeros(len(xx)),xx,np.full(len(xx),h)]*1e-3,I)*1e4
        f,ax=plt.subplots(2,2,figsize=(11,8),layout='constrained')
        mag=np.linalg.norm(B,axis=2);im=ax[0,0].contourf(X,Z,mag,levels=25,cmap='viridis');f.colorbar(im,ax=ax[0,0],label='|B| (G)')
        ax[0,0].streamplot(x,z,B[:,:,0],B[:,:,2],color='white',density=.9,linewidth=.6);ax[0,0].plot(0,h,'rx',ms=7);ax[0,0].set(xlabel='x (mm)',ylabel='z (mm)',title='Full routed coil field: y = 0')
        for k,label in enumerate(['Bx','By','Bz']):ax[0,1].plot(zz,az[:,k],label=label)
        ax[0,1].axvline(h,color='gray',ls=':');ax[0,1].set(xlabel='z (mm)',ylabel='B (G)',title='On optical axis');ax[0,1].legend()
        ax[1,0].plot(xx,rx[:,0],label='Bx along x');ax[1,0].plot(xx,ry[:,1],label='By along y');ax[1,0].set(xlabel='transverse displacement (mm)',ylabel='B (G)',title=f'Transverse cuts at z = {h:g} mm');ax[1,0].legend()
        ax[1,1].plot(zz,np.gradient(az[:,2],zz)*10,label='dBz/dz');ax[1,1].axvline(h,color='gray',ls=':');ax[1,1].set(xlabel='z (mm)',ylabel='gradient (G/cm)',title='Axial gradient');ax[1,1].legend()
        for a in ax.flat:a.grid(alpha=.2)
        f.suptitle(f'{h:g} mm design | I = {I:.4f} A | approximately 20 G/cm axial\nFree-space centerline model; actual copper, vias and backside links; illustrative external return',fontsize=12)
        f.savefig(out/'field_plots_20G.png',dpi=170);plt.close(f)
        a=allax.flat[idx];im=a.contourf(X,Z,mag,levels=np.linspace(0,12,25),extend='max',cmap='viridis');a.streamplot(x,z,B[:,:,0],B[:,:,2],color='white',density=.65,linewidth=.6);a.plot(0,h,'rx');a.set(xlabel='x (mm)',ylabel='z (mm)',title=f'h = {h:g} mm, I = {I:.3f} A');figall.colorbar(im,ax=a,label='|B| (G)')
        np.savez_compressed(out/'field_maps_20G.npz',x_mm=x,z_mm=z,B_G=B,current_A=I)
        np.savetxt(out/'axis_field_20G.csv',np.c_[zz,az,np.gradient(az[:,2],zz)*10],delimiter=',',header='z_mm,Bx_G,By_G,Bz_G,dBz_dz_G_per_cm',comments='')
        (out/'metadata.json').write_text(json.dumps(meta,indent=2))
        print(name,json.dumps(check),flush=True)
    figall.suptitle('Routed PCB MOT coils: free-space DC field, y = 0\n20 G/cm axial setting; no mirror, ambient field or compensation coils',fontsize=13)
    figall.savefig(ROOT/'field_comparison.png',dpi=170)
    (ROOT/'validation.json').write_text(json.dumps(checks,indent=2))
if __name__=='__main__':main()
