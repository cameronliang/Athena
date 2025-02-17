import matplotlib
matplotlib.use('Agg')

import numpy as np

import h5py
import progressbar
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm



# setup latex fonts
#rcParams['text.usetex']=True
#rcParams['text.latex.unicode']=True


#matplotlib.rc('font', family='serif', serif='cm10')
#matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

# The full quantities' names are:
#['rho' 'vel1' 'vel2' 'vel3' 'pgas' 'B1' 'B2' 'B3' 'Er' 'Fr1' 'Fr2' 'Fr3'
# 'Pr11' 'Pr22' 'Pr33' 'Pr12' 'Pr13' 'Pr23' 'Er0' 'Fr01' 'Fr02' 'Fr03'
# 'Sigma_s' 'Sigma_a']


def ReduceData(filename):

# Read attributes and data
  f=h5py.File(filename, 'r')
  level = f.attrs['MaxLevel']
  block_size = f.attrs['MeshBlockSize']
  root_grid_size = f.attrs['RootGridSize']
  nx1 = root_grid_size[0] * 2**level
  nx2 = root_grid_size[1] * 2**level if root_grid_size[1] > 1 else 1
  nx3 = root_grid_size[2] * 2**level if root_grid_size[2] > 1 else 1

  if nx3 > 1:
    dim = 3
  elif nx2 > 1:
    dim = 2
  else:
    dim = 1


  quantities = f.attrs['VariableNames'][:]
  quantities = [str(q) for q in quantities \
        if q != 'x1f' and q != 'x2f' and q != 'x3f']

  data={}

# the coordinates
  for d in range(1,4):
    nx = (nx1,nx2,nx3)[d-1]
    xmin = f.attrs['RootGridX'+repr(d)][0]
    xmax = f.attrs['RootGridX'+repr(d)][1]
    xrat_root = f.attrs['RootGridX'+repr(d)][2]
    if (xrat_root == 1.0):
      data['x'+repr(d)+'f'] = np.linspace(xmin, xmax, nx+1)
    else:
      xrat = xrat_root ** (1.0 / 2**level)
      data['x'+repr(d)+'f'] = \
         xmin + (1.0-xrat**np.arange(nx+1)) / (1.0-xrat**nx) * (xmax-xmin)

# Get metadata describing file layout
  num_blocks = f.attrs['NumMeshBlocks']
  dataset_names = f.attrs['DatasetNames'][:]
  dataset_sizes = f.attrs['NumVariables'][:]
  dataset_sizes_cumulative = np.cumsum(dataset_sizes)
  variable_names = f.attrs['VariableNames'][:]
  levels = f['Levels'][:]
  logical_locations = f['LogicalLocations'][:]
  quantity_datasets = []
  quantity_indices = []
  spec_datasets = []
  spec_indices = []
  for q in quantities:
    var_num = np.where(variable_names == q)[0][0]
    dataset_num = np.where(dataset_sizes_cumulative > var_num)[0][0]
    if dataset_num == 0:
      dataset_index = var_num
    else:
      dataset_index = var_num - dataset_sizes_cumulative[dataset_num-1]
    quantity_datasets.append(dataset_names[dataset_num])
    quantity_indices.append(dataset_index)
    if q == 'rho' or q=='vel1' or q=='vel2' or q=='vel3' or q=='pgas' or q=='B1' or \
                q=='B2' or q=='B3' or q=='Er' or q=='Sigma_s' or q=='Sigma_a' or q=='Fr01' or q=='Fr02':
      spec_datasets.append(dataset_names[dataset_num])
      spec_indices.append(dataset_index)

# get rho, v1, v2, v3, B1, B2, B3, kappa_s, kappa_a


# now add the derived quantities
  quantities.append('Maxwell')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('Reynolds')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('BrBphi')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('BthetaBphi')
  quantity_datasets.append('None')
  quantity_indices.append(0)
  
  quantities.append('rhovrvphi')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('rhovthetavphi')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('vrEr')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('vthetaEr')
  quantity_datasets.append('None')
  quantity_indices.append(0)
  
  quantities.append('rhoPB')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('rhosq')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('PB1')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('PB2')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('PB3')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('PB')
  quantity_datasets.append('None')
  quantity_indices.append(0)


  quantities.append('PBsq')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('kappa_s')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('kappa_a')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('Radacc')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('RhoVr')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('RhoVtheta')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('RhoVphi')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('RhoVout')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('RhoVin')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('Ekin1')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('Ekin2')
  quantity_datasets.append('None')
  quantity_indices.append(0)
  
  quantities.append('Ekin3')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('tgas')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('Fr01Sigma')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('Fr02Sigma')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('Fr01kappa')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('Fr02kappa')
  quantity_datasets.append('None')
  quantity_indices.append(0)

# the quantities
# get the azimuthal averaged data

  for q in quantities:
    data[q] = np.zeros((nx3,nx2,nx1))



  bar = progressbar.ProgressBar(maxval=num_blocks, \
           widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()]).start()


# Go through blocks in data file
  for block_num in range(num_blocks):
    # Extract location information
    block_level = levels[block_num]
    block_location = logical_locations[block_num,:]
   # Calculate scale (number of copies per dimension)
    s= 2**(level-block_level)
  
  # the block size
    radius=f['x1f'][block_num]
    theta=f['x2f'][block_num]
    phi=f['x3f'][block_num]
  
  # get cell center coordinates
    x1v=np.zeros(block_size[0])
    x2v=np.zeros(block_size[1])
    x3v=np.zeros(block_size[2])
  
    for ni in range(block_size[0]):
      x1v[ni]=0.75*(radius[ni+1]**4.0 - radius[ni]**4.0)/(radius[ni+1]**3.0-radius[ni]**3.0);
  
    for nj in range(block_size[1]):
      x2v[nj]=((np.sin(theta[nj+1]) - theta[nj+1] * np.cos(theta[nj+1])) \
            -(np.sin(theta[nj]) - theta[nj] * np.cos(theta[nj]))) \
              / (np.cos(theta[nj]) - np.cos(theta[nj+1]));

    for nk in range(block_size[2]):
      x3v[nk] = 0.5 * (phi[nk+1]+phi[nk])
  
    grid_phi,grid_theta,grid_r=np.meshgrid(x3v,x2v,x1v,indexing='ij')
  

# Calculate fine-level begin indices
    il = block_location[0] * block_size[0] * s
    jl = block_location[1] * block_size[1] * s if dim >= 2 else 0
    kl = block_location[2] * block_size[2] * s if dim >= 3 else 0

# Calculate fine-level end indices
    iu = il + block_size[0] * s
    ju = jl + block_size[1] * s if dim >= 2 else 1
    ku = kl + block_size[2] * s if dim >= 3 else 1

# Calculate fine-level offsets
    io_vals = range(s)
    jo_vals = range(s) if dim >= 2 else (0,)
    ko_vals = range(s) if dim >= 3 else (0,)

          
    rho_data=f[spec_datasets[0]][spec_indices[0],block_num,:]
    v1_data=f[spec_datasets[1]][spec_indices[1],block_num,:]
    vout_data=v1_data.clip(min=0.0)
    vin_data=v1_data.clip(max=0.0)
    v2_data=f[spec_datasets[2]][spec_indices[2],block_num,:]
    v3_data=f[spec_datasets[3]][spec_indices[3],block_num,:]
    pgas_data=f[spec_datasets[4]][spec_indices[4],block_num,:]
    B1_data=f[spec_datasets[5]][spec_indices[5],block_num,:]
    B2_data=f[spec_datasets[6]][spec_indices[6],block_num,:]
    B3_data=f[spec_datasets[7]][spec_indices[7],block_num,:]
    Er_data=f[spec_datasets[8]][spec_indices[8],block_num,:]
    Fr01_data=f[spec_datasets[9]][spec_indices[9],block_num,:]
    Fr02_data=f[spec_datasets[10]][spec_indices[10],block_num,:]
    sigma_s_data=f[spec_datasets[11]][spec_indices[11],block_num,:]
    sigma_a_data=f[spec_datasets[12]][spec_indices[12],block_num,:]
    PB_data=0.5*(np.multiply(B1_data,B1_data)+np.multiply(B2_data,B2_data)+np.multiply(B3_data,B3_data))
    PB1_data=0.5*np.multiply(B1_data,B1_data)   
    PB2_data=0.5*np.multiply(B2_data,B2_data)
    PB3_data=0.5*np.multiply(B3_data,B3_data)
    rhovphi_data=np.multiply(rho_data,v3_data)



# Assign values
    for q,dataset,index in zip(quantities,quantity_datasets,quantity_indices):
      if q=='rhosq':
         oridata=np.multiply(rho_data,rho_data)
      elif q=='PB':
         oridata=PB_data
      elif q=='PB1':
         oridata=PB1_data
      elif q=='PB2':
         oridata=PB2_data
      elif q=='PB3':
         oridata=PB3_data
      elif q=='PBsq':
         oridata=np.multiply(PB_data,PB_data)
      elif q=='kappa_s':
         oridata=np.divide(sigma_s_data,rho_data)
      elif q=='kappa_a':
         oridata=np.divide(sigma_a_data,rho_data)
      elif q=='Maxwell':
         bxdata=np.multiply(B2_data,np.cos(grid_theta))+np.multiply(B1_data,np.sin(grid_theta))
         oridata=-np.multiply(bxdata,B3_data)
      elif q=='Reynolds':
         vxdata=np.multiply(v2_data,np.cos(grid_theta))+np.multiply(v1_data,np.sin(grid_theta))
         oridata=np.multiply(vxdata,np.multiply(v3_data,rho_data))
      elif q=='Radacc':
         oridata=np.divide(np.multiply(Fr01_data,(sigma_s_data+sigma_a_data)),rho_data)
      elif q=='RhoVr':
         oridata=np.multiply(v1_data,rho_data)
      elif q=='RhoVtheta':
         oridata=np.multiply(v2_data,rho_data)
      elif q=='RhoVphi':
         oridata=np.multiply(v3_data,rho_data)
      elif q=='RhoVout':
         oridata=np.multiply(vout_data,rho_data)
      elif q=='RhoVin':
         oridata=np.multiply(vin_data,rho_data)
      elif q=='Ekin1':
         oridata=0.5*np.multiply(v1_data, v1_data)
         oridata=np.multiply(oridata,rho_data)
      elif q=='Ekin2':
         oridata=0.5*np.multiply(v2_data, v2_data)
         oridata=np.multiply(oridata,rho_data)
      elif q=='Ekin3':
         oridata=0.5*np.multiply(v3_data, v3_data)
         oridata=np.multiply(oridata,rho_data)
      elif q=='BrBphi':
         oridata=np.multiply(B1_data,B3_data)
      elif q=='BthetaBphi':
         oridata=np.multiply(B2_data,B3_data)
      elif q=='rhovrvphi':
         oridata=np.multiply(v1_data,rhovphi_data)
      elif q=='rhovthetavphi':
         oridata=np.multiply(v2_data,rhovphi_data)
      elif q=='vrEr':
         oridata=np.multiply(v1_data,Er_data)
      elif q=='vthetaEr':
         oridata=np.multiply(v2_data,Er_data)
      elif q=='rhoPB':
         oridata=np.multiply(rho_data,PB_data)
      elif q=='tgas':
         oridata=np.divide(pgas_data,rho_data)
      elif q=='Fr01Sigma':
         oridata=np.multiply(Fr01_data,(sigma_s_data+sigma_a_data))
      elif q=='Fr02Sigma':
         oridata=np.multiply(Fr02_data,(sigma_s_data+sigma_a_data))
      elif q=='Fr01kappa':
         oridata=np.divide(np.multiply(Fr01_data,(sigma_s_data+sigma_a_data)),rho_data)
      elif q=='Fr02kappa':
         oridata=np.divide(np.multiply(Fr02_data,(sigma_s_data+sigma_a_data)),rho_data)
      else:
         oridata=f[dataset][index,block_num,:]
          
      for ko in ko_vals:
        for jo in jo_vals:
          for io in io_vals:
            data[q][kl+ko:ku+ko:s,jl+jo:ju+jo:s,il+io:iu+io:s] = oridata

    bar.update(block_num+1)


  bar.finish()

# get the cell center coordinates
  # get cell center coordinates
  data['x1v']=np.zeros(nx1)
  data['x2v']=np.zeros(nx2)
  data['x3v']=np.zeros(nx3)
  
  for ni in range(nx1):
    data['x1v'][ni]=0.75*(data['x1f'][ni+1]**4.0 - data['x1f'][ni]**4.0)/(data['x1f'][ni+1]**3.0-data['x1f'][ni]**3.0);
  
  for nj in range(nx2):
    data['x2v'][nj]=((np.sin(data['x2f'][nj+1]) - data['x2f'][nj+1] * np.cos(data['x2f'][nj+1])) \
           -(np.sin(data['x2f'][nj]) - data['x2f'][nj] * np.cos(data['x2f'][nj]))) \
              / (np.cos(data['x2f'][nj]) - np.cos(data['x2f'][nj+1]));

  for nk in range(nx3):
    data['x3v'][nk]=0.5*(data['x3f'][nk+1] + data['x3f'][nk]);



  #add time
  data['Time']=f.attrs['Time']

  f.close()

  return data





def MakeRhoVSlice(data_cart, vx_cart, vy_cart, minval, maxval, vlim1, vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname,vel=0,logscale=1, ledge=0.15, redge=0.78,title=None):
    plots, axes = plt.subplots(figsize=(10,8.4),dpi=300)
    plt.xlabel('$x/r_s$', size = 30)
    plt.ylabel('$y/r_s$', size = 30)
    plt.subplots_adjust(left=ledge,right=redge,top=0.85,bottom=0.1)
    if title is not None:
      plt.title(title,size=25,y=1.05)
    
    if vel>0:
      speed=np.sqrt(vx_cart**2+vy_cart**2)
      if vlim2 < vlim1:
         vlim2=speed.max()
      speed=np.clip(speed,vlim1,vlim2)
      logspeed=np.log10(speed)
      vcolor=axes.imshow(speed,cmap='RdYlBu_r',norm=LogNorm(vmin=vlim1,vmax=vlim2),origin='lower', extent=[xmin,xmax,ymin,ymax])
      velaxes = plots.add_axes([0.24, 0.89, 0.5, 0.03])
      velcbar = plots.colorbar(vcolor, cax = velaxes, orientation = 'horizontal')
      velcbar.set_label(label2, size=30, labelpad=-90)
      velcbar.ax.tick_params(labelsize=25)

      velfield=axes.streamplot(xgrid,ygrid,vx_cart,vy_cart,cmap='RdYlBu_r',density=1,color=logspeed)
        
      axes.set_ylim([ymin,ymax])

    if logscale > 0:
      im = axes.imshow(data_cart,cmap='RdYlBu_r', norm=LogNorm(vmin = minval, vmax=maxval), \
                       origin='lower', extent=[xmin,xmax,ymin,ymax])
    else:
      im = axes.imshow(data_cart,cmap='RdYlBu_r', vmin=minval, vmax=maxval, \
                       origin='lower', extent=[xmin,xmax,ymin,ymax])
    
    
    
    cbaxes = plots.add_axes([0.8,0.15,0.03,0.6])
    cbar=plots.colorbar(im,cax=cbaxes)
    cbar.set_label(label1, size=30)
    cbar.ax.tick_params(labelsize=25)
#    axes.set_xticks([xmin, width/4.0, width/2.0, width*3/4, width])
    axes.yaxis.set_tick_params(labelsize=20)
    axes.xaxis.set_tick_params(labelsize=20)
    
    axes.set_aspect('auto')
    plt.savefig(outputname)
    plt.close(plots)



######################################################################################################

crat=5694.76

num=11574
filename='disk.out1.'+'{:05d}'.format(num)+'.athdf'
data=ReduceData(filename)
quantities=data.keys()
meanrhovphi=np.mean(data['RhoVphi'],axis=0)
meanrho=np.mean(data['rho'],axis=0)
meanrhovr=np.mean(data['RhoVr'],axis=0)
meanvphi=meanrhovphi/meanrho


x2v=data['x2v']
x1v=data['x1v']
x3v=data['x3v']

ntheta=np.size(x2v)
nr=np.size(x1v)
nphi=np.size(x3v)

rho_1D=np.zeros(nr*nphi)
drhovrvphi_1D=np.zeros(nr*nphi)
max_1D=np.zeros(nr*nphi)
PB_1D=np.zeros(nr*nphi)
Er_1D=np.zeros(nr*nphi)

radius=np.zeros(nr*nphi)
phi=np.zeros(nr*nphi)


nloc=ntheta/2

for k in range(nphi):
  for i in range(nr):
    rho_1D[k*nr+i]=data['rho'][k,nloc,i]
    drhovrvphi_1D[k*nr+i]=data['rhovrvphi'][k,nloc,i]-meanrhovr[nloc,i]*meanvphi[nloc,i]
    max_1D[k*nr+i]=-data['B1'][k,nloc,i]*data['B3'][k,nloc,i]
    PB_1D[k*nr+i]=data['PB'][k,nloc,i]
    Er_1D[k*nr+i]=data['Er'][k,nloc,i]
    radius[k*nr+i]=x1v[i]
    phi[k*nr+i]=x3v[k]

# convert to cartesian coordinate

xcoord=radius*np.cos(phi)
ycoord=radius*np.sin(phi)

width=80
height=80
nx=256
ny=int(height*nx/width)

xmin=-width/2
# xmax=np.max(rcoord)
xmax=width/2
ymin=-height/2
ymax=height/2
xgrid=np.linspace(xmin, xmax, nx)
ygrid=np.linspace(ymin, ymax, ny)

xmesh,ymesh=np.meshgrid(xgrid,ygrid)

rho_cart=griddata(np.c_[xcoord,ycoord],rho_1D,(xmesh,ymesh),method='nearest')
reyn_cart=griddata(np.c_[xcoord,ycoord],drhovrvphi_1D,(xmesh,ymesh),method='nearest')
max_cart=griddata(np.c_[xcoord,ycoord],max_1D,(xmesh,ymesh),method='nearest')
PB_cart=griddata(np.c_[xcoord,ycoord],PB_1D,(xmesh,ymesh),method='nearest')
Er_cart=griddata(np.c_[xcoord,ycoord],Er_1D,(xmesh,ymesh),method='nearest')

title='${\\rm time}='+"%4.2f"%(data['Time']*crat)+'{\ r_s/c}$'

outputname='rho.png'

label1='$\\rho/\\rho_0$'
label2=''


MakeRhoVSlice(rho_cart, None, None, 1.e-3, 0.1, 1.e-2, 1.e2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname,0,logscale=1,title=title)




outputname='maxwell.png'

label1='$-B_xB_{\\phi}$'
label2=''

indices = max_cart < 1.e-10

max_cart[indices] = 1.e-10


MakeRhoVSlice(max_cart, None, None, 1.e-3, 5, 1.e-2, 1.e2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname,0,logscale=1,title=title)


outputname='reynolds.png'

label1='$\\rho v_r v_{\\phi} - <\\rho v_r><v_{\\phi}>$'
label2=''

indices = reyn_cart < 1.e-10

reyn_cart[indices] = 1.e-10


MakeRhoVSlice(reyn_cart, None, None, 1.0, 3000, 1.e-2, 1.e2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname,0,logscale=1,title=title)



outputname='PB.png'

label1='$P_B/P_0$'
label2=''


MakeRhoVSlice(PB_cart, None, None, 1.e-3, 100, 1.e-2, 1.e2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname,0,logscale=1,title=title)

outputname='Er.png'

label1='$E_r/a_rT_0^4$'
label2=''


MakeRhoVSlice(Er_cart, None, None, 10.0, 2000, 1.e-2, 1.e2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname,0,logscale=1,title=title)

