import numpy as np
import cartopy as cp
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import netCDF4 as nc
import matplotlib.pyplot as plt
from cartopy import config
import cartopy.feature as cfeature
from matplotlib.axes import Axes
import matplotlib.gridspec as gridspec
from cartopy.mpl.geoaxes import GeoAxes
import cartopy.crs as ccrs
import cartopy.geodesic as cgeo
import os

directory = "C:/Users/Luisv/Desktop/Situaciones2"
archivos = []
for nombre in os.listdir(directory):
    archivos = archivos +[os.path.join(directory, nombre)]
leer = nc.Dataset(archivos[5])
time = leer.dimensions['Time']
lat = leer.dimensions['south_north']
lon = leer.dimensions['west_east']
lats = leer ['XLAT'][0,:,:]
lons = leer ['XLONG'][0,:,:]
time_var = leer.variables['Times'][:]
temp_veces = time.size - 1
temp_max = np.zeros([lat.size,lon.size])
temp_min = np.zeros([lat.size,lon.size])
vientos_max = np.zeros([lat.size,lon.size])
prec_max = np.zeros([lat.size,lon.size])
sumador = np.zeros([lat.size,lon.size])
text_kwargs = dict(family='serif', size=40, color='black')
momentotempfrio = 0
mediamax = 273
mediamin = 273
vientos = 0
i = 0
u1 = 0
v1 = 12
h = 5
print(time.size)
prec_max = np.zeros([lat.size,lon.size])
sumador = np.zeros([lat.size,lon.size])
i = time.size-1
while i >= 12:
    prec_suma = np.zeros([lat.size,lon.size])
    prec_resta = np.zeros([lat.size,lon.size])
    
    prec_c = leer ['RAINC'][i,:,:]
    prec_nc = leer ['RAINNC'][i,:,:]
    prec_suma = prec_c + prec_nc
    prec_c = leer ['RAINC'][i-12,:,:]
    prec_nc = leer ['RAINNC'][i-12,:,:]
    prec_resta = prec_c + prec_nc
    
    prec_12 = prec_suma - prec_resta
    if prec_12.mean() >= prec_max.mean():
        prec_max = prec_12
        prec_12 = np.zeros([lat.size,lon.size])
        print("Me disparo en" ,i)
        inicio_pre = i-12
        final_pre = i
    i = i-1

    
for i in range (0, temp_veces):
    temp = leer ['T2'][i,:,:]
    vientoU = leer ['V10'][i,:,:]
    vientoV = leer ['U10'][i,:,:]
    media_t = np.ma.max(temp)
    media_frio = np.ma.min(temp)
    media_v = np.ma.average(np.sqrt(vientoU*vientoU + vientoV*vientoV))
    if media_t > mediamax:
        mediamax = media_t
        temp_max = temp
        momentotempcalor = i
        print ("TEMPERATURA MAX", i)
    elif media_frio < mediamin:
        mediamin = media_frio
        temp_min = temp
        print ("TEMPERATURA MIN", i)
        momentotempfrio = i
    elif media_v > vientos:
        vientos = media_v
        vientos_maxU = vientoU
        vientos_maxV = vientoV
        print ("TVIENTO", i)
    elif mediamin == 273:
        temp_min = temp
    i = i + 1
leer.close()

levels_max = [302, 305, 307, np.ma.max(temp_max)+30]
levels_min = [240,265, 269, 272]
levels_lluvia = [40, 80, 120, np.ma.max(prec_max)+300]
levels_viento = [10, 20, 30]
pluv = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
avisos=['                            Bajo',
        '                            Medio',
        '                            Alto','']
avisosr=['                            Alto',
        '                            Medio',
        '                            Bajo','']

colores = [(255/255, 255/255, 1/255),
           (255/255, 153/255, 51/255),
           (190/255, 1/255, 1/255),
           (0/255, 0/255, 0/255)]

colores2 = [(0/255, 0/255, 0/255),
            (190/255, 1/255, 1/255),
           (255/255, 153/255, 51/255),
           (255/255, 255/255, 1/255)]



fechacal = time_var[momentotempcalor].tobytes
fechafrio = time_var[momentotempfrio].tobytes
fechapre_in = time_var[inicio_pre].tobytes
fechapre_fin = time_var[final_pre].tobytes

titulocalor= fechacal("utf-8")
titulofrio = fechafrio("utf-8")
titulo_pre_inicio = fechapre_in("utf-8")
titulo_pre_final = fechapre_fin("utf-8")


my_cmap = LinearSegmentedColormap.from_list('basico', colores)
my_cmap_r = LinearSegmentedColormap.from_list('basico', colores2)



#mapa cal#
jjj = plt.figure(figsize=(40, 20))
ax1 = plt.subplot(2, 1, 1, projection=ccrs.PlateCarree())
ax1.coastlines()
ax1.set_extent([-10, 4, 41, 44], crs=ccrs.PlateCarree())
tt = ax1.contourf(lons, lats, temp_max, 40, levels=levels_max, cmap = my_cmap)
cbar = plt.colorbar(tt, orientation="horizontal", extendrect=True, ticks=levels_max)
cbar.ax.set_title('Nivel de alerta', fontsize = 40)
cbar.ax.set_xticklabels(avisos)
cbar.ax.tick_params(labelsize=50)
gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
gl.xlabel_style = {'size': 40, 'color': 'gray'}
gl.ylabel_style = {'size': 40, 'color': 'gray'}
gl.xlabels_bottom = False
gl.ylabels_right = False
ax1.add_feature(cp.feature.OCEAN, zorder=100, edgecolor='k')
ax1.add_feature(cp.feature.LAND)
plt.title('Avisos', fontsize=50, pad=20)

ax2 = plt.subplot(2, 1, 2, projection=ccrs.PlateCarree())
ax2.coastlines()
ax2.set_extent([-10, 4, 41, 44], crs=ccrs.PlateCarree())
temp_max= temp_max -273
levels_max = 0
tt = ax2.contourf(lons, lats, temp_max, 40, cmap = my_cmap)
cbar = plt.colorbar(tt, orientation="horizontal", extendrect=True)
cbar.ax.set_title('Temperatura (ºC)', fontsize = 40)
cbar.ax.tick_params(labelsize=50)
gl = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
gl.xlabel_style = {'size': 40, 'color': 'gray'}
gl.ylabel_style = {'size': 40, 'color': 'gray'}
gl.xlabels_bottom = False
gl.ylabels_right = False
ax2.add_feature(cp.feature.OCEAN, zorder=100, edgecolor='k')
ax2.add_feature(cp.feature.LAND)
plt.title('Temperatura', fontsize=50, pad=20,y=1.2)
plt.suptitle('Altas temperaturas para\n' +' '+titulocalor.decode(), fontsize=60, y=1.02)
jjj=plt.tight_layout(pad=3.0)

#mapa frio#
jjj = plt.figure(figsize=(40, 20))
ax1 = plt.subplot(2, 1, 1, projection=ccrs.PlateCarree())
ax1.coastlines()
ax1.set_extent([-10, 4, 41, 44], crs=ccrs.PlateCarree())
tt = ax1.contourf(lons, lats, temp_min, 40, levels=levels_min, cmap = my_cmap_r)
cbar = plt.colorbar(tt, orientation="horizontal", extendrect=True, ticks=levels_min)
cbar.ax.set_title('Nivel de alerta', fontsize = 40)
cbar.ax.set_xticklabels(avisosr)
cbar.ax.tick_params(labelsize=50)
gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
gl.xlabel_style = {'size': 40, 'color': 'gray'}
gl.ylabel_style = {'size': 40, 'color': 'gray'}
gl.xlabels_bottom = False
gl.ylabels_right = False
ax1.add_feature(cp.feature.OCEAN, zorder=100, edgecolor='k')
ax1.add_feature(cp.feature.LAND)
plt.title('Avisos', fontsize=50, pad=20)

ax2 = plt.subplot(2, 1, 2, projection=ccrs.PlateCarree())
ax2.coastlines()
ax2.set_extent([-10, 4, 41, 44], crs=ccrs.PlateCarree())
temp_min = temp_min-273
levels_min=0
tt = ax2.contourf(lons, lats, temp_min, 40, cmap = 'BuPu_r' ,extend='min')
cbar = plt.colorbar(tt, orientation="horizontal", extendrect=True)
cbar.ax.set_title('Temperatura (ºC)', fontsize = 40)
cbar.ax.tick_params(labelsize=50)
gl= ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
gl.xlabel_style = {'size': 40, 'color': 'gray'}
gl.ylabel_style = {'size': 40, 'color': 'gray'}
gl.xlabels_bottom = False
gl.ylabels_right = False
ax2.add_feature(cp.feature.OCEAN, zorder=100, edgecolor='k')
plt.title('Temperatura', fontsize=50, pad=20, y=1.2)
plt.suptitle('Bajas temperaturas para\n' + ' '+ titulofrio.decode() , fontsize=60, y=1.02)
jjj=plt.tight_layout(pad=3.0)

#mapa prec#
jjj= plt.figure(figsize=(40, 20))
plt.subplots_adjust(wspace=0.5,bottom=0.1)
ax1 = plt.subplot(2, 1, 1, projection=ccrs.PlateCarree())
ax1.coastlines()
ax1.set_extent([-10, 4, 41, 44], crs=ccrs.PlateCarree())
tt = ax1.contourf(lons, lats, prec_max, 40, levels=levels_lluvia, cmap = my_cmap)
cbar = plt.colorbar(tt, orientation="horizontal", extendrect=True)
cbar.ax.set_title('Nivel de alerta', fontsize = 40)
cbar.ax.set_xticklabels(avisos)
cbar.ax.tick_params(labelsize=50)

gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
gl.xlabel_style = {'size': 40, 'color': 'gray'}
gl.ylabel_style = {'size': 40, 'color': 'gray'}
gl.xlabels_bottom = False
gl.ylabels_right = False
ax1.add_feature(cp.feature.OCEAN, edgecolor='k')
ax1.add_feature(cp.feature.LAND)
plt.title('Avisos', fontsize=50, pad=20)



ax2 = plt.subplot(2, 1, 2, projection=ccrs.PlateCarree())
ax2.coastlines()
ax2.set_extent([-10, 4, 41, 44], crs=ccrs.PlateCarree())
tt = ax2.contourf(lons, lats, prec_max, 40, levels=pluv, cmap = 'Blues')
cbar = plt.colorbar(tt, orientation="horizontal", extendrect=True)
cbar.ax.set_title('Precipitación acumulada (mm)', fontsize = 40)
cbar.ax.tick_params(labelsize=50)
gl= ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
gl.xlabel_style = {'size': 40, 'color': 'gray'}
gl.ylabel_style = {'size': 40, 'color': 'gray'}
gl.xlabels_bottom = False
gl.ylabels_right = False
ax2.add_feature(cp.feature.OCEAN, edgecolor='k')
ax2.add_feature(cp.feature.LAND)
plt.title('Precipitación', fontsize=50, pad=20, y=1.2)
plt.suptitle('Lluvias entre' + ' ' + titulo_pre_inicio.decode()+' '+ 'y\n' +' '+titulo_pre_final.decode(), fontsize=60, y=1.02)
jjj=plt.tight_layout(pad=3.0)


