# Código realizado por JPI Ingeniería e Innovación SAC
# Realizado para el Curso:
# Análisis Sísmico de Edificaciones en Python

import openseespy.opensees as ops
import openseespy.postprocessing.ops_vis as opsv
import openseespy.postprocessing.Get_Rendering as opsplt
import matplotlib.pyplot as plt
from OSP_tools import *
import numpy as np
import pandas as pd

# Creamos el modelo
ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)

# Generamos la malla
RigidDiaphragm = 'ON'
dx, dy, alt  = 4, 4, 3 # m
nx, ny, nz = 3, 2, 3
if nx*dx>ny*dy:
    limx1 = 0
    limx2 = nx*dx
    limy1 = (ny*dy-limx2)/2
    limy2 = limx2 + limy1
else:
    limy1 = 0
    limy2 = ny*dy
    limx1 = (nx*dx-limy2)/2
    limx2 = limy2 + limx1

Nodes, Elems, Diap = GeoModel(dx,dy,alt,nx,ny,nz)

# Asignamos los nodos
for Ni in Nodes:
    ops.node(int(Ni[0]), *Ni[1:])

# Definimos diafragma rígidos
if RigidDiaphragm == 'ON':
    dirDia = 3 # perpendicular al plano del diafragma
    for Nd in Diap:
        # Nd[2]=Nd[2]+dy*ny*0.05 # excentricidad para analisis en X
        ops.node(int(Nd[0]), *Nd[1:]*m)
        ops.fix(int(Nd[0]),*[0,0,1,1,1,0])
        NodesDi = []
        for Ni in Nodes:
            if Ni[3]==Nd[3]:
                NodesDi.append(int(Ni[0]))
        ops.rigidDiaphragm(dirDia,int(Nd[0]),*NodesDi)

# Restricciones
ops.fixZ(0.0, *[1,1,1,1,1,1], '-tol', 1e-10)

# Definimos la sección
ops.section('Elastic', 1, E, Ac, Izc, Iyc, G, Jxxc)
ops.section('Elastic', 2, E, A, Iz, Iy, G, Jxx)
Nint = 5 #Number of integration points
ops.beamIntegration('Lobatto', 1, 1, Nint)
ops.beamIntegration('Lobatto', 2, 2, Nint)
CTransf = 1
BTransf = 2
ops.geomTransf('PDelta', CTransf, 1, 0, 0)
ops.geomTransf('Linear', BTransf, 0, 0, 1)

# Creamos los elementos
for Ele in Elems:
    if int(Ele[3]) == 1:# Column
        ops.element('forceBeamColumn', int(Ele[0]), int(Ele[1]), int(Ele[2]), int(Ele[3]), 1, '-mass', ρlc)
    else: # Beam
        ops.element('forceBeamColumn', int(Ele[0]), int(Ele[1]), int(Ele[2]), int(Ele[3]), 2, '-mass', ρl*(dx-a)/dx)

opsv.plot_model(fig_wi_he=(30., 25.),az_el=(-50,10),fig_lbrt = (0, 0, 1, 1))
plt.axis([limx1,limx2,limy1,limy2])
plt.savefig('./imagenes/Modelo_numerico.png',dpi=200)

# Aplicando Cargas vivas y muertas
wLive = 250*kg/m**2
wLosa = 300*kg/m**2
wAcab = 100*kg/m**2
wTabi = 150*kg/m**2
wTotal = 1.0*(wLosa+wAcab+wTabi)+0.25*wLive
#
Carga = wTotal*dx*dy*m**2
for Ni in Nodes:
    ops.mass(int(Ni[0]),Ni[4]*Carga,Ni[4]*Carga,0.0)
# print('Masa en el Nodo 36:',ops.nodeMass(36))

# ploteamos el modelo
ele_shapes = {}
#
for i in range(int(Ele[0])+1):
    if int(Elems[i,3])==1: # Column
        ele_shapes[i] = ['rect', [a, a]]
    else: # Beam
        ele_shapes[i] = ['rect', [b, h]]

opsv.plot_extruded_shapes_3d(ele_shapes, fig_wi_he=(30.0, 25.0), 
                        az_el=(-50,10),fig_lbrt = (0, 0, 1, 1))
plt.axis([limx1,limx2,limy1,limy2])
# plt.show()
plt.savefig('./imagenes/Modelo_3D.png',dpi=200)

# Obtenemos los modos
Nmodes = 12
Tmodes = ops.eigen(Nmodes)

# ploteamos el modo 1
opsplt.plot_modeshape(1, 500)
plt.savefig('./imagenes/Modo_1.png',dpi=200) # comentar plt.show() de código fuente
opsplt.plot_modeshape(2, 500)
plt.savefig('./imagenes/Modo_2.png',dpi=200) 
opsplt.plot_modeshape(3, 500)
plt.savefig('./imagenes/Modo_3.png',dpi=200)

Tmodes = ops.eigen(Nmodes)
for i in range(Nmodes):
    Tmodes[i] = 2*np.pi/Tmodes[i]**0.5

# guardando las formas modales
NT = 3*nz
Tags = ops.getNodeTags()
modo = np.zeros((Nmodes,NT))
for j in range(1,Nmodes+1):
    ind = 0
    for i in Tags[-nz:]:
        temp = ops.nodeEigenvector(i,j)
        modo[j-1,[ind,ind+1,ind+2]] = temp[0],temp[1],temp[-1]
        ind = ind + 3

# Realizamos un análisis para obtener la matriz de Masas
ops.system('FullGeneral')
ops.numberer("Plain")
ops.constraints('Transformation') 
ops.algorithm('Linear')
ops.analysis('Transient')
ops.integrator('GimmeMCK',1.0,0.0,0.0)
ops.analyze(1,0.0) 

# Obtenemos la matriz de Masas
N = ops.systemSize()         # Número de Grados de Libertad
Mmatrix = ops.printA('-ret')
Mmatrix = np.array(Mmatrix)
Mmatrix.shape = (N,N)
MF = Mmatrix[N-NT:,N-NT:]

# Obtenemos las pseudo aceleraciones
Ro=8.
E030 = espectro_E030(Tmodes,Z=0.45,U=1.0,S=1.0,Tp=0.4,Tl=2.5,R=Ro)
#
x = np.arange(0,3.0,0.02)
y = espectro_E030(x,Z=0.45,U=1.0,S=1.0,Tp=0.4,Tl=2.5,R=Ro)
plt.close('all')
plt.plot(x,y,'k-')
plt.xlabel('Periodo (s)')
plt.ylabel('Sa (g)')
plt.axis([0,3.0,0,1.1*y.max()])
plt.savefig('./imagenes/Espectro_E030.png')
#
H = np.arange(1,nz+1)*alt
P = sum(MF[0::3,0::3]) # Masa por nivel
print('Masa por Nivel:\n',P)
#
DDx, ΔDx, VDx, DDy, ΔDy, VDy, df1, df2 = getCombo(E030,MF,modo,NT,Tmodes)
print('\nFACTOR DE MASAS PARTICIPATIVAS')
print(df1.round(4))
#
print('\nANÁLISIS DINÁMICO SIN ESCALAR')
print(df2.round(4))

# Obteniendo fuerzas estáticas
F, k = get_static_loads(E030[0],P,H,Tmodes[0])
texto1 = 'Las fuerzas distribuidas en altura del análisis estático se \
obtinen a partir del coeficiente basal (%.4f) y el exponente k \
relacionado al periodo fundamental (%.4f). '%(E030[0],k)
texto1 = texto1 + 'A continuación, se muestran las fuerzas \
y los desplazamientos obtenidos del análisis estático en la \
dirección X e Y.'

# Definición de caso de carga estático
tsTag = 1
ops.timeSeries('Linear', tsTag, '-factor', 1.0)
ops.pattern('Plain', 1, tsTag, '-fact', 1.0)

# Aplicando fuerzas estáticas en X
Le = ny*dy*0.05
for i in range(nz):
    ops.load(int(Diap[i][0]),F[i]*9.81,0.,0.,0.,0.,F[i]*9.81*Le)

# Ejecución del análisis estático
ops.wipeAnalysis()
ops.system('FullGeneral')
ops.numberer('Plain')
ops.constraints('Transformation') 
ops.integrator('LoadControl',1)
ops.algorithm('Linear')
ops.analysis('Static')
ops.analyze(1)
ops.loadConst('-time', 0.0)	# para reestablecer el tiempo a 0s
#
VS = np.cumsum(F[::-1])[::-1]
print('\nANÁLISIS ESTÁTICO EN X')
df3 = pd.DataFrame(columns=['Nivel','Vx(ton)','UxMax(cm)','UyMax(cm)','DriftX(‰)','DriftY(‰)'])
tempX, tempY = 0., 0.
for i in range(nz):
    desX = ops.nodeDisp(int(Diap[i][0]))[0]
    desY = ops.nodeDisp(int(Diap[i][0]))[1]
    rotZ = ops.nodeDisp(int(Diap[i][0]))[5]
    desX = (desX + abs(rotZ*ny*dy/2))*0.75*Ro
    desY = (desY + abs(rotZ*nx*dx/2))*0.75*Ro
    driftX=1000.*(desX-tempX)/alt
    driftY=1000.*(desY-tempY)/alt 
    tempX, tempY = desX, desY
    df3 = df3.append({'Nivel':i+1,'Vx(ton)':VS[i]/1000,'UxMax(cm)':desX*100,'UyMax(cm)':desY*100,
                    'DriftX(‰)':driftX,'DriftY(‰)':driftY}, ignore_index=True)
print(df3.round(4))

# Definición de caso de carga estático
ops.remove('timeSeries',1)
ops.remove('loadPattern',1)
tsTag = 1
ops.timeSeries('Linear', tsTag, '-factor', 1.0)
ops.pattern('Plain', 1, tsTag, '-fact', 1.0)

# Aplicando fuerzas estáticas en Y
Le = nx*dx*0.05
for i in range(nz):
    ops.load(int(Diap[i][0]),0.,F[i]*9.81,0.,0.,0.,F[i]*9.81*Le)

# Ejecución del análisis estático
ops.wipeAnalysis()
ops.system('FullGeneral')
ops.numberer('Plain')
ops.constraints('Transformation') 
ops.integrator('LoadControl',1)
ops.algorithm('Linear')
ops.analysis('Static')
ops.analyze(1)
#
VS = np.cumsum(F[::-1])[::-1]
print('\nANÁLISIS ESTÁTICO EN Y')
df4 = pd.DataFrame(columns=['Nivel','Vy(ton)','UxMax(cm)','UyMax(cm)','DriftX(‰)','DriftY(‰)'])
tempX, tempY = 0., 0.
for i in range(nz):
    desX = ops.nodeDisp(int(Diap[i][0]))[0]
    desY = ops.nodeDisp(int(Diap[i][0]))[1]
    rotZ = ops.nodeDisp(int(Diap[i][0]))[5]
    desX = (desX + abs(rotZ*ny*dy/2))*0.75*Ro
    desY = (desY + abs(rotZ*nx*dx/2))*0.75*Ro
    driftX=1000.*(desX-tempX)/alt
    driftY=1000.*(desY-tempY)/alt 
    tempX, tempY = desX, desY
    df4 = df4.append({'Nivel':i+1,'Vy(ton)':VS[i]/1000,'UxMax(cm)':desX*100,'UyMax(cm)':desY*100,
                    'DriftX(‰)':driftX,'DriftY(‰)':driftY}, ignore_index=True)
print(df4.round(4))

# Escalamiento de los resultados del análisis dinámico
if VDx[0::3][0]<0.80*VS[0]:
    FSx  = 0.90*VS[0]/VDx[0::3][0]
    msjx = 'SI es necesario aplicar un Factor de Escala en X: %.4f'%FSx
else:
    FSx = 1.
    msjx = 'NO es necesario escalar en X'

if VDy[0::3][0]<0.80*VS[0]:
    FSy  = 0.90*VS[0]/VDy[0::3][0]
    msjy = 'SI es necesario aplicar un Factor de Escala en Y: %.4f'%FSy
else:
    FSy = 1.
    msjy = 'NO es necesario escalar en Y'

texto2 = '\nAl comparar la cortante basal obtenida en el análisis dinámico en X \
(%.2f ton) y el 80%% de la cortante basal del análisis estático en X (%.2f ton), \
se obtiene que %s. '%(VDx[0::3][0]/1000,0.80*VS[0]/1000,msjx)
texto2 = texto2 + 'En la dirección Y, la cortante basal obtenida en el análisis \
dinámico es %.2f ton y el 80%% de la cortante basal del análisis estático es %.2f ton. \
Por lo que %s.'%(VDy[0::3][0]/1000,0.80*VS[0]/1000,msjy)

# Se aplican los Factores de Escala
print('\nANÁLISIS DINÁMICO FINAL')
df5 = pd.DataFrame(columns=['Nivel','Vx(ton)','Vy(ton)','Ux(cm)','Uy(cm)','Δx(‰)','Δy(‰)'])
for i in range(nz):
    Δx = FSx*0.75*Ro*ΔDx[0::3][i]/alt
    Δy = FSy*0.75*Ro*ΔDy[1::3][i]/alt
    #
    df5 = df5.append({'Nivel':i+1, 'Vx(ton)':FSx*VDx[0::3][i]/1000,
        'Vy(ton)':FSy*VDy[1::3][i]/1000,'Ux(cm)':FSx*0.75*Ro*DDx[0::3][i]*100,
        'Uy(cm)':FSy*0.75*Ro*DDy[1::3][i]*100,'Δx(‰)':Δx*1000,'Δy(‰)':Δy*1000}, ignore_index=True)
df5 = df5.astype({'Nivel':int})
print(df5.round(4))
#
vecX = np.array(df5.loc[:,'Δx(‰)'])
vecY = np.array(df5.loc[:,'Δy(‰)'])
lim = 1.1*max(vecX.max(),vecY.max())
plt.close('all')
plt.plot(np.insert(vecX,0,0),np.arange(nz+1),'bo--',label='drift X',lw = 0.8)
plt.plot(np.insert(vecY,0,0),np.arange(nz+1),'ro--',label='drift Y',lw = 0.8)
plt.legend()
plt.xlabel('Distorsión (‰)')
plt.ylabel('Nivel')
plt.axis([-0.05,lim,-0.05,nz+0.05])
plt.yticks(np.arange(0, nz+0.05, 1))
plt.savefig('./imagenes/distorsion_din.png')

genReport(df1,df2.iloc[::-1],df3.iloc[::-1],df4.iloc[::-1],df5.iloc[::-1],texto1,texto2)
