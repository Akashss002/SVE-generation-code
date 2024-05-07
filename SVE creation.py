from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import regionToolset
import random
import numpy as np
import csv
import time


start_time = time.time()
#Parameters (Length in mm)
a = 0.005 #gap from bottom
b = 0.005 #gap from left
g = 0.008 #gap between fiber centres
r = 0.0035
tol = 0.1 #10% tolerance
r1 = (1-tol)*r
r2 = (1+tol)*r
n = 10 #no. of fibers in a row
m = 10 #no. of fibers in a column
sp = 0.01 #Strain percentage
fiber_prop = np.array([241e+03,0.38])
matrix_prop = np.array([3.12e+03, 0.2])
num = 50
mesh_s = 0.0035

def parts_creation(modelname,r_n):
    mdb.Model(modelType=STANDARD_EXPLICIT, name=modelname)
    mdb.models[modelname].ConstrainedSketch(name='__profile__', sheetSize=1.0)
    l = 0
    for i in range(0,n):
        for j in range(0,m):
            k = r_n[l]
            mdb.models[modelname].sketches['__profile__'].CircleByCenterPerimeter(center=(
                a+i*g, b+j*g), point1=(a+i*g+k, b+j*g))
            l+=1
    mdb.models[modelname].Part(dimensionality=THREE_D, name='Fibers', type=
        DEFORMABLE_BODY)
    mdb.models[modelname].parts['Fibers'].BaseSolidExtrude(depth=2*a+(n-1)*g, sketch=
        mdb.models[modelname].sketches['__profile__'])

    mdb.models[modelname].ConstrainedSketch(name='__profile__2', sheetSize=1.0)
    mdb.models[modelname].sketches['__profile__2'].rectangle(point1=(0.0, 0.0), 
        point2=(2*a+(n-1)*g, 2*b + (m-1)*g))
    mdb.models[modelname].Part(dimensionality=THREE_D, name='Matrix', type=
        DEFORMABLE_BODY)
    mdb.models[modelname].parts['Matrix'].BaseSolidExtrude(depth=2*a+(n-1)*g, sketch=
        mdb.models[modelname].sketches['__profile__2'])
    # return arr

def assembly_creation(modelname):
    mdb.models[modelname].rootAssembly.DatumCsysByDefault(CARTESIAN)
    mdb.models[modelname].rootAssembly.Instance(dependent=ON, name='Fibers-1', 
        part=mdb.models[modelname].parts['Fibers'])
    mdb.models[modelname].rootAssembly.Instance(dependent=ON, name='Matrix-1', 
        part=mdb.models[modelname].parts['Matrix'])
    mdb.models[modelname].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
        instances=(mdb.models[modelname].rootAssembly.instances['Matrix-1'], 
        mdb.models[modelname].rootAssembly.instances['Fibers-1']), 
        keepIntersections=ON, name='RVE', originalInstances=SUPPRESS)
    del mdb.models[modelname].parts['Matrix']
    del mdb.models[modelname].parts['Fibers']

def material_creation(modelname):
    mdb.models[modelname].Material(name='Matrix')
    mdb.models[modelname].materials['Matrix'].Elastic(table=(matrix_prop, ))
    mdb.models[modelname].HomogeneousSolidSection(material='Matrix', name='MatrixSec', thickness=None)

    mdb.models[modelname].Material(name='Fiber')
    mdb.models[modelname].materials['Fiber'].Elastic(table=(fiber_prop, ))
    mdb.models[modelname].HomogeneousSolidSection(material='Fiber', name='FiberSec', thickness=None)

def material_assign(modelname):
    c = mdb.models[modelname].parts['RVE'].cells
    N = n*m
    cells_mat = c[N:N+1]
    region = mdb.models[modelname].parts['RVE'].Set(cells=cells_mat, name='Set-1')
    mdb.models[modelname].parts['RVE'].SectionAssignment(offset=0.0, offsetField=''
        , offsetType=MIDDLE_SURFACE, region=region, sectionName='MatrixSec', thicknessAssignment=FROM_SECTION)
    cells_fib = c[0:N]
    region = mdb.models[modelname].parts['RVE'].Set(cells=cells_fib, name='Set-2')
    mdb.models[modelname].parts['RVE'].SectionAssignment(offset=0.0, offsetField=''
        , offsetType=MIDDLE_SURFACE, region=region, sectionName='FiberSec', thicknessAssignment=FROM_SECTION)

def steps(modelname):
    mdb.models[modelname].StaticStep(name='Step-1', previous='Initial')
    mdb.models[modelname].fieldOutputRequests['F-Output-1'].setValues(variables=(
        'S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'CSTRESS', 'CDISP', 'IVOL'))
    # mdb.models[modelname].FieldOutputRequest(createStepName='Step-3', name=
    #     'F-Output-3', variables=('S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 
    #     'CSTRESS', 'CDISP', 'IVOL'))
    # mdb.models[modelname].fieldOutputRequests['F-Output-2'].deactivate('Step-3')
    # mdb.models[modelname].fieldOutputRequests['F-Output-3'].deactivate('Step-4')
    # mdb.models[modelname].FieldOutputRequest(createStepName='Step-4', name=
    #     'F-Output-4', variables=('S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 
    #     'CSTRESS', 'CDISP', 'IVOL'))
    # mdb.models[modelname].FieldOutputRequest(createStepName='Step-5', name=
    #     'F-Output-5', variables=('S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 
    #     'CSTRESS', 'CDISP', 'IVOL'))
    # mdb.models[modelname].fieldOutputRequests['F-Output-4'].deactivate('Step-5')

def facesets(modelname):
    f = mdb.models[modelname].rootAssembly.instances['RVE-1'].faces
    N = len(f)
    XFront = f[3*n*m+2:3*n*m+3]
    XBack = f[3*n*m:3*n*m+1]
    YTop = f[3*n*m+1:3*n*m+2]
    YBottom = f[3*n*m+3:3*n*m+4]
    mdb.models[modelname].rootAssembly.Set(faces=XFront, name='XFront')
    mdb.models[modelname].rootAssembly.Set(faces=YTop, name='YTop')
    mdb.models[modelname].rootAssembly.Set(faces=XBack, name='XBack')
    mdb.models[modelname].rootAssembly.Set(faces=YBottom, name='YBottom')
    ZFront1 = f[n*m:2*n*m]
    ZFront2 = f[N-2:N-1]
    mdb.models[modelname].rootAssembly.Set(faces=ZFront1+ZFront2, name='ZFront')
    ZBack1 = f[2*n*m:3*n*m]
    ZBack2 = f[N-1:N]
    mdb.models[modelname].rootAssembly.Set(faces=ZBack1+ZBack2, name='ZBack')

    mdb.models[modelname].rootAssembly.ReferencePoint(point=(0.0, 0.0, 1.5*(2*a+(n-1)*g)))
    RP_position1 = mdb.models[modelname].rootAssembly.referencePoints.findAt((0.0, 0.0, 1.5*(2*a+(n-1)*g)),)
    RP1 = (RP_position1,)
    mdb.models[modelname].rootAssembly.Set(referencePoints=RP1, name = 'RP_set1')

def loads_XTens(modelname):
    mdb.models[modelname].Equation(name='XConstraint', terms=((1.0, 'XFront', 1), (
    -1.0, 'RP_set1', 1)))
    #initial BCs
    mdb.models[modelname].ZsymmBC(createStepName='Initial', localCsys=None, name=
    'SymZBack', region=mdb.models[modelname].rootAssembly.sets['ZBack'])
    mdb.models[modelname].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name=
        'RollerYBottom', region=mdb.models[modelname].rootAssembly.sets['YBottom'], 
        u1=UNSET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    mdb.models[modelname].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name='RollerXBack', 
        region=mdb.models[modelname].rootAssembly.sets['XBack'], u1=SET, u2=UNSET, 
        u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    
    #step1 BC
    mdb.models[modelname].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
    'DisplaceX', region=mdb.models[modelname].rootAssembly.sets['RP_set1'], u1=
    sp*(2*a+(n-1)*g), u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

def loads_YTens(modelname):
    mdb.models[modelname].Equation(name='YConstraint', terms=((1.0, 'YTop', 2), (
    -1.0, 'RP_set1', 2)))
    #initial BCs
    mdb.models[modelname].ZsymmBC(createStepName='Initial', localCsys=None, name=
    'SymZBack', region=mdb.models[modelname].rootAssembly.sets['ZBack'])
    mdb.models[modelname].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name=
        'RollerYBottom', region=mdb.models[modelname].rootAssembly.sets['YBottom'], 
        u1=UNSET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    mdb.models[modelname].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name='RollerXBack', 
        region=mdb.models[modelname].rootAssembly.sets['XBack'], u1=SET, u2=UNSET, 
        u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    
    #Step1 BC
    mdb.models[modelname].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
    'DisplaceY', region=mdb.models[modelname].rootAssembly.sets['RP_set1'], u1=
    UNSET, u2=sp*(2*a+(n-1)*g), u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

def loads_ZTens(modelname):
    mdb.models[modelname].Equation(name='ZConstraint', terms=((1.0, 'ZFront', 3), (
    -1.0, 'RP_set1', 3)))
    #initial BCs
    mdb.models[modelname].ZsymmBC(createStepName='Initial', localCsys=None, name=
    'SymZBack', region=mdb.models[modelname].rootAssembly.sets['ZBack'])
    mdb.models[modelname].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name=
        'RollerYBottom', region=mdb.models[modelname].rootAssembly.sets['YBottom'], 
        u1=UNSET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    mdb.models[modelname].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name='RollerXBack', 
        region=mdb.models[modelname].rootAssembly.sets['XBack'], u1=SET, u2=UNSET, 
        u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

    #Step1 BC
    mdb.models[modelname].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
    'DisplaceZ', region=mdb.models[modelname].rootAssembly.sets['RP_set1'], u1=UNSET,
      u2=UNSET, u3=sp*(2*a+(n-1)*g), ur1=UNSET, ur2=UNSET, ur3=UNSET)

def loads_XYShear(modelname):
    mdb.models[modelname].Equation(name='XYConstraint', terms=((1.0, 'XFront', 2), (
    -1.0, 'RP_set1', 2)))
    #initial BCs
    # mdb.models[modelname].ZsymmBC(createStepName='Initial', localCsys=None, name=
    # 'SymZBack', region=mdb.models[modelname].rootAssembly.sets['ZBack'])
    mdb.models[modelname].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name=
        'RollerXFront', region=mdb.models[modelname].rootAssembly.sets['XFront'], 
        u1=SET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    mdb.models[modelname].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name='FixedXBack', 
        region=mdb.models[modelname].rootAssembly.sets['XBack'], u1=0, u2=0, 
        u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET)

    #Step1 BC
    mdb.models[modelname].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
    'DisplaceY2', region=mdb.models[modelname].rootAssembly.sets['RP_set1'], u1=UNSET,
      u2=sp*(2*a+(n-1)*g), u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    
def loads_YZShear(modelname):
    mdb.models[modelname].Equation(name='YZConstraint', terms=((1.0, 'YTop', 3), (
    -1.0, 'RP_set1', 3)))
    #initial BCs
    # mdb.models[modelname].ZsymmBC(createStepName='Initial', localCsys=None, name=
    # 'SymZBack', region=mdb.models[modelname].rootAssembly.sets['ZBack'])
    mdb.models[modelname].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name=
        'RollerYTop', region=mdb.models[modelname].rootAssembly.sets['YTop'], 
        u1=UNSET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    mdb.models[modelname].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name='FixedYBottom', 
        region=mdb.models[modelname].rootAssembly.sets['YBottom'], u1=0, u2=0, 
        u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET)

    #Step1 BC
    mdb.models[modelname].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
    'DisplaceZ2', region=mdb.models[modelname].rootAssembly.sets['RP_set1'], u1=UNSET,
      u2=UNSET, u3=sp*(2*a+(n-1)*g), ur1=UNSET, ur2=UNSET, ur3=UNSET)

def loads_ZXShear(modelname):
    mdb.models[modelname].Equation(name='ZXConstraint', terms=((1.0, 'ZFront', 1), (
    -1.0, 'RP_set1', 1)))
    #initial BCs
    # mdb.models[modelname].ZsymmBC(createStepName='Initial', localCsys=None, name=
    # 'SymZBack', region=mdb.models[modelname].rootAssembly.sets['ZBack'])
    mdb.models[modelname].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name=
        'RollerZFront', region=mdb.models[modelname].rootAssembly.sets['ZFront'], 
        u1=UNSET, u2=UNSET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    
    #Step1 BC
    mdb.models[modelname].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name='FixedZBack', 
        region=mdb.models[modelname].rootAssembly.sets['ZBack'], u1=0, u2=0, 
        u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    
    mdb.models[modelname].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
    'DisplaceX2', region=mdb.models[modelname].rootAssembly.sets['RP_set1'], u1=sp*(2*a+(n-1)*g),
      u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

def mesh(modelname):
    mdb.models[modelname].parts['RVE'].seedPart(deviationFactor=0.1, minSizeFactor=
    0.1, size=mesh_s)
    mdb.models[modelname].parts['RVE'].generateMesh()

def job(modelname):
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model=modelname, modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name=str(n*m)+'fibres_'+str(modelname), nodalOutputPrecision=SINGLE, 
    numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
    ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
    mdb.jobs[str(n*m)+'fibres_'+str(modelname)].submit(consistencyChecking=OFF)
    mdb.jobs[str(n*m)+'fibres_'+str(modelname)].waitForCompletion()

def results_E11(modelname):
    f = openOdb(path=str(n*m)+'fibres_'+str(modelname)+'.odb')
    Stress = f.steps['Step-1'].frames[-1].fieldOutputs['S'].values
    Strain = f.steps['Step-1'].frames[-1].fieldOutputs['E'].values
    Volume = f.steps['Step-1'].frames[-1].fieldOutputs['IVOL'].values

    V = 0
    for i in range(len(Volume)):
        V = V + Volume[i].data

    S33_eff = 0.0
    for i in range(len(Volume)):
        S33_eff = S33_eff + Stress[i].data[2]*Volume[i].data/V

    Str33_eff = 0.0
    Str11_eff = 0.0
    Str22_eff = 0.0
    for i in range(len(Volume)):
        Str33_eff = Str33_eff + Strain[i].data[2]*Volume[i].data/V
        Str11_eff = Str11_eff + Strain[i].data[0]*Volume[i].data/V
        Str22_eff = Str22_eff + Strain[i].data[1]*Volume[i].data/V
    E11 = S33_eff/Str33_eff
    v12 = -Str22_eff/Str33_eff
    v13 = -Str11_eff/Str33_eff
    return E11,v12,v13

def results_E22(modelname):
    f = openOdb(path=str(n*m)+'fibres_'+str(modelname)+'.odb')
    Stress = f.steps['Step-1'].frames[-1].fieldOutputs['S'].values
    Strain = f.steps['Step-1'].frames[-1].fieldOutputs['E'].values
    Volume = f.steps['Step-1'].frames[-1].fieldOutputs['IVOL'].values

    V = 0
    for i in range(len(Volume)):
        V = V + Volume[i].data

    S22_eff = 0.0
    for i in range(len(Volume)):
        S22_eff = S22_eff + Stress[i].data[1]*Volume[i].data/V
    
    Str22_eff = 0.0
    Str11_eff = 0.0
    Str33_eff = 0.0
    for i in range(len(Volume)):
        Str22_eff = Str22_eff + Strain[i].data[1]*Volume[i].data/V
        Str33_eff = Str33_eff + Strain[i].data[2]*Volume[i].data/V
        Str11_eff = Str11_eff + Strain[i].data[0]*Volume[i].data/V
    
    E22 = S22_eff/Str22_eff
    v23 = -Str11_eff/Str22_eff
    v21 = -Str33_eff/Str22_eff
    return E22,v21,v23

def results_E33(modelname):
    f = openOdb(path=str(n*m)+'fibres_'+str(modelname)+'.odb')
    Stress = f.steps['Step-1'].frames[-1].fieldOutputs['S'].values
    Strain = f.steps['Step-1'].frames[-1].fieldOutputs['E'].values
    Volume = f.steps['Step-1'].frames[-1].fieldOutputs['IVOL'].values

    V = 0
    for i in range(len(Volume)):
        V = V + Volume[i].data

    S11_eff = 0.0
    for i in range(len(Volume)):
        S11_eff = S11_eff + Stress[i].data[0]*Volume[i].data/V
    
    Str11_eff = 0.0
    Str22_eff = 0.0
    Str33_eff = 0.0
    
    for i in range(len(Volume)):
        Str11_eff = Str11_eff + Strain[i].data[0]*Volume[i].data/V
        Str22_eff = Str22_eff + Strain[i].data[3]*Volume[i].data/V
        Str33_eff = Str33_eff + Strain[i].data[4]*Volume[i].data/V
        
    E33 = S11_eff/Str11_eff
    v32 = Str22_eff/Str11_eff
    v31 = Str33_eff/Str11_eff
    
    return E33,v31,v32

def results_G23(modelname):
    f = openOdb(path=str(n*m)+'fibres_'+str(modelname)+'.odb')
    Stress = f.steps['Step-1'].frames[-1].fieldOutputs['S'].values
    Strain = f.steps['Step-1'].frames[-1].fieldOutputs['E'].values
    Volume = f.steps['Step-1'].frames[-1].fieldOutputs['IVOL'].values

    V = 0
    for i in range(len(Volume)):
        V = V + Volume[i].data

    S12_eff = 0.0
    for i in range(len(Volume)):
        S12_eff = S12_eff + Stress[i].data[3]*Volume[i].data/V
    
    Str12_eff = 0.0
    for i in range(len(Volume)):
        Str12_eff = Str12_eff + Strain[i].data[3]*Volume[i].data/V
    
    G23 = S12_eff/Str12_eff
    return G23

def results_G13(modelname):
    f = openOdb(path=str(n*m)+'fibres_'+str(modelname)+'.odb')
    Stress = f.steps['Step-1'].frames[-1].fieldOutputs['S'].values
    Strain = f.steps['Step-1'].frames[-1].fieldOutputs['E'].values
    Volume = f.steps['Step-1'].frames[-1].fieldOutputs['IVOL'].values

    V = 0
    for i in range(len(Volume)):
        V = V + Volume[i].data

    S13_eff = 0.0
    for i in range(len(Volume)):
        S13_eff = S13_eff + Stress[i].data[4]*Volume[i].data/V
    
    Str13_eff = 0.0
    for i in range(len(Volume)):
        Str13_eff = Str13_eff + Strain[i].data[4]*Volume[i].data/V
    
    G13 = S13_eff/Str13_eff
    return G13

def results_G12(modelname):
    f = openOdb(path=str(n*m)+'fibres_'+str(modelname)+'.odb')
    Stress = f.steps['Step-1'].frames[-1].fieldOutputs['S'].values
    Strain = f.steps['Step-1'].frames[-1].fieldOutputs['E'].values
    Volume = f.steps['Step-1'].frames[-1].fieldOutputs['IVOL'].values

    V = 0
    for i in range(len(Volume)):
        V = V + Volume[i].data

    S23_eff = 0.0
    for i in range(len(Volume)):
        S23_eff = S23_eff + Stress[i].data[5]*Volume[i].data/V
    
    Str23_eff = 0.0
    for i in range(len(Volume)):
        Str23_eff = Str23_eff + Strain[i].data[5]*Volume[i].data/V
    
    G12 = S23_eff/Str23_eff
    return G12

def XTens(r_n,num):
    modelname = 'XTens'+str(num)  
    parts_creation(modelname,r_n)
    assembly_creation(modelname)
    material_creation(modelname)
    material_assign(modelname)
    steps(modelname)
    facesets(modelname)
    loads_XTens(modelname)
    mesh(modelname)
    job(modelname)
    E33,v31,v32 = results_E33(modelname)
    return E33,v31,v32

def YTens(r_n,num):
    modelname = 'YTens'+str(num)   
    parts_creation(modelname,r_n)
    assembly_creation(modelname)
    material_creation(modelname)
    material_assign(modelname)
    steps(modelname)
    facesets(modelname)
    loads_YTens(modelname)
    mesh(modelname)
    job(modelname)
    E22,v21,v23 = results_E22(modelname)
    return E22,v21,v23

def ZTens(r_n,num):
    modelname = 'ZTens'+str(num) 
    parts_creation(modelname,r_n)
    assembly_creation(modelname)
    material_creation(modelname)
    material_assign(modelname)
    steps(modelname)
    facesets(modelname)
    loads_ZTens(modelname)
    mesh(modelname)
    job(modelname)
    E11,v12,v13 = results_E11(modelname)
    return E11,v12,v13

def XYShear(r_n,num):
    modelname = 'XYShear'+str(num)   
    parts_creation(modelname,r_n)
    assembly_creation(modelname)
    material_creation(modelname)
    material_assign(modelname)
    steps(modelname)
    facesets(modelname)
    loads_XYShear(modelname)
    mesh(modelname)
    job(modelname)
    G23 = results_G23(modelname)
    return G23

def YZShear(r_n,num):
    modelname = 'YZShear'+str(num)   
    parts_creation(modelname,r_n)
    assembly_creation(modelname)
    material_creation(modelname)
    material_assign(modelname)
    steps(modelname)
    facesets(modelname)
    loads_YZShear(modelname)
    mesh(modelname)
    job(modelname)
    G12= results_G12(modelname)
    return G12

def ZXShear(r_n,num):
    modelname = 'ZXShear'+str(num)   
    parts_creation(modelname,r_n)
    assembly_creation(modelname)
    material_creation(modelname)
    material_assign(modelname)
    steps(modelname)
    facesets(modelname)
    loads_ZXShear(modelname)
    mesh(modelname)
    job(modelname)
    G13 = results_G13(modelname)
    return G13

def RVE(num):
    r_n =[]
    for i in range(0,n):
        for j in range(0,m):
            k = random.uniform(r1, r2)
            r_n.append(k)
    mu = np.mean(r_n)
    sigma = np.std(r_n)
    N = m*n
    E11,v12,v13 = ZTens(r_n,num)
    E22,v21,v23 = YTens(r_n,num)
    E33,v31,v32 = XTens(r_n,num)
    G12 = YZShear(r_n,num)
    G13 = ZXShear(r_n,num)
    G23 = XYShear(r_n,num)
    row = [num,mu,sigma,N,E11,E22,E33,G12,G13,G23,v12,v13,v21,v23,v31,v32]
    with open('data.csv', 'a') as f:
        writer = csv.writer(f)
        writer.writerow(row)
        f.close()

def results():
    print(results_E11('ZTens'))
    print(results_E22('YTens'))
    print(results_E33('XTens'))
    print(results_G23('XYShear'))
    print(results_G13('ZXShear'))
    print(results_G12('YZShear'))

for i in range(0,num):
    RVE(i)
# results()

end_time = time.time()
execution_time = end_time - start_time
print("Time taken to run the code:", execution_time, "seconds")