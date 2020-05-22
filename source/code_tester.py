#Code tester for HEAT
import CADClass
import MHDClass
import toolsClass
import heatfluxClass
import time

t0 = time.time()

machine = 'nstx'
#For NSTXU
if machine == 'nstx':
    infile = 'NSTXU_input.csv'
    PartsFile = 'NSTXUparts.csv'
    IntersectFile = 'NSTXUintersects.csv'
    #t=1004
    t=851
#For ST40
if machine == 'st40':
    infile = 'ST40_input.csv'
    PartsFile = 'ST40parts.csv'
    IntersectFile = 'ST40intersects.csv'
    t=1

tools = toolsClass.tools()

HFpointcd = True
print("-"*70)
print("CAD MODULE INITIALIZED")
CAD = CADClass.CAD()

#Strip STP file if necessary
if False:
    print("Stripping STP file")
    rawSTP = '/u/tlooby/ST40/CAD/20190904_IVC2_divertor.stp'
    outSTP = '/u/tlooby/ST40/CAD/IVC2_PFCs_TL.stp'
    CAD.stripSTPfile(PartsFile,rawSTP,outfile=outSTP)
    input()


#Import CAD
if True:
    print("Import CAD and generate ROI meshes")
    import numpy as np
    tools.initializeInput(CAD, mode='CLI', infile=infile)
    CAD.loadSTEP()
    CAD.getROIfromfile(PartsFile)
    CAD.getROImeshes()
    CAD.writeMesh2file(CAD.ROImeshes, CAD.ROI, path=CAD.STLpath)
    print("CAD Completed.  Time Elapsed: {:f}".format(time.time() - t0))

#Import MHD EFIT
if True:
    print('\n')
    print("-"*70)
    print("MHD MODULE INITIALIZED")
    MHD = MHDClass.MHD()
    tools.initializeInput(MHD, mode='CLI', infile=infile)

    #For NSTXU
    if machine=='nstx':
        gfile = '/home/tlooby/results/baseline_mafot_pwrblnce/g116313.00851'
        MHD.get_mhd_inputs('nstx',gfile)

    #For ST40
    if machine=='st40':
        gfile = '/u/tlooby/ST40/g000002.00001'
        MHD.get_mhd_inputs('st40',gfile)

    try:
        MHD.make1EFITobject(MHD.shot, MHD.timesteps[0])
    except:
        MHD.make1EFITobject(MHD.shot, MHD.timesteps)


    NCPUs = 4
    controlfile = '_lamCTL.dat'
    controlfilePath = MHD.dataPath + '/' + '{:06d}/'.format(t)
    gridfile = MHD.dataPath + '/' + '{:06d}/grid.dat'.format(t)
    outputFile = controlfilePath + 'lam.dat'
    print("MHD Load Completed.  Time Elapsed: {:f}".format(time.time() - t0))

#Calculate Bfield point cloud at tile surface
if True:
    ctrs = CAD.ROIctrs[0]/1000.0
    R,Z,phi = tools.xyz2cyl(ctrs[:,0],ctrs[:,1],ctrs[:,2])
    Bxyz = MHD.Bfield_pointcloud(MHD.ep, R, Z, phi)
    MHD.write_B_pointcloud(ctrs,Bxyz,controlfilePath)
    print("Bfield Completed.  Time Elapsed: {:f}".format(time.time() - t0))

#Calculate Bfield in other places
if False:
    ctrs = CAD.ROIctrs[0]
    R,Z,phi = tools.xyz2cyl(ctrs[:,0],ctrs[:,1],ctrs[:,2])
    #R = np.linspace(np.amin(R),np.amax(R),1000)
    R = np.linspace(750.0,900.0,1000)
    Z = np.ones(R.shape)*(0.0)
    phi = np.average(phi) * np.ones(R.shape) + 180.0
    x,y,z = tools.cyl2xyz(R,Z,phi)
    ctrs = np.concatenate((x,y,z)).reshape((-1, 3), order='F')
    Bxyz = MHD.Bfield_pointcloud(MHD.ep, R, Z, phi)
    MHD.write_B_pointcloud(ctrs,Bxyz,controlfilePath)
    print("Bfield Completed.  Time Elapsed: {:f}".format(time.time() - t0))
    #test point
#    x = 360.0
#    y = 533.0
#    z = -871.0
#    r,z,phi = tools.xyz2cyl(x,y,z)
#    print(r)
#    print(z)
#    print(phi)
#    Bxyz = MHD.Bfield_pointcloud(MHD.ep, r, z, phi)
#    print(Bxyz)


#Mesh potential intersections by name
if HFpointcd:
    print('\n')
    print("-"*70)
    print("MESH MODULE INITIALIZED")
    intersects = CAD.readIntersects(IntersectFile, CAD.CADparts)
    print("Potential intersects on these tiles:")
    for i in range(len(CAD.CADparts)):
        if intersects[i] == 1.0:
            print(CAD.CADparts[i].Label)
    res=CAD.gridRes
    targetMeshes, targetLabels = CAD.meshPotentialIntersects(intersects, CAD.ROIparts,CAD.CADparts, resolution=res)
    CAD.writeMesh2file(targetMeshes, targetLabels, path=CAD.STLpath,resolution=res)
    print("Mesh Generation Completed.  Time Elapsed: {:f}".format(time.time() - t0))

#Mesh Potential Intersections using internal scripts / tests
if False:
    print('\n')
    print("-"*70)
    print("MESH MODULE INITIALIZED")
    deltaPhis = CAD.findRelativePhi(CAD.ROIparts, CAD.CADparts)
    COMs = CAD.getCOMs(CAD.ROIparts)
    R,Z,phi = tools.xyz2cyl(COMs[:,0],COMs[:,1],COMs[:,2])
    sourceBts = MHD.ep.BtFunc.ev(COMs[:,0]/1000.0,COMs[:,1]/1000.0)
    intersects = CAD.findPotentialIntersectParts(deltaPhis,sourceBts,CAD.ROIparts,CAD.CADparts)
    print("Potential intersects on these tiles:")
    for i in range(len(CAD.CADparts)):
        if intersects[0,i] == 1.0:
            print(CAD.CADparts[i].Label)
    res=CAD.gridRes
    targetMeshes, targetLabels = CAD.meshPotentialIntersects(intersects, CAD.ROIparts,CAD.CADparts, resolution=res)
    CAD.writeMesh2file(targetMeshes, targetLabels, path=CAD.STLpath,resolution=res)
    #for i in range(len(deltaPhis.flatten())):
    #    print("{:f}\t{:f}\t".format(deltaPhis.flatten()[i],intersects.flatten()[i])+CAD.getLabels(CAD.CADparts)[i])
    print("Mesh Generation Completed.  Time Elapsed: {:f}".format(time.time() - t0))

#Trace Bfield
if False:
    #COMs = CAD.getCOMs(CAD.ROIparts)
    #x = COMs[:,0][0]
    #y = COMs[:,1][0]
    #z = COMs[:,2][0]
    x = 360.0
    y = 533.0
    z = -870.0
    #x = x[-1]
    #y = y[-1]
    #z = z[-1]
    xyz = np.array([x,y,z])
    controlfile = '_structCTL.dat'
    t = 1
    dphi = 1.0
    MHD.ittStruct = 360
    controlfilePath = MHD.dataPath + '/' + '{:06d}/'.format(t)
    gridfile = MHD.dataPath + '/' + '{:06d}/struct_grid.dat'.format(t)
    MHD.writeControlFile(controlfile, t, mode='struct')
    MHD.writeMAFOTpointfile(xyz,gridfile)
    MHD.getFieldpath(dphi, gridfile, controlfilePath, controlfile, paraview_mask=True)
    print("Traced B field")


#Initialize HF Object
#if True:
if HFpointcd:
    print('\n')
    print("-"*70)
    print("HEAT FLUX MODULE INITIALIZED")
    HF = heatfluxClass.heatFlux()
    tools.initializeInput(HF, mode='CLI', infile=infile)
    HF.makePFCs(MHD, CAD.ROImeshes, CAD.ROIctrs, CAD.ROInorms, CAD.ROIareas)
    PFC = HF.PFCs[0]
    HF.write_shadow_pointcloud(PFC.centers,PFC.shadowed_mask,controlfilePath)



#Check for intersections with MAFOT struct
if HFpointcd:
    print('\n')
    print("-"*70)
    print("INTERSECTION CHECK INITIALIZED")
    HF.findShadows_structure(MHD,PFC,targetMeshes)
    print("Intersection Test Completed.  Time Elapsed: {:f}".format(time.time() - t0))

#    gridfile = MHD.dataPath + '/' + '{:06d}/struct_grid.dat'.format(t)
#    MHD.writeMAFOTpointfile(PFC.centers,gridfile)



#Run MAFOT Laminar
if HFpointcd:
    print('\n')
    print("-"*70)
    print("MAFOT LAMINAR MODULE INITIALIZED")

    MHD.writeControlFile(controlfile, t, mode='laminar')
    use = np.where(PFC.shadowed_mask != 1)[0]
    MHD.writeMAFOTpointfile(PFC.centers[use],gridfile)
    MHD.runMAFOTlaminar(gridfile,controlfilePath,controlfile,NCPUs)


    #unique, counts = np.unique(PFC.psimin, return_counts=True)
    #print(unique[-1])
    #print(counts[-1])
    #print(np.sort(PFC.psimin))


#Create Heat Flux Profile
if HFpointcd:
    print('\n')
    print("-"*70)
    print("HEAT FLUX CALCULATION")
    P = 3.0 # [MW]
    lq = 3.0 # [mm]
    S = 3.0 # [mm]
    qBG = 0.0 #[MW/m^2]
    HF.readMAFOTLaminarOutput(PFC,outputFile)
    q = HF.getHFprofile(PFC, P, lq, S, qBG)
    qDiv = HF.q_div(PFC, MHD, q)
    print('Maximum heat load on tile: {:f}'.format(max(qDiv)))


#Print Stuff
if HFpointcd:
    print('Input Power = {:f}'.format(P))
    #print('Analytically Integrated Power = {:f}'.format(P_old))
    print('Tessellated Total Power = {:f}'.format(HF.power_sum_mesh(PFC)))





#Create pointclouds for paraview
if HFpointcd:
    print("Paraview Point Clouds...")
    R,Z,phi = tools.xyz2cyl(PFC.centers[:,0],PFC.centers[:,1],PFC.centers[:,2])
    Bxyz = MHD.Bfield_pointcloud(MHD.ep, R, Z, phi)
    MHD.write_B_pointcloud(PFC.centers,Bxyz,controlfilePath)
    HF.write_shadow_pointcloud(PFC.centers,PFC.shadowed_mask,controlfilePath)
    #CAD.write_normal_pointcloud(PFC.centers,PFC.norms,controlfilePath)
    HF.write_heatflux_pointcloud(PFC.centers,qDiv,controlfilePath)
    #structOutfile = MHD.dataPath + '/' + '{:06d}/struct.csv'.format(PFC.t)
    #HF.PointCloudfromStructOutput(structOutfile)

print("All Completed.  Time Elapsed: {:f}".format(time.time() - t0))
print('Done')
