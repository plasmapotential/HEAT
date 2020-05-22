import sys
path = '/opt/paraview/ParaView-5.7.0-MPI-Linux-Python3.7-64bit/lib/python3.7/site-packages'
sys.path.append(path)
from paraview import simple

pc = simple.CSVReader(FileName=['/home/tlooby/source/HEAT/rev4/data/204118/001004/B_pointcloud.csv'])
t2p = simple.TableToPoints(Input=pc)
t2p.XColumn = '# X'
t2p.YColumn = 'Y'
t2p.ZColumn = 'Z'
writer = simple.CreateWriter("/u/tlooby/source/HEAT/rev4/data/204118/001004/B_pointcloud.vtk", t2p)
writer.UpdatePipeline()
