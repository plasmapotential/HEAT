
import CADClass

def run(obj):
    CAD = CADClass.CAD()
    CAD.STPfile = obj.STPfilebox.text()
    CAD.STLpath = obj.STLfilebox.text()
    CAD.FreeCADPath = obj.FCfilebox.text()
    CAD.ROIGridRes = obj.ROIGridbox.text()
    CAD.gridRes = obj.gridbox.text()
    CAD.permute_mask = obj.permutebox.isChecked()
    CAD.unitConvert = float(obj.unitbox.text())
    CAD.assembly_mask = obj.assybox.isChecked()
    CAD.rootDir = obj.rootfilebox.text()
    #self.CAD = Import.open(self.STPfile)
    #self.CADdoc = FreeCAD.ActiveDocument
