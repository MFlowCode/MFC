# state file generated using paraview version 5.11.0
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
import os, re
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [2328, 1162]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [0.5, 0.5, 0.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [0.5, 0.5, 2.7320508075688776]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.7071067811865476
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(2328, 1162)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'VisItSiloReader'
case_dir=os.getcwd()
print(case_dir)
# create a new 'VisItSiloReader'
files = glob.glob(f"{case_dir}/silo_hdf5/root/*")
print(files)
sorted_files = sorted(files, key=lambda x: int(x.rsplit('_', 1)[-1][:-5]))
print(sorted_files)
collection_0silo = VisItSiloReader(registrationName='collection_0.silo*', FileName=sorted_files)
collection_0silo.MeshStatus = ['rectilinear_grid']
collection_0silo.CellArrayStatus = ['alpha1', 'alpha2', 'alpha_rho1', 'alpha_rho2', 'pres', 'vel1', 'vel2']

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from collection_0silo
collection_0siloDisplay = Show(collection_0silo, renderView1, 'UniformGridRepresentation')

# get 2D transfer function for 'vtkBlockColors'
vtkBlockColorsTF2D = GetTransferFunction2D('vtkBlockColors')

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')
vtkBlockColorsLUT.InterpretValuesAsCategories = 1
vtkBlockColorsLUT.AnnotationsInitialized = 1
vtkBlockColorsLUT.TransferFunction2D = vtkBlockColorsTF2D
vtkBlockColorsLUT.Annotations = ['0', '0', '1', '1', '2', '2', '3', '3', '4', '4', '5', '5', '6', '6', '7', '7', '8', '8', '9', '9', '10', '10', '11', '11']
vtkBlockColorsLUT.ActiveAnnotatedValues = ['0', '1', '2', '3']
vtkBlockColorsLUT.IndexedColors = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.63, 0.63, 1.0, 0.67, 0.5, 0.33, 1.0, 0.5, 0.75, 0.53, 0.35, 0.7, 1.0, 0.75, 0.5]

# get opacity transfer function/opacity map for 'vtkBlockColors'
vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')

# trace defaults for the display properties.
collection_0siloDisplay.Representation = 'Surface'
collection_0siloDisplay.ColorArrayName = ['FIELD', 'vtkBlockColors']
collection_0siloDisplay.LookupTable = vtkBlockColorsLUT
collection_0siloDisplay.SelectTCoordArray = 'None'
collection_0siloDisplay.SelectNormalArray = 'None'
collection_0siloDisplay.SelectTangentArray = 'None'
collection_0siloDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
collection_0siloDisplay.SelectOrientationVectors = 'None'
collection_0siloDisplay.ScaleFactor = 0.1
collection_0siloDisplay.SelectScaleArray = 'None'
collection_0siloDisplay.GlyphType = 'Arrow'
collection_0siloDisplay.GlyphTableIndexArray = 'None'
collection_0siloDisplay.GaussianRadius = 0.005
collection_0siloDisplay.SetScaleArray = [None, '']
collection_0siloDisplay.ScaleTransferFunction = 'PiecewiseFunction'
collection_0siloDisplay.OpacityArray = [None, '']
collection_0siloDisplay.OpacityTransferFunction = 'PiecewiseFunction'
collection_0siloDisplay.DataAxesGrid = 'GridAxesRepresentation'
collection_0siloDisplay.PolarAxes = 'PolarAxesRepresentation'
collection_0siloDisplay.ScalarOpacityUnitDistance = 0.06564197879454707
collection_0siloDisplay.ScalarOpacityFunction = vtkBlockColorsPWF
collection_0siloDisplay.TransferFunction2D = vtkBlockColorsTF2D
collection_0siloDisplay.OpacityArrayName = ['CELLS', 'alpha1']
collection_0siloDisplay.ColorArray2Name = ['CELLS', 'alpha1']
collection_0siloDisplay.SliceFunction = 'Plane'
collection_0siloDisplay.SelectInputVectors = [None, '']
collection_0siloDisplay.WriteLog = ''

# init the 'Plane' selected for 'SliceFunction'
collection_0siloDisplay.SliceFunction.Origin = [0.5, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for vtkBlockColorsLUT in view renderView1
vtkBlockColorsLUTColorBar = GetScalarBar(vtkBlockColorsLUT, renderView1)
vtkBlockColorsLUTColorBar.Title = 'vtkBlockColors'
vtkBlockColorsLUTColorBar.ComponentTitle = ''

# set color bar visibility
vtkBlockColorsLUTColorBar.Visibility = 1

# show color legend
collection_0siloDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# restore active source
SetActiveSource(collection_0silo)
# ----------------------------------------------------------------

# Ensure all time steps are considered
timeKeeper = GetTimeKeeper()
timeSteps = timeKeeper.TimestepValues

animationScene = GetAnimationScene()

directory_path = f"{case_dir}/render"
os.makedirs(directory_path, exist_ok=True)

i = 0
# Save all timesteps
for t in timeSteps:
    animationScene.AnimationTime = t
    SaveScreenshot(f"{case_dir}/render/pic.{i:04d}.png", renderView1, ImageResolution=[2110,722])
    print(i)
    i = i + 1

