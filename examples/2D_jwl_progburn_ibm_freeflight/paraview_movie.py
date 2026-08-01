# ParaView batch renderer for 2D_jwl_progburn_ibm_freeflight.
#
# Renders the pressure field (log-scaled, inferno) over the whole run with the
# moving immersed cylinder overlaid, then encodes an mp4. Reads this case's own
# silo_hdf5/collection.silo.series.
#
#   /Applications/ParaView-5.11.0.app/Contents/bin/pvbatch paraview_movie.py
#
# IMPORTANT: post-process the run SERIALLY (`./mfc.sh run case.py -n 1 -t
# post_process`). MFC's parallel post_process (-n>1) blanks the upper MPI ranks
# in the silo output; the simulation itself may be run in parallel.
import os
import glob
import subprocess

from paraview.simple import *

HERE = os.path.dirname(os.path.abspath(__file__))
SERIES = os.path.join(HERE, "silo_hdf5", "collection.silo.series")
FRAMES = os.path.join(HERE, "pv_frames")
OUT = os.path.join(HERE, "jwl_progburn_freeflight_pv.mp4")
PMIN, PMAX = 2.0e5, 2.0e8  # fixed log range so frames are comparable

os.makedirs(FRAMES, exist_ok=True)

# --- flow field: rectilinear grid coloured by pressure ---
field = OpenDataFile(SERIES)
field.MeshStatus = ["rectilinear_grid"]
field.CellArrayStatus = ["pres"]

field.UpdatePipeline()
bx0, bx1, by0, by1, _, _ = field.GetDataInformation().GetBounds()
cx, cy = 0.5 * (bx0 + bx1), 0.5 * (by0 + by1)
half_y = 0.5 * (by1 - by0)
aspect = (bx1 - bx0) / (by1 - by0)
vh = 760
vw = int(vh * aspect) + 180  # + room for the colour bar

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [vw, vh]
view.OrientationAxesVisibility = 0
view.CameraParallelProjection = 1
view.UseColorPaletteForBackground = 0
view.BackgroundColorMode = "Single Color"
view.Background = [1, 1, 1]

disp = Show(field, view)
ColorBy(disp, ("CELLS", "pres"))
pres_lut = GetColorTransferFunction("pres")
pres_lut.ApplyPreset("Inferno (matplotlib)", True)
pres_lut.UseLogScale = 1
pres_lut.RescaleTransferFunction(PMIN, PMAX)
disp.SetScalarBarVisibility(view, True)
bar = GetScalarBar(pres_lut, view)
bar.Title = "Pressure (Pa)"
bar.ComponentTitle = ""
bar.TitleColor = [0, 0, 0]
bar.LabelColor = [0, 0, 0]

# --- immersed body: the real IB footprint from the ib_markers field ---
# ib_markers = 1 in solid cells, 0 in fluid. Threshold gives the actual moving
# cylinder as a filled body (honest cell-resolution shape), drawn opaque on top
# of the pressure field so its interior does not read as a pressure void.
solid = OpenDataFile(SERIES)
solid.MeshStatus = ["rectilinear_grid"]
solid.CellArrayStatus = ["ib_markers"]
body = Threshold(Input=solid)
body.Scalars = ["CELLS", "ib_markers"]
body.LowerThreshold = 0.5
body.UpperThreshold = 1.0
bdisp = Show(body, view)
bdisp.ColorArrayName = ["CELLS", ""]  # solid colour, not scalar-mapped
bdisp.DiffuseColor = [0.82, 0.82, 0.85]  # neutral metallic grey
bdisp.Representation = "Surface"

# --- 2D top-down camera on the [-0.1, 0.1] domain ---
view.InteractionMode = "2D"
view.CameraPosition = [cx, cy, 0.5]
view.CameraFocalPoint = [cx, cy, 0.0]
view.CameraViewUp = [0, 1, 0]
view.CameraParallelScale = half_y

# Time annotation. The silo series "time" is the save index (0,1,2,...); physical
# time = index * t_step_save * dt = index * 20 * 2.5e-8 s = index * 0.5 microseconds.
tann = AnnotateTimeFilter(field)
tann.Scale = 0.5
tann.Format = "t = {time:.1f} us"
tdisp = Show(tann, view)
tdisp.Color = [0, 0, 0]
tdisp.FontSize = 16
tdisp.WindowLocation = "Upper Center"

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
times = GetTimeKeeper().TimestepValues
times = list(times) if times is not None else [0.0]

for i, t in enumerate(times):
    view.ViewTime = t
    UpdatePipeline(t, field)
    UpdatePipeline(t, body)
    Render(view)
    SaveScreenshot(os.path.join(FRAMES, f"frame_{i:04d}.png"), view,
                   ImageResolution=[860, 760], TransparentBackground=0)

print("rendered", len(times), "frames")
subprocess.run([
    "ffmpeg", "-y", "-framerate", "20",
    "-i", os.path.join(FRAMES, "frame_%04d.png"),
    "-c:v", "libx264", "-pix_fmt", "yuv420p",
    "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", OUT,
], check=True)
print("wrote", OUT)
