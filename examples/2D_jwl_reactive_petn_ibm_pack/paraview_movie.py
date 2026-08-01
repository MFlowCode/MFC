# ParaView batch renderer for 2D_jwl_reactive_petn_ibm_pack.
#
# Renders the PETN detonation and the launched IB pack as a numerical schlieren
# (|grad p|, log inferno on black) with the four moving cylinders overlaid, then
# encodes an mp4. Schlieren, not raw pressure: the blast decays toward ambient over
# the run, so any pressure map goes black in the flight phase. |grad p| keeps every
# wave front glowing, the CJ detonation AND each sphere's bow shock and wake as it
# flies, so the physics stays visible the whole time.
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
OUT = os.path.join(HERE, "petn_pack_reactive_pv.mp4")
# |grad p| range (Pa/m), log-scaled. Tuned so flight-phase wakes/bow shocks
# (~1e9-1e10) sit in the bright half of inferno and the CJ front (~1e12) saturates.
GMIN, GMAX = 1.0e7, 5.0e10
FPS = 5  # slow playback so the detonation and the flight are both readable

os.makedirs(FRAMES, exist_ok=True)

# --- flow field: rectilinear grid, schlieren = |grad(pres)| ---
field = OpenDataFile(SERIES)
field.MeshStatus = ["rectilinear_grid"]
field.CellArrayStatus = ["pres"]

grad = Gradient(Input=field)
grad.ScalarArray = ["CELLS", "pres"]
grad.ResultArrayName = "pgrad"

grad.UpdatePipeline()
bx0, bx1, by0, by1, _, _ = grad.GetDataInformation().GetBounds()
cx, cy = 0.5 * (bx0 + bx1), 0.5 * (by0 + by1)
half_y = 0.5 * (by1 - by0)
aspect = (bx1 - bx0) / (by1 - by0)
vh = 760
vw = int(vh * aspect) + 200  # + room for the colour bar

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [vw, vh]
view.OrientationAxesVisibility = 0
view.CameraParallelProjection = 1
view.UseColorPaletteForBackground = 0
view.BackgroundColorMode = "Single Color"
view.Background = [0, 0, 0]  # black: the bright blast and blue pack pop against it

disp = Show(grad, view)
ColorBy(disp, ("CELLS", "pgrad"))  # coloured by vector magnitude
pres_lut = GetColorTransferFunction("pgrad")
pres_lut.ApplyPreset("Inferno (matplotlib)", True)
pres_lut.UseLogScale = 1
pres_lut.RescaleTransferFunction(GMIN, GMAX)
disp.SetScalarBarVisibility(view, True)
bar = GetScalarBar(pres_lut, view)
bar.Title = "|grad p| (Pa/m)"
bar.ComponentTitle = ""
bar.TitleColor = [1, 1, 1]
bar.LabelColor = [1, 1, 1]

# --- domain outline so the launched pack reads as leaving the box, not floating ---
frame = Outline(Input=field)
fdisp = Show(frame, view)
fdisp.ColorArrayName = ["POINTS", ""]
fdisp.DiffuseColor = [0.5, 0.5, 0.5]
fdisp.LineWidth = 1.5

# --- immersed pack: ib_bodies point mesh drawn as smooth sphere glyphs ---
# The gridded ib_markers field does NOT track moving bodies in the post output
# (it collapses to a fixed artifact for t>0), but the ib_bodies point mesh carries
# the correct moving centroids. Sphere glyphs (high theta/phi resolution) give a
# clean round particle instead of the low-facet 2D-circle octagon.
bodies = OpenDataFile(SERIES)
bodies.MeshStatus = ["ib_bodies"]
bodies.PointArrayStatus = ["ib_diameter"]
glyph = Glyph(Input=bodies, GlyphType="Sphere")
glyph.GlyphType.Radius = 0.5           # unit diameter; scaled by ib_diameter below
glyph.GlyphType.ThetaResolution = 48
glyph.GlyphType.PhiResolution = 48
glyph.ScaleArray = ["POINTS", "ib_diameter"]
glyph.ScaleFactor = 1.0
glyph.GlyphMode = "All Points"
bdisp = Show(glyph, view)
bdisp.ColorArrayName = ["POINTS", ""]  # solid colour, not scalar-mapped
bdisp.DiffuseColor = [0.35, 0.62, 0.95]  # bright steel blue, pops on black
bdisp.Representation = "Surface"

# --- 2D top-down camera ---
view.InteractionMode = "2D"
view.CameraPosition = [cx, cy, 0.5]
view.CameraFocalPoint = [cx, cy, 0.0]
view.CameraViewUp = [0, 1, 0]
view.CameraParallelScale = half_y

# Time annotation. The silo series "time" is the save index (0,1,2,...); physical
# time = index * t_save = index * 2.0e-6 s = index * 2.0 microseconds.
tann = AnnotateTimeFilter(field)
tann.Scale = 2.0
tann.Format = "t = {time:.1f} us"
tdisp = Show(tann, view)
tdisp.Color = [1, 1, 1]
tdisp.FontSize = 18
tdisp.WindowLocation = "Upper Center"

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
times = GetTimeKeeper().TimestepValues
times = list(times) if times is not None else [0.0]

for i, t in enumerate(times):
    view.ViewTime = t
    UpdatePipeline(t, field)
    UpdatePipeline(t, glyph)
    Render(view)
    SaveScreenshot(os.path.join(FRAMES, f"frame_{i:04d}.png"), view,
                   ImageResolution=[vw, vh], TransparentBackground=0)

print("rendered", len(times), "frames")
subprocess.run([
    "ffmpeg", "-y", "-framerate", str(FPS),
    "-i", os.path.join(FRAMES, "frame_%04d.png"),
    "-c:v", "libx264", "-pix_fmt", "yuv420p",
    "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", OUT,
], check=True)
print("wrote", OUT)
