'''
The visualization module automates some visualization procedures for
post processing data.

An example of using this is as follows:

.. code-block:: python

    import files
    import visualization as vis

    vis.ovito_xyz_to_gif(files.read_xyz("CNH_HCN.xyz"), "/fs/home/hch54/tmp", renderer='Tachyon')

- :func:`ovito_xyz_to_image`
- :func:`ovito_xyz_to_gif`

------------

'''
import os
import subprocess
from squid.files import write_xyz
from squid.files.misc import which, close_pipes


def get_ovito_obj(version="2.9.0"):
    '''
    This function returns the ovito object.  Note, currently the code below
    only works on version 2.9.0.
    '''
    ovito_path = which("ovitos")

    # Determine version
    ovito_pipe = subprocess.Popen(
        [ovito_path, "-v"], shell=False,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout = str(ovito_pipe.stdout.read().decode("utf-8").strip())

    assert "Ovito" in stdout,\
        "Error - Unable to access Ovito.  Please ensure it is in your PATH \
environment variable!"
    assert version in stdout,\
        "Error - Incorrect Ovito version!  It should be %s, but is %s."\
        % (version, stdout.strip().split()[-1])

    close_pipes(ovito_pipe)

    return ovito_path


def ovito_xyz_to_image(
        xyz, scratch,
        fname="image",
        camera_pos=(10, 0, 0), camera_dir=(-1, 0, 0),
        size=(800, 600),
        renderer="OpenGLRenderer",
        display_cell=False,
        renderer_settings={}):
    '''
    This function will, using the ovito python api, generate a png image of an
    xyz file.

    **Parameters**

        xyz: *str*
            A path to an xyz file.
        scratch: *str*
            A directory you want to have each image saved to.
        fname: *str, optional*
            The prefix for the image names.
        camera_pos: *tuple, float, optional*
            A tuple of x, y, and z coordinates for the camera to be positioned.
        camera_dir: *tuple, float, optional*
            The direction the camera is facing.
        size: *tuple, int, optional*
            Image size (width, height).
        delay: *int, optional*
            In the event of a gif, how long it should play for.
        renderer: *str, optional*
            What kind of renderer you wish to use: OpenGL or Tachyon.
        display_cell: *bool, optional*
            Whether to display the box around the system or not.
        renderer_settings: *dict, optional*
            Here you can change specific renderer settings.

    **Returns**

        None
    '''

    ovitos_path = get_ovito_obj("2.9.0")
    assert ovitos_path is not None,\
        "Error - Cannot find ovitos in the PATH env var."

    # Default settings for the renderers
    tach_rend_set = {
        'ambient_occlusion': False,
        'ambient_occlusion_brightness': 0.8,
        'ambient_occlusion_samples': 12,
        'antialiasing': True,
        'antialiasing_samples': 12,
        'direct_light': True,
        'direct_light_intensity': 0.9,
        'shadows': True
    }

    open_rend_set = {
        'antialiasing_level': 3
    }

    # General script to visualize things
    script = '''
from ovito.vis import *
from ovito.io import import_file

vp = Viewport()
vp.type = Viewport.Type.PERSPECTIVE
vp.camera_pos = $CAMERA_POS
vp.camera_dir = $CAMERA_DIR

settings = $SETTINGS
rs = RenderSettings(size=$SIZE, filename="$SCRATCH$FNAME.png", renderer=$RENDERER(**settings))

node = import_file("$XYZ", columns=["Particle Type", "Position.X", "Position.Y", "Position.Z"])
node.add_to_scene()

cell = node.source.cell
cell.display.enabled = $DISPLAY_CELL

vp.render(rs)
'''
    if not scratch.endswith("/"):
        scratch += "/"

    # Setup our renderer
    renderer = renderer.lower()
    if renderer.startswith("tach"):
        renderer = "TachyonRenderer"
        settings = tach_rend_set
        settings.update(renderer_settings)
    else:
        renderer = "OpenGLRenderer"
        settings = open_rend_set
        settings.update(renderer_settings)

    # Set all our variables
    holders = [
        ("$XYZ", str(xyz)),
        ("$CAMERA_POS", str(camera_pos)),
        ("$CAMERA_DIR", str(camera_dir)),
        ("$SIZE", str(size)),
        ("$SCRATCH", str(scratch)),
        ("$FNAME", str(fname)),
        ("$SETTINGS", str(settings)),
        ("$RENDERER", renderer),
        ("$DISPLAY_CELL", str(display_cell))
    ]
    for s_id, val in holders:
        script = script.replace(s_id, val)

    # Save the python script
    fptr = open("tmp_imgGen.py", 'w')
    fptr.write(script)
    fptr.close()

    # Call the ovitos python interpreter
    os.system("%s %s" % (ovitos_path, "tmp_imgGen.py"))
    os.system("rm tmp_imgGen.py")


def ovito_xyz_to_gif(
        frames, scratch,
        fname="image",
        camera_pos=(10, 0, 0), camera_dir=(-1, 0, 0),
        size=(800, 600),
        delay=10,
        display_cell=False,
        renderer="OpenGLRenderer",
        renderer_settings={},
        overwrite=False):
    '''
    This function will, using the ovito python api, generate either a single
    image or a gif of the input frames.  Note, a gif is only generated when
    more than one frame exists.

    **Parameters**

        frames: *str* or *list,* :class:`squid.structures.atom.Atom`
            A list of frames you wish to generate an image for, or a path to
            an xyz file.
        scratch: *str*
            A directory you want to have each image saved to.
        fname: *str, optional*
            The prefix for the image names.
        camera_pos: *tuple, float, optional*
            A tuple of x, y, and z coordinates for the camera to be positioned.
        camera_dir: *tuple, float, optional*
            The direction the camera is facing.
        size: *tuple, int, optional*
            Image size (width, height).
        delay: *int, optional*
            In the event of a gif, how long it should play for.
        display_cell: *bool, optional*
            Whether to display the box around the system or not.
        renderer: *str, optional*
            What kind of renderer you wish to use: OpenGL or Tachyon.
        renderer_settings: *dict, optional*
            Here you can change specific renderer settings.
        overwrite: *bool, optional*
            Whether to delete any files already existing in the scratch dir.

    **Returns**

        None
    '''

    convert_path = which("convert")
    assert convert_path is not None,\
        "Error - Cannot find convert in the PATH env var."

    if fname.endswith(".gif"):
        fname.replace(".gif", "")

    # First ensure we have frames and things in the correct format
    if isinstance(frames, str):
        frames = open(frames)
    if not isinstance(frames[0], list):
        frames = [frames]
    if not scratch.endswith("/"):
        scratch += "/"

    # Next, ensure scratch exists
    if not os.path.exists(scratch):
        os.system("mkdir -p %s" % scratch)
    elif len(os.listdir(scratch)) > 0:
        if overwrite:
            os.system("rm %s/*.png" % scratch)
        else:
            raise Exception("Error - Scratch directory is not empty!")

    # For each frame, generate an image
    for i, frame in enumerate(frames):
        write_xyz(frame, "tmp.xyz")
        ovito_xyz_to_image(
            "tmp.xyz", scratch,
            fname="%04d" % i,
            camera_pos=camera_pos,
            camera_dir=camera_dir,
            size=size,
            renderer=renderer,
            renderer_settings=renderer_settings,
            display_cell=display_cell
        )
        os.system("rm tmp.xyz")

    # If more than one frame exists, compile to gif
    if len(frames) > 1:
        cmd = "convert -delay $DELAY -loop 0 $(ls -v $PATH/*.png) output.gif"
        holders = [
            ("$PATH", str(scratch)),
            ("$DELAY", str(delay))
        ]
        for s_id, val in holders:
            cmd = cmd.replace(s_id, val)
        os.system(cmd)
        os.rename("output.gif", fname + ".gif")
