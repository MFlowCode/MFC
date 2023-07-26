import math, argparse

from PIL import Image, ImageFilter

parser = argparse.ArgumentParser(prog='img2stl', description='Turn an image into an STL file')
parser.add_argument('input',  type=str, help='Path to the image file')
parser.add_argument('output', type=str, help='Path to the output OBJ file')

VARS = vars(parser.parse_args())

img = Image.open(VARS['input'])

W, H = img.size

pixels = img.load()
for y in range(H):
    for x in range(W):
        if pixels[x,y][3] < 255:
            pixels[x,y] = (0,0,0,255)

img.save('temp.png')
img = img.convert('L')

THRESHOLD = 0.1

with open(VARS['output'], 'w') as f:
    i = 0
    pixels = img.load()
    for y in range(H):
        for x in range(W):
            c = pixels[x,y]

            if c < 255 * THRESHOLD: continue

            area = c / 255.0
            d    = (1.0 - math.sqrt(area)) / 2.0

            def tx(_x: int, _f: bool):
                if _f: _x = _x + d*(1 - _x*2)
                return -1.0+2.0 * (x + _x)/float(W)

            def ty(_y: int, _f: bool):
                if _f: _y = _y + d*(1 - _y*2)
                return  1.0-2.0 * ((1-H/float(W))/2 + (y+_y)/float(W))

            def tz(_z: int, _f: bool):
                if _f: _z = _z + d*(1 - _z*2)
                return max(1.0 / float(W), 1.0 / float(H)) * (-1.0+2.0 * _z)

            def p1(_x: int, _y: int, _z: int):
                return tx(_x, True), ty(_y, True), tz(_z, False)

            for _x, _y, _z in [p1(0,0,0), p1(1,0,0), p1(0,1,0), p1(1,1,0),
                               p1(0,0,1), p1(1,0,1), p1(0,1,1), p1(1,1,1)]:
                f.write(f'v {_x} {_y} {_z}\n')
            
            f.write(f'f {i+1} {i+2} {i+3}\n'); f.write(f'f {i+2} {i+3} {i+4}\n')
            f.write(f'f {i+5} {i+6} {i+7}\n'); f.write(f'f {i+6} {i+7} {i+8}\n')

            i = i + 8

            def p2(_x: int, _y: int, _z: int):
                return tx(_x, False), ty(_y, True), tz(_z, True)

            for _x, _y, _z in [p2(0,0,0), p2(1,0,0), p2(0,1,0), p2(1,1,0),
                               p2(0,0,1), p2(1,0,1), p2(0,1,1), p2(1,1,1)]:
                f.write(f'v {_x} {_y} {_z}\n')

            f.write(f'f {i+1} {i+5} {i+7}\n'); f.write(f'f {i+1} {i+3} {i+7}\n')
            f.write(f'f {i+2} {i+6} {i+8}\n'); f.write(f'f {i+2} {i+4} {i+8}\n')

            i = i + 8

            def p3(_x: int, _y: int, _z: int):
                return tx(_x, True), ty(_y, False), tz(_z, True)

            for _x, _y, _z in [p3(0,0,0), p3(1,0,0), p3(0,1,0), p3(1,1,0),
                               p3(0,0,1), p3(1,0,1), p3(0,1,1), p3(1,1,1)]:
                f.write(f'v {_x} {_y} {_z}\n')

            f.write(f'f {i+1} {i+2} {i+5}\n'); f.write(f'f {i+2} {i+5} {i+6}\n')
            f.write(f'f {i+3} {i+4} {i+7}\n'); f.write(f'f {i+4} {i+7} {i+8}\n')

            i = i + 8
