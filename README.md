# cimp
Simple command line program for editing images in the vein of [GIMP](https://www.gimp.org/)

## About
This started a project for one my classes in college, along with my classmates Alex Ma and Amy Wang. I have since added much more functionality and refactored the previously present image processing functions to be more efficient. The actual driver for the program remains largely the same.

## Building
Currently this only runs on linux systems, but any major linux distribution will work
To run on windows, install [wsl](https://learn.microsoft.com/en-us/windows/wsl/) in powershell or cmd:
```
wsl --install
```

This will install Ubuntu by default, to specify a distro run
```
wsl --install -d <DistroName>
```

To build, first clone the repo:
```
git clone https://github.com/ansht2000/cimp.git
```

Then cd into the created directory, it should be called cimp, and simply run:
```
make cimp
```

This will create the cimp executable, which you can run with ./cimp <input_file> <output_file> <operation> \[<args>\]

## Commands
Currently there are 10 commands available:
grayscale
binarize
crop
flip
rotate
gradient
seam
blend
pointilism
dither
